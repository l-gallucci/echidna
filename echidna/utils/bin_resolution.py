"""Quality-prioritized overlap resolution with optional cluster replacements.

All changes are exported to a new directory. Original genome FASTAs are read-only.
"""
from collections import Counter, defaultdict
import csv
import json
import math
import os
from pathlib import Path
import tempfile

from .bin_overlap import analyze_bins, discover_bins, read_bin, sequence_digest, write_reports, _write_tsv


def genome_stem(name):
    for suffix in ('.gzip', '.gz'):
        if name.lower().endswith(suffix):
            name = name[:-len(suffix)]
            break
    for suffix in ('.fasta', '.fna', '.fa'):
        if name.lower().endswith(suffix):
            return name[:-len(suffix)]
    return Path(name).stem


def read_table(path):
    with Path(path).open(encoding='utf-8-sig', newline='') as stream:
        header = stream.readline()
        stream.seek(0)
        reader = csv.DictReader(stream, delimiter='\t' if '\t' in header else ',')
        fields = [str(k).strip().lower() for k in (reader.fieldnames or [])]
        if not fields or len(fields) != len(set(fields)):
            raise ValueError(f'Missing or duplicate table headers in {path}')
        rows = []
        for row in reader:
            if None in row or any(v is None for v in row.values()):
                raise ValueError(f'Malformed table row in {path}')
            clean = {k.strip().lower(): v.strip() for k, v in row.items()}
            if any(clean.values()):
                rows.append(clean)
    return fields, rows


def genome_column(fields, path):
    for field in ('genome', 'name', 'bin'):
        if field in fields:
            return field
    raise ValueError(f'{path} needs a genome, Name or bin column')


class Inventory:
    def __init__(self, files, before, candidate_roots=None, recursive=False, extensions='fa,fna,fasta'):
        self.files = dict(files)
        self.primary = set(files)
        self.aliases = defaultdict(set)
        self.cache = {key: Counter() for key in files}
        for record in before['records']:
            self.cache[record.bin_id][record.digest] += 1
        for key, path in files.items():
            self.add_aliases(key, path, key)
        if candidate_roots:
            for label, path in discover_bins(candidate_roots, recursive, extensions).items():
                key = 'candidate:' + label
                while key in self.files:
                    key = 'candidate:' + key
                same_name = [k for k in self.primary if files[k].name == path.name]
                counts = None
                if len(same_name) == 1:
                    counts = self.read_counts(path)
                    if counts == self.cache[same_name[0]]:
                        # dRep winners are usually copies of the original candidates.
                        self.add_aliases(same_name[0], path, key)
                        continue
                self.files[key] = path
                if counts is not None:
                    self.cache[key] = counts
                self.add_aliases(key, path, key)

    @staticmethod
    def read_counts(path):
        return Counter(sequence_digest(seq) for _, seq in read_bin(path))

    def counts(self, key):
        if key not in self.cache:
            self.cache[key] = self.read_counts(self.files[key])
        return self.cache[key]

    def add_aliases(self, key, path, label):
        for alias in (label, str(path), str(path.resolve()), path.name, genome_stem(path.name)):
            if alias:
                self.aliases[alias].add(key)

    def resolve(self, name):
        if name in self.files:
            return name
        matches = self.aliases.get(name, set())
        if not matches and Path(name).name != name:
            matches = self.aliases.get(Path(name).name, set())
        if len(matches) > 1:
            raise ValueError(f'Ambiguous genome {name}: use an exact bin label from {sorted(matches)}')
        return next(iter(matches)) if matches else None


def load_quality(path, inventory):
    fields, rows = read_table(path)
    name_field = genome_column(fields, path)
    if not {'completeness', 'contamination'}.issubset(fields):
        raise ValueError(f'{path} needs completeness and contamination columns')
    quality = {}
    for row in rows:
        if not row[name_field]:
            raise ValueError(f'Missing genome identifier in {path}')
        key = inventory.resolve(row[name_field])
        if key is None:
            continue  # A full CheckM2 report may also contain unrelated genomes.
        if key in quality:
            raise ValueError(f'Duplicate quality entry for {key} in {path}')
        try:
            comp, cont = [float(row[field].replace('%', '')) for field in ('completeness', 'contamination')]
        except ValueError:
            raise ValueError(f'Invalid quality values for {key} in {path}') from None
        if not math.isfinite(comp) or not math.isfinite(cont) or not 0 <= comp <= 100 or cont < 0:
            raise ValueError(f'Invalid quality values for {key} in {path}')
        quality[key] = (comp, cont)
    return quality


def load_clusters(path, inventory):
    fields, rows = read_table(path)
    name_field = genome_column(fields, path)
    cluster_field = next((f for f in ('secondary_cluster', 'cluster', 'cluster_id') if f in fields), None)
    if cluster_field is None:
        raise ValueError(f'{path} needs secondary_cluster or cluster; primary_cluster alone is not supported')
    clusters = {}
    unavailable = []
    seen_names = {}
    for row in rows:
        name, cluster = row[name_field], row[cluster_field]
        if not name or cluster.lower() in ('', 'nan', 'none', 'null', 'na'):
            raise ValueError(f'Missing genome or cluster in {path}')
        if name in seen_names:
            if seen_names[name] != cluster:
                raise ValueError(f'Genome {name} assigned to multiple clusters in {path}')
            continue
        seen_names[name] = cluster
        key = inventory.resolve(name)
        if key is None:
            unavailable.append({'genome': name, 'cluster': cluster})
        elif key in clusters and clusters[key] != cluster:
            raise ValueError(f'Genome {key} assigned to multiple clusters in {path}')
        else:
            clusters[key] = cluster
    missing = inventory.primary - clusters.keys()
    if missing:
        raise ValueError(f'Cluster membership missing for selected MAGs: {sorted(missing)}')
    selected = {}
    for key in sorted(inventory.primary):
        cluster = clusters[key]
        if cluster in selected:
            raise ValueError(f'Multiple selected MAGs in cluster {cluster}: {selected[cluster]}, {key}; '
                             'replacement expects one input representative per cluster')
        selected[cluster] = key
    unavailable = [row for row in unavailable if row['cluster'] in selected]
    eligible = {key for key, cluster in clusters.items()
                if cluster in selected and key not in inventory.primary}
    return clusters, eligible, unavailable


def plan_resolution(inventory, quality, mode, rank_by='completeness', clusters=None, eligible=None,
                    taxonomy=None):
    if mode not in ('exclude-mags', 'remove-contigs') or rank_by not in ('completeness', 'quality-score'):
        raise ValueError('Unknown resolution mode or ranking')
    clusters, eligible = clusters or {}, eligible or set()
    missing = (inventory.primary | eligible) - quality.keys()
    if missing:
        raise ValueError(f'Quality missing for MAGs: {sorted(missing)}')

    def rank(key):
        comp, cont = quality[key]
        return (-(comp - 5 * cont), -comp, cont, key) if rank_by == 'quality-score' else (-comp, cont, key)

    ordered = sorted(inventory.primary, key=rank)
    plans, retained, owners = {}, {}, {}
    for key in ordered:
        keys = set(inventory.counts(key))
        shared = keys & owners.keys()
        blockers = sorted({owners[digest] for digest in shared})
        action = 'retained'
        remaining = keys
        if shared:
            remaining = set() if mode == 'exclude-mags' else keys - shared
            action = 'trimmed' if remaining else 'excluded'
        tax = taxonomy.assess(key, remaining, 'initial') if taxonomy and remaining else {}
        if tax and tax['taxonomy_status'] != 'coherent':
            remaining, action = set(), 'excluded'
        plans[key] = {'action': action, 'removed': shared if action == 'trimmed' else set(),
                      'blockers': blockers, 'replacement_for': '', 'replaced_by': '',
                      'exclusion_reason': tax.get('taxonomy_status', '') if action == 'excluded' and tax else (
                          'shared_contigs' if action == 'excluded' else '')}
        if remaining:
            retained[key] = remaining
            owners.update({digest: key for digest in remaining})
    affected = [key for key in ordered if plans[key]['action'] != 'retained']
    attempts, unresolved = [], []
    if clusters:
        alternatives_by_cluster = defaultdict(set)
        for key in eligible:
            alternatives_by_cluster[clusters[key]].add(key)
        ranked_candidates = {original: sorted(alternatives_by_cluster[clusters[original]] | {original}, key=rank)
                             for original in affected}
        pending = list(affected)
        round_no = 0
        remaining = []
        while pending:
            round_no += 1
            remaining = []
            for original in pending:
                cluster = clusters[original]
                # The original may become compatible after another cluster's replacement.
                for candidate in ranked_candidates[original]:
                    keys = set(inventory.counts(candidate))
                    blockers = sorted({owners[digest] for digest in keys
                                       if digest in owners and owners[digest] != original})
                    tax = taxonomy.assess(candidate, keys, 'replacement') if taxonomy and not blockers else {}
                    decision = 'shared_contigs' if blockers else tax.get('taxonomy_status', 'coherent')
                    attempts.append({'original': original, 'cluster': cluster, 'candidate': candidate,
                                     'round': round_no,
                                     'decision': 'accepted' if decision == 'coherent' else decision,
                                     'blocking_mags': json.dumps(blockers),
                                     'completeness': quality[candidate][0], 'contamination': quality[candidate][1],
                                     'taxonomy_status': tax.get('taxonomy_status', 'not_checked' if taxonomy else 'not_requested'),
                                     'taxonomy_resolved_fraction_bp': tax.get('taxonomy_resolved_fraction_bp', ''),
                                     'taxonomy_foreign_fraction': tax.get('taxonomy_foreign_fraction', '')})
                    if decision != 'coherent':
                        continue
                    for digest in retained.pop(original, set()):
                        del owners[digest]
                    retained[candidate] = keys
                    owners.update({digest: candidate for digest in keys})
                    if candidate == original:
                        plans[original].update(action='restored', removed=set(), exclusion_reason='')
                    else:
                        plans[original].update(action='replaced', removed=set(), replaced_by=candidate)
                        plans[candidate] = {'action': 'replacement', 'removed': set(), 'blockers': [],
                                            'replacement_for': original, 'replaced_by': '', 'exclusion_reason': ''}
                    break
                else:
                    remaining.append(original)
            if len(remaining) == len(pending):
                break  # No accepted change can unblock more clusters.
            pending = remaining
        for original in remaining:
            cluster = clusters[original]
            alternatives = alternatives_by_cluster[cluster]
            unresolved.append({'original': original, 'cluster': cluster,
                               'reason': 'no_compatible_candidate' if alternatives else 'no_available_candidate',
                               'final_action': plans[original]['action']})
    if mode == 'remove-contigs':
        # Whole replacements may free individual contigs even when the original
        # representative as a whole is still incompatible. Recover every such
        # contig in priority order, including nonempty remnants of excluded MAGs.
        for original in ordered:
            decision = plans[original]
            if decision['action'] in ('retained', 'restored', 'replaced'):
                continue
            keys = set(inventory.counts(original))
            available = {digest for digest in keys if digest not in owners or owners[digest] == original}
            if not available:
                continue
            tax = taxonomy.assess(original, available, 'contig_recovery') if taxonomy else {}
            if tax and tax['taxonomy_status'] != 'coherent':
                continue
            retained[original] = available
            owners.update({digest: original for digest in available})
            removed = keys - available
            decision.update(action='trimmed' if removed else 'restored', removed=removed,
                            blockers=sorted({owners[digest] for digest in removed}), exclusion_reason='')
        unresolved = [dict(row, final_action=plans[row['original']]['action']) for row in unresolved
                      if plans[row['original']]['action'] in ('trimmed', 'excluded')]
    # A final global assertion catches interactions between different replacements.
    final_owners = {}
    for key, keys in retained.items():
        if keys & final_owners.keys():
            raise RuntimeError('Internal resolution error: selected MAGs still share contigs')
        final_owners.update({digest: key for digest in keys})
    return {'plans': plans, 'retained': retained, 'attempts': attempts, 'unresolved': unresolved,
            'rank_by': rank_by, 'mode': mode, 'taxonomy': taxonomy}


def export_genome(path, target, removed, expected, removed_target=None):
    """Preserve complete headers and sequence orientation; save removed records."""
    observed = Counter()
    kept_count = removed_count = removed_bp = 0
    with target.open('w', encoding='utf-8') as kept:
        removed_stream = removed_target.open('w', encoding='utf-8') if removed_target else None
        try:
            for header, sequence in read_bin(path, full_header=True):
                digest = sequence_digest(sequence)
                observed[digest] += 1
                if digest in removed:
                    removed_count += 1
                    removed_bp += len(sequence)
                    stream = removed_stream
                else:
                    kept_count += 1
                    stream = kept
                if stream is not None:
                    stream.write('>' + header + '\n')
                    for start in range(0, len(sequence), 80):
                        stream.write(sequence[start:start + 80] + '\n')
        finally:
            if removed_stream:
                removed_stream.close()
    if observed != expected:
        raise ValueError(f'Input FASTA changed during resolution: {path}')
    return kept_count, removed_count, removed_bp


def write_resolution(before, inventory, quality, resolution, output, clusters=None, unavailable=None, metadata=None):
    """Publish only after validating the actual exported FASTAs in a temporary directory."""
    output = output.resolve()
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        raise ValueError(f'Output must be a new or empty directory: {output}')
    output.parent.mkdir(parents=True, exist_ok=True)
    plans = resolution['plans']
    taxonomy = resolution.get('taxonomy')
    clusters = clusters or {}
    names = {}
    counts = Counter(genome_stem(inventory.files[key].name) + '.fa' for key in plans)
    reserved, used = set(counts), set()
    for index, key in enumerate(sorted(plans), 1):
        name = genome_stem(inventory.files[key].name) + '.fa'
        if counts[name] > 1:
            name = f'{index:06d}__{name}'
            while name in reserved or name in used:
                name = 'resolved_' + name
        names[key] = name
        used.add(name)
    with tempfile.TemporaryDirectory(prefix='.echidna-resolution-', dir=str(output.parent)) as temp:
        root = Path(temp)
        summary = write_reports(before, root)
        for directory in ('retained', 'excluded_mags', 'removed_contigs'):
            (root / directory).mkdir()
        rows, retained_files = [], {}
        for key, decision in sorted(plans.items()):
            selected = key in resolution['retained']
            relative = Path('retained' if selected else 'excluded_mags') / names[key]
            removed_target = root / 'removed_contigs' / names[key] if decision['removed'] else None
            kept, removed_count, removed_bp = export_genome(
                inventory.files[key], root / relative, decision['removed'], inventory.counts(key), removed_target)
            comp, cont = quality[key]
            status = 'requires_reassessment' if decision['action'] == 'trimmed' else (
                'unchanged_input_estimate' if selected else 'not_retained')
            if selected:
                retained_files[key] = root / relative
            tax = taxonomy.evidence.assess_file(root / relative) if taxonomy else {}
            if taxonomy:
                taxonomy.checks.append(dict(genome=key, stage='exported_retained' if selected else 'exported_excluded', **tax))
                if selected and tax['taxonomy_status'] != 'coherent':
                    raise RuntimeError(f'Exported MAG {key} failed the taxonomy recheck; output not published')
            rows.append({
                'genome': key, 'source_file': str(inventory.files[key]), 'cluster': clusters.get(key, ''),
                'action': decision['action'], 'replacement_for': decision['replacement_for'],
                'exclusion_reason': decision.get('exclusion_reason', ''),
                'replaced_by': decision['replaced_by'], 'blocking_mags': json.dumps(decision['blockers']),
                'output_file': relative.as_posix(), 'input_completeness': comp, 'input_contamination': cont,
                'retained_contigs': kept if selected else 0, 'removed_contigs': removed_count,
                'removed_bp': removed_bp, 'quality_status': status,
                'requires_downstream_qc': decision['action'] in ('trimmed', 'replacement'),
                'final_completeness': comp if selected and status == 'unchanged_input_estimate' else '',
                'final_contamination': cont if selected and status == 'unchanged_input_estimate' else '',
                'taxonomy_status': tax.get('taxonomy_status', 'not_requested'),
                'taxonomy_rank': tax.get('taxonomy_rank', ''),
                'taxonomy_resolved_fraction_bp': tax.get('taxonomy_resolved_fraction_bp', ''),
                'taxonomy_foreign_fraction': tax.get('taxonomy_foreign_fraction', ''),
                'taxonomy_majority': tax.get('taxonomy_majority', ''),
            })
        after = analyze_bins(retained_files, allow_empty=True) if not retained_files else analyze_bins(retained_files)
        if after['shared']:
            raise RuntimeError('Exported MAGs still share contigs; resolution output not published')
        for key, row in after['bins'].items():
            row['file_path'] = str(output / 'retained' / names[key])
        after_summary = write_reports(after, root / 'after')
        _write_tsv(root / 'resolution.tsv', list(rows[0]), rows)
        _write_tsv(root / 'replacement_attempts.tsv',
                   ['original', 'cluster', 'candidate', 'round', 'decision', 'blocking_mags', 'completeness', 'contamination',
                    'taxonomy_status', 'taxonomy_resolved_fraction_bp', 'taxonomy_foreign_fraction'],
                   resolution['attempts'])
        if taxonomy:
            _write_tsv(root / 'taxonomy_checks.tsv', list(taxonomy.checks[0]), taxonomy.checks)
        _write_tsv(root / 'unresolved_clusters.tsv', ['original', 'cluster', 'reason', 'final_action'], resolution['unresolved'])
        _write_tsv(root / 'unavailable_cluster_members.tsv', ['genome', 'cluster'], unavailable or [])
        summary['resolution'] = {
            'mode': resolution['mode'], 'rank_by': resolution['rank_by'],
            'n_retained_mags': len(retained_files),
            'n_excluded_original_mags': sum(d['action'] in ('excluded', 'replaced') for d in plans.values()),
            'n_trimmed_mags': sum(d['action'] == 'trimmed' for d in plans.values()),
            'n_replacements': sum(d['action'] == 'replacement' for d in plans.values()),
            'n_restored_originals': sum(d['action'] == 'restored' for d in plans.values()),
            'n_unresolved_clusters': len(resolution['unresolved']),
            'n_unavailable_cluster_members': len(unavailable or []),
            'final_n_shared_sequences': after_summary['n_shared_sequences'],
            'final_overlap_check': 'passed',
            'final_taxonomy_check': ('passed' if retained_files else 'no_retained_mags') if taxonomy else 'not_requested',
            'taxonomy_retained_mags_checked': len(retained_files) if taxonomy else 0,
            'taxonomy_rank': taxonomy.evidence.rank if taxonomy else None,
            'min_taxonomy_resolved_fraction_bp': taxonomy.evidence.min_resolved_frac if taxonomy else None,
            'min_foreign_fraction': taxonomy.evidence.min_foreign_frac if taxonomy else None,
            'quality_reassessment_required': any(d['action'] == 'trimmed' for d in plans.values()),
            'downstream_qc_required': any(d['action'] in ('trimmed', 'replacement') for d in plans.values()),
            'metadata_sources': metadata or {},
            'selection': 'Deterministic greedy priority; not a globally optimal subset or rerun of dereplication',
        }
        (root / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
        # Do not expose partial results or replace any user files on a validation failure.
        if output.exists():
            if not output.is_dir() or any(output.iterdir()):
                raise ValueError(f'Output must still be empty: {output}')
            output.rmdir()
        os.rename(root, output)
    return summary
