"""Exact full-contig sharing between genome FASTAs; no alignments or refinement."""
import csv
from dataclasses import dataclass
import hashlib
from itertools import combinations
import json
from pathlib import Path
from typing import Dict, Iterator, List, Tuple

from .fasta import open_file

DNA = frozenset('ACGTRYSWKMBDHVN')
COMPLEMENT = str.maketrans('ACGTRYSWKMBDHVN', 'TGCAYRSWMKVHDBN')


@dataclass(frozen=True)
class Contig:
    bin_id: str
    contig_id: str
    length: int
    digest: str


def sequence_digest(sequence: str) -> str:
    """SHA256 of validated upper-case DNA in a canonical strand orientation."""
    reverse = sequence.translate(COMPLEMENT)[::-1]
    return hashlib.sha256(min(sequence, reverse).encode('ascii')).hexdigest()


def read_bin(path: Path) -> Iterator[Tuple[str, str]]:
    """Strict FASTA reader: reject empty records, bad DNA and duplicate IDs."""
    current = None
    parts = []
    seen = set()
    with open_file(path, 'rt') as stream:
        for line_no, line in enumerate(stream, 1):
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current is not None:
                    if not parts:
                        raise ValueError(f'Empty sequence for {current} in {path}')
                    yield current, ''.join(parts)
                header = line[1:].split()
                if not header:
                    raise ValueError(f'Empty FASTA header at {path}:{line_no}')
                current = header[0]
                if current in seen:
                    raise ValueError(f'Duplicate contig ID {current} within {path}')
                seen.add(current)
                parts = []
            else:
                if current is None:
                    raise ValueError(f'Sequence before first FASTA header at {path}:{line_no}')
                sequence = ''.join(line.split()).upper()
                invalid = set(sequence) - DNA
                if invalid:
                    raise ValueError(f'Invalid DNA characters {sorted(invalid)} at {path}:{line_no}')
                parts.append(sequence)
    if current is None:
        raise ValueError(f'No FASTA records in {path}')
    if not parts:
        raise ValueError(f'Empty sequence for {current} in {path}')
    yield current, ''.join(parts)


def discover_bins(roots: List[Path], recursive: bool = False,
                  extensions: str = 'fa,fna,fasta') -> Dict[str, Path]:
    """Stable file labels; reject a path included twice through overlapping roots."""
    suffixes = {ext.strip().lower().lstrip('.') for ext in extensions.split(',') if ext.strip().lstrip('.')}
    if not suffixes or any('/' in ext or '\\' in ext for ext in suffixes):
        raise ValueError('Provide FASTA extensions, for example fa,fna,fasta')
    files = {}
    seen = set()
    for root_no, root in enumerate(roots, 1):
        root = Path(root).resolve()
        if not root.is_dir():
            raise ValueError(f'Bins directory not found: {root}')
        found = 0
        for path in sorted(root.rglob('*') if recursive else root.iterdir()):
            if not path.is_file():
                continue
            name = path.name.lower()
            for compression in ('.gzip', '.gz'):
                if name.endswith(compression):
                    name = name[:-len(compression)]
                    break
            if not any(name.endswith('.' + ext) for ext in suffixes):
                continue
            if path in seen:
                raise ValueError(f'Input FASTA included more than once: {path}')
            seen.add(path)
            relative = path.relative_to(root).as_posix()
            label = f'set{root_no}/{relative}' if len(roots) > 1 else relative
            files[label] = path
            found += 1
        if not found:
            raise ValueError(f'No matching FASTA files in {root}; check --extensions or --recursive')
    return dict(sorted(files.items()))


def analyze_bins(files: Dict[str, Path]) -> dict:
    """Index sequence hashes and IDs, then enumerate only candidate bin pairs.

    Sequences are discarded after hashing. Metadata scales with contig count;
    pair generation enumerates only bins sharing an exact sequence.
    """
    if not files:
        raise ValueError('No bins supplied')
    records = []
    by_sequence = {}
    by_id = {}
    summaries = {}
    unique_bp = {}
    unique_counts = {}
    for bin_id, path in sorted(files.items()):
        start = len(records)
        unique = {}
        for contig_id, sequence in read_bin(path):
            digest = sequence_digest(sequence)
            index = len(records)
            records.append(Contig(bin_id, contig_id, len(sequence), digest))
            by_sequence.setdefault(digest, []).append(index)
            by_id.setdefault(contig_id, []).append(index)
            unique[digest] = len(sequence)
        own = records[start:]
        summaries[bin_id] = {
            'bin': bin_id, 'file_path': str(path), 'n_contigs': len(own),
            'total_bp': sum(rec.length for rec in own), 'n_unique_sequences': len(unique),
            'n_shared_contigs': 0, 'shared_bp': 0, 'shared_fraction_bp': 0.0,
            'n_conflicting_ids': 0,
        }
        unique_bp[bin_id] = sum(unique.values())
        unique_counts[bin_id] = len(unique)
    shared = {}
    conflicts = {}
    pairs = {}

    def pair(a, b):
        key = tuple(sorted((a, b)))
        if key not in pairs:
            pairs[key] = {'bin_a': key[0], 'bin_b': key[1],
                          'n_shared_sequences': 0, 'shared_unique_bp': 0}
        return pairs[key]

    for digest, indices in sorted(by_sequence.items()):
        bins = sorted({records[i].bin_id for i in indices})
        if len(bins) < 2:
            continue
        shared[digest] = indices
        for i in indices:
            rec = records[i]
            summaries[rec.bin_id]['n_shared_contigs'] += 1
            summaries[rec.bin_id]['shared_bp'] += rec.length
        length = records[indices[0]].length
        for a, b in combinations(bins, 2):
            row = pair(a, b)
            row['n_shared_sequences'] += 1
            row['shared_unique_bp'] += length
    for contig_id, indices in sorted(by_id.items()):
        if len({records[i].bin_id for i in indices}) < 2:
            continue
        if len({records[i].digest for i in indices}) < 2:
            continue
        conflicts[contig_id] = indices
        for i in indices:
            summaries[records[i].bin_id]['n_conflicting_ids'] += 1
    for row in summaries.values():
        row['shared_fraction_bp'] = row['shared_bp'] / row['total_bp']
    for (a, b), row in pairs.items():
        row['unique_bp_a'] = unique_bp[a]
        row['unique_bp_b'] = unique_bp[b]
        row['shared_fraction_a'] = row['shared_unique_bp'] / unique_bp[a]
        row['shared_fraction_b'] = row['shared_unique_bp'] / unique_bp[b]
        row['sequence_jaccard'] = row['n_shared_sequences'] / (
            unique_counts[a] + unique_counts[b] - row['n_shared_sequences'])
    return {'records': records, 'shared': shared, 'conflicts': conflicts,
            'bins': summaries, 'pairs': dict(sorted(pairs.items()))}


def _write_tsv(path, fields, rows):
    with path.open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        for row in rows:
            writer.writerow({key: (f'{value:.6f}' if isinstance(value, float) else value)
                             for key, value in row.items()})


def write_reports(result: dict, output: Path) -> dict:
    """Write report tables after all inputs have been validated; never alter bins."""
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        raise ValueError(f'Output must be a new or empty directory: {output}')
    output.mkdir(parents=True, exist_ok=True)
    records = result['records']

    def memberships(groups, key_name):
        for key, indices in groups.items():
            n_bins = len({records[i].bin_id for i in indices})
            for i in indices:
                rec = records[i]
                row = {key_name: key, 'n_bins': n_bins, 'bin': rec.bin_id,
                       'contig_id': rec.contig_id, 'length_bp': rec.length}
                if key_name == 'conflicting_id':
                    row['sequence_sha256'] = rec.digest
                yield row

    _write_tsv(output / 'shared_contigs.tsv',
               ['sequence_sha256', 'n_bins', 'bin', 'contig_id', 'length_bp'],
               memberships(result['shared'], 'sequence_sha256'))
    _write_tsv(output / 'id_conflicts.tsv',
               ['conflicting_id', 'n_bins', 'bin', 'contig_id', 'length_bp', 'sequence_sha256'],
               memberships(result['conflicts'], 'conflicting_id'))
    _write_tsv(output / 'bin_pairs.tsv',
               ['bin_a', 'bin_b', 'n_shared_sequences', 'shared_unique_bp',
                'unique_bp_a', 'unique_bp_b', 'shared_fraction_a', 'shared_fraction_b', 'sequence_jaccard'],
               result['pairs'].values())
    _write_tsv(output / 'bin_summary.tsv',
               ['bin', 'file_path', 'n_contigs', 'total_bp', 'n_unique_sequences',
                'n_shared_contigs', 'shared_bp', 'shared_fraction_bp', 'n_conflicting_ids'],
               result['bins'].values())
    summary = {
        'n_bins': len(result['bins']), 'n_contigs': len(records),
        'total_bp': sum(row['total_bp'] for row in result['bins'].values()),
        'n_shared_sequences': len(result['shared']),
        'n_bins_with_shared_contigs': sum(row['n_shared_contigs'] > 0 for row in result['bins'].values()),
        'n_pairs_with_shared_contigs': sum(row['n_shared_sequences'] > 0 for row in result['pairs'].values()),
        'n_conflicting_ids': len(result['conflicts']),
        'matching': 'SHA256 of full-length, upper-case, whitespace-free IUPAC DNA; reverse-complement invariant',
        'contig_ids': 'Exact first token of each FASTA header; conflicting sequences reported separately',
        'pair_bp_denominator': 'Sum of lengths of distinct canonical sequences in each bin',
        'bin_bp_denominator': 'Sum of lengths of all contig records in that bin',
    }
    (output / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
    return summary
