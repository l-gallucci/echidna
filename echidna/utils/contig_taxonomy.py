"""
Per-contig taxonomy utilities (independent of GUNC).

A bin whose contigs resolve to two or more taxa at a chosen rank is a
taxonomic chimera: contigs of different taxonomy do not belong in one MAG at
any rank. This mirrors the GUNC per-contig signal but uses a contig-level
classifier (Kraken2, Kaiju or MMseqs2), so it is a fully independent
corroborating signal.

Two ways to resolve a contig to the requested rank:
- lineage: an MMseqs2 TSV written with `--tax-lineage 1` (or any GTDB lineage
  string) carries `d_/p_/c_...`, so the rank is read straight from it - no
  nodes.dmp, and ranks stay in the DB's system (GTDB). Preferred.
- taxid + nodes.dmp: a Kraken2/Kaiju (or lineage-less MMseqs2) taxid is walked UP
  the NCBI tree to the rank. A contig classified only above the rank stays
  unresolved (conservative: it cannot prove a split).

Then, per bin, the distinct rank-level taxa are counted (>= 2 flags a chimera)
and the minority contigs are named as the likely foreign/contaminant ones.

No third-party Python dependencies; nodes.dmp/names.dmp ship with every Kraken2
or Kaiju database, and the lineage path needs neither.
"""

import re
from pathlib import Path
from typing import Callable, Dict, List, Optional, Set, Tuple

from .fasta import read_fasta_streaming
from .logging import log_info, log_warning


# NCBI calls the top rank 'superkingdom'; GTDB (and our other commands) call it
# 'domain'. Accept the GTDB word and translate so --rank domain works too.
_RANK_ALIASES = {'domain': 'superkingdom'}

VALID_RANKS = ('domain', 'superkingdom', 'phylum', 'class', 'order',
               'family', 'genus', 'species')

# Single rank letter used in semicolon lineage strings (MMseqs `--tax-lineage 1`
# and GTDB both use d/p/c/o/f/g/s; MMseqs writes one underscore, GTDB two).
_RANK_LETTER = {'domain': 'd', 'superkingdom': 'd', 'phylum': 'p', 'class': 'c',
                'order': 'o', 'family': 'f', 'genus': 'g', 'species': 's'}


def parse_mmseqs_taxonomy(tsv_file: Path) -> Dict[str, int]:
    """
    Parse an MMseqs2 taxonomy TSV (``mmseqs createtsv`` on a ``mmseqs taxonomy``
    result). The stable leading columns are:

        query_id    taxid    rank    name    [lineage / fragment columns...]

    Only query_id and taxid are used, so the walk to the requested rank is done
    with nodes.dmp here (works for an NCBI or a GTDB-derived database alike).
    Returns {contig: taxid}; unassigned contigs (taxid 0) are kept as 0.
    """
    contig_taxid: Dict[str, int] = {}
    with open(tsv_file) as f:
        for lineno, line in enumerate(f, start=1):
            if not line.strip():
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 2:
                continue
            query = fields[0].strip()
            try:
                taxid = int(fields[1].strip())
            except ValueError:
                log_warning(f"Non-integer taxid on MMseqs2 line {lineno}, skipping")
                continue
            contig_taxid[query] = taxid
    log_info(f"Loaded MMseqs2 taxonomy for {len(contig_taxid):,} contigs "
             f"from {tsv_file.name}")
    return contig_taxid


def parse_mmseqs_lineage(tsv_file: Path) -> Dict[str, str]:
    """
    Parse the lineage column of an MMseqs2 taxonomy TSV produced with
    ``mmseqs taxonomy --tax-lineage 1``. The leading columns are
    ``query_id  taxid  rank  name`` and the ';'-separated lineage
    (``d_Bacteria;p_...;c_...``) is appended LAST - but MMseqs puts a varying
    number of count/support columns between ``name`` and it, so the lineage is
    located by content (the last field containing ';'), not a fixed index.

    Using the lineage (instead of a taxid + nodes.dmp walk) is what lets this
    work with a GTDB MMseqs database, which carries no NCBI nodes.dmp, and keeps
    the ranks in GTDB's own system. Returns {contig: lineage_string}; contigs
    without a lineage are omitted.
    """
    contig_lineage: Dict[str, str] = {}
    with open(tsv_file) as f:
        for line in f:
            if not line.strip():
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 2:
                continue
            contig = fields[0].strip()
            # lineage = last field that looks like a ';'-separated path
            lineage = ''
            for fld in reversed(fields[1:]):
                if ';' in fld:
                    lineage = fld.strip()
                    break
            if contig and lineage:
                contig_lineage[contig] = lineage
    if contig_lineage:
        log_info(f"Loaded MMseqs2 lineage for {len(contig_lineage):,} contigs "
                 f"from {tsv_file.name}")
    return contig_lineage


def extract_rank_from_lineage(lineage: str, rank: str) -> Optional[str]:
    """
    Pull one rank's name out of a semicolon lineage string. Tolerant to one or
    two underscores after the rank letter (MMseqs `c_Gammaproteobacteria`, GTDB
    `c__Gammaproteobacteria`). Returns None when the rank is absent (contig
    classified only above it) so it is left unresolved, not counted as a split.
    """
    letter = _RANK_LETTER.get(rank)
    if not letter:
        return None
    pat = re.compile(r'^' + letter + r'_+(.+)$')
    for part in lineage.split(';'):
        m = pat.match(part.strip())
        if m:
            name = m.group(1).strip()
            return name or None
    return None


def parse_nodes_dmp(nodes_file: Path) -> Dict[int, Tuple[int, str]]:
    """
    Parse an NCBI taxonomy nodes.dmp into {taxid: (parent_taxid, rank)}.

    Format (fields separated by '\\t|\\t', rows end with '\\t|'):
        taxid | parent_taxid | rank | ...
    """
    nodes: Dict[int, Tuple[int, str]] = {}
    with open(nodes_file) as f:
        for line in f:
            if not line.strip():
                continue
            fields = [c.strip() for c in line.rstrip('\n').rstrip('|').split('|')]
            if len(fields) < 3:
                continue
            try:
                taxid = int(fields[0])
                parent = int(fields[1])
            except ValueError:
                continue
            nodes[taxid] = (parent, fields[2])
    log_info(f"Loaded taxonomy tree for {len(nodes):,} nodes from {nodes_file.name}")
    return nodes


def resolve_at_rank(taxid: int, target_rank: str,
                    nodes: Dict[int, Tuple[int, str]]) -> Optional[int]:
    """
    Walk a taxid up the tree until a node of ``target_rank`` is reached; return
    that ancestor's taxid, or None when the rank is not on the lineage (the
    taxid is already above the rank, or is unknown/unclassified).

    Guards against missing parents and cycles (root is its own parent).
    """
    if taxid <= 0:
        return None
    rank = _RANK_ALIASES.get(target_rank, target_rank)
    seen: Set[int] = set()
    cur = taxid
    while cur and cur not in seen:
        seen.add(cur)
        node = nodes.get(cur)
        if node is None:
            return None
        parent, node_rank = node
        if node_rank == rank:
            return cur
        if parent == cur:          # reached root without hitting the rank
            return None
        cur = parent
    return None


def assign_contigs_to_bins(bins_dir: Path,
                           extensions: List[str]) -> Dict[str, List[str]]:
    """
    Map each bin to its contig names (header first word) by reading the bin
    FASTA files. Returns {bin_name: [contig, ...]}.
    """
    bin_files: List[Path] = []
    for ext in extensions:
        bin_files.extend(bins_dir.rglob(f"*.{ext}"))
    bin_files = sorted(set(bin_files))

    bins: Dict[str, List[str]] = {}
    for bin_file in bin_files:
        contigs = [name for name, _ in read_fasta_streaming(bin_file)]
        bins[bin_file.stem] = contigs
    log_info(f"Mapped {sum(len(c) for c in bins.values()):,} contigs "
             f"across {len(bins):,} bins")
    return bins


def _taxon_name(taxid: int, names: Optional[Dict[int, Tuple[str, str]]]) -> str:
    if names and taxid in names:
        return names[taxid][1]
    return f"taxid:{taxid}"


def _summarize_bins(
    bin_contigs: Dict[str, List[str]],
    contig_key: Dict[str, object],
    label_fn: Callable[[object], str],
    target_rank: str,
) -> Dict[str, Dict]:
    """
    Shared core. ``contig_key`` maps each contig to a hashable rank-level key
    (or None when the contig cannot be resolved at the rank). The key may be a
    taxid (nodes.dmp path) or a taxon-name string (lineage path); ``label_fn``
    turns a key into its display name.

    Like MAGpurify2's taxonomy module, the bin's MAJORITY rank-taxon (most
    contigs) is its identity and the remaining resolved contigs are reported as
    foreign/minority contigs - the likely contaminants - so a single misplaced
    contig is named, not just counted.
    """
    out: Dict[str, Dict] = {}
    for bin_name, contigs in bin_contigs.items():
        per_contig: Dict[str, Tuple[Optional[object], str]] = {}
        keys: Set[object] = set()
        key_count: Dict[object, int] = {}
        for contig in contigs:
            key = contig_key.get(contig)
            if key is None:
                per_contig[contig] = (None, 'unresolved')
            else:
                keys.add(key)
                key_count[key] = key_count.get(key, 0) + 1
                per_contig[contig] = (key, label_fn(key))

        # Majority key = the bin's identity (ties broken deterministically on the
        # display label). Everything else resolved is foreign.
        majority_key: Optional[object] = None
        if key_count:
            majority_key = max(sorted(key_count, key=label_fn),
                               key=key_count.get)
        foreign_contigs = sorted(
            (c, v[1]) for c, v in per_contig.items()
            if v[0] is not None and v[0] != majority_key
        )

        out[bin_name] = {
            'n_contigs': len(contigs),
            'n_contigs_resolved': sum(1 for v in per_contig.values() if v[0] is not None),
            'n_taxa': len(keys),
            'taxa': sorted(label_fn(k) for k in keys),
            'chimera': len(keys) >= 2,
            'majority_taxon': majority_key,
            'majority_name': (label_fn(majority_key)
                              if majority_key is not None else 'N/A'),
            'foreign_contigs': foreign_contigs,
            'per_contig': per_contig,
        }
    n_chimeric = sum(1 for r in out.values() if r['chimera'])
    log_info(f"Per-contig taxonomy at rank '{target_rank}': "
             f"{n_chimeric:,}/{len(out):,} bins span >= 2 taxa")
    return out


def summarize_bin_taxa(
    bin_contigs: Dict[str, List[str]],
    contig_taxid: Dict[str, int],
    target_rank: str,
    nodes: Dict[int, Tuple[int, str]],
    names: Optional[Dict[int, Tuple[str, str]]] = None,
) -> Dict[str, Dict]:
    """
    Resolve each contig's taxid up to ``target_rank`` via nodes.dmp (Kraken2 /
    Kaiju / MMseqs-without-lineage path), then summarize per bin. ``majority_taxon``
    is a taxid. See _summarize_bins for the shape.
    """
    contig_key = {
        contig: resolve_at_rank(contig_taxid.get(contig, 0), target_rank, nodes)
        for contigs in bin_contigs.values() for contig in contigs
    }
    return _summarize_bins(
        bin_contigs, contig_key, lambda t: _taxon_name(t, names), target_rank)


def summarize_bin_lineage(
    bin_contigs: Dict[str, List[str]],
    contig_lineage: Dict[str, str],
    target_rank: str,
) -> Dict[str, Dict]:
    """
    Resolve each contig to its ``target_rank`` taxon name from a lineage string
    (MMseqs `--tax-lineage 1` / GTDB), then summarize per bin. No nodes.dmp
    needed; ``majority_taxon`` is the rank name string. See _summarize_bins.
    """
    contig_key = {
        contig: extract_rank_from_lineage(contig_lineage.get(contig, ''), target_rank)
        for contigs in bin_contigs.values() for contig in contigs
    }
    return _summarize_bins(bin_contigs, contig_key, lambda s: s, target_rank)
