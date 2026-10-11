"""
Per-contig taxonomy utilities (independent of GUNC).

A bin whose contigs resolve to two or more taxa at a chosen rank is a
taxonomic chimera: contigs of different taxonomy do not belong in one MAG at
any rank. This mirrors the GUNC per-contig signal but uses a contig-level
classifier (Kraken2 or Kaiju) and the NCBI taxonomy tree, so it is a fully
independent corroborating signal.

Flow
----
1. A contig-level classifier assigns one taxid per contig (Kraken2/Kaiju output
   where the "sequence id" is the contig name; parsing is shared with the
   read-level code since the file shape is identical).
2. nodes.dmp gives each taxid its parent and rank, so a contig's LCA taxid is
   walked UP to the requested rank (e.g. order). A contig classified only above
   the rank cannot be resolved there and is left out (conservative: it cannot
   prove a split).
3. Per bin, count the distinct rank-level taxa across its contigs. >= 2 is a
   chimera flag.

No third-party Python dependencies; nodes.dmp/names.dmp ship with every Kraken2
or Kaiju database.
"""

from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from .fasta import read_fasta_streaming
from .logging import log_info, log_warning


# NCBI calls the top rank 'superkingdom'; GTDB (and our other commands) call it
# 'domain'. Accept the GTDB word and translate so --rank domain works too.
_RANK_ALIASES = {'domain': 'superkingdom'}

VALID_RANKS = ('domain', 'superkingdom', 'phylum', 'class', 'order',
               'family', 'genus', 'species')


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


def summarize_bin_taxa(
    bin_contigs: Dict[str, List[str]],
    contig_taxid: Dict[str, int],
    target_rank: str,
    nodes: Dict[int, Tuple[int, str]],
    names: Optional[Dict[int, Tuple[str, str]]] = None,
) -> Dict[str, Dict]:
    """
    For each bin, resolve every contig's taxid up to ``target_rank`` and count
    the distinct rank-level taxa.

    Like MAGpurify2's taxonomy module, the bin's MAJORITY rank-taxon (most
    contigs) is treated as the bin's identity and the remaining resolved contigs
    are reported as foreign/minority contigs - the likely contaminants - so a
    single misplaced contig in an otherwise clean MAG is named, not just counted.

    Returns {bin: {n_contigs, n_contigs_resolved, n_taxa, taxa (sorted names),
    chimera (bool, n_taxa >= 2), majority_taxon (taxid|None), majority_name,
    foreign_contigs [(contig, rank_name), ...], per_contig {contig:
    (rank_taxid|None, rank_name|'unresolved')}}}.
    """
    out: Dict[str, Dict] = {}
    for bin_name, contigs in bin_contigs.items():
        per_contig: Dict[str, Tuple[Optional[int], str]] = {}
        rank_taxids: Set[int] = set()
        taxid_contig_count: Dict[int, int] = {}
        for contig in contigs:
            taxid = contig_taxid.get(contig, 0)
            rank_taxid = resolve_at_rank(taxid, target_rank, nodes)
            if rank_taxid is None:
                per_contig[contig] = (None, 'unresolved')
            else:
                rank_taxids.add(rank_taxid)
                taxid_contig_count[rank_taxid] = taxid_contig_count.get(rank_taxid, 0) + 1
                per_contig[contig] = (rank_taxid, _taxon_name(rank_taxid, names))

        # Majority rank-taxon = the bin's identity (ties broken by taxid for a
        # stable, reproducible pick). Everything else resolved is foreign.
        majority_taxid: Optional[int] = None
        if taxid_contig_count:
            majority_taxid = max(sorted(taxid_contig_count),
                                 key=taxid_contig_count.get)
        foreign_contigs = sorted(
            (c, v[1]) for c, v in per_contig.items()
            if v[0] is not None and v[0] != majority_taxid
        )

        taxa = sorted(_taxon_name(t, names) for t in rank_taxids)
        out[bin_name] = {
            'n_contigs': len(contigs),
            'n_contigs_resolved': sum(1 for v in per_contig.values() if v[0] is not None),
            'n_taxa': len(rank_taxids),
            'taxa': taxa,
            'chimera': len(rank_taxids) >= 2,
            'majority_taxon': majority_taxid,
            'majority_name': (_taxon_name(majority_taxid, names)
                              if majority_taxid is not None else 'N/A'),
            'foreign_contigs': foreign_contigs,
            'per_contig': per_contig,
        }
    n_chimeric = sum(1 for r in out.values() if r['chimera'])
    log_info(f"Per-contig taxonomy at rank '{target_rank}': "
             f"{n_chimeric:,}/{len(out):,} bins span >= 2 taxa")
    return out
