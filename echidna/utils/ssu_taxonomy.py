"""
16S/SSU rRNA taxonomy utilities.

An orthogonal, independent chimera/misbinning signal: classify a MAG's 16S
(SSU) rRNA gene and check it agrees with the genome-level GTDB-Tk taxonomy. A
16S that resolves to a different taxon than the genome means misbinning or a
contaminating rRNA from another organism.

Flow
----
1. barrnap (run-mimag-rna --outseq) writes the rRNA genes per genome; we keep
   the 16S (or 18S) record.
2. blastn the 16S against an SSU reference DB (GTDB ssu_all, or SILVA).
3. Read the best hit's taxon at a chosen rank and compare it to the genome's
   GTDB-Tk taxon at the same rank.

Limitations (honest): only MAGs with a recovered ~full-length 16S get checked
(rRNA is often fragmented/unbinned), and it is confirmatory, not exhaustive - a
chimeric bin may carry only the dominant member's 16S. With a SILVA DB the
taxon names differ from GTDB (e.g. Proteobacteria vs Pseudomonadota), so the
agreement is best-effort; a GTDB SSU DB compares 1:1.

Pure parsers here (blastn is shelled out in the command); no third-party deps.
"""

import re
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from .fasta import read_fasta_streaming
from .contig_taxonomy import extract_rank_from_lineage
from .logging import log_info, log_warning, log_error
from .validation import check_dependencies


# barrnap --outseq headers look like "16S_rRNA::contig:123-1600(+)".
_SSU_PREFIXES = ('16S_rRNA', '18S_rRNA')

# blastn tabular columns we request (stitle LAST so it can contain spaces/;).
BLAST_OUTFMT = '6 qseqid sseqid pident length slen bitscore stitle'

# SILVA lineage is positional (no rank prefixes): Domain;Phylum;Class;Order;...
_SILVA_RANK_INDEX = {'domain': 0, 'superkingdom': 0, 'phylum': 1, 'class': 2,
                     'order': 3, 'family': 4, 'genus': 5, 'species': 6}


def extract_ssu_records(rrna_fasta: Path) -> List[Tuple[str, str]]:
    """Return the (header, seq) of SSU (16S/18S) rRNA genes in a barrnap
    --outseq FASTA. Empty when none found (or file missing)."""
    if not Path(rrna_fasta).exists():
        return []
    out = []
    for name, seq in read_fasta_streaming(rrna_fasta):
        if name.startswith(_SSU_PREFIXES):
            out.append((name, seq))
    return out


def longest_ssu(records: List[Tuple[str, str]]) -> Optional[Tuple[str, str]]:
    """Pick the longest SSU record (the most complete copy) or None."""
    if not records:
        return None
    return max(records, key=lambda r: len(r[1]))


def check_blastn(exe: str = 'blastn') -> bool:
    deps = check_dependencies([exe])
    if not deps.get(exe):
        log_error(f"{exe} not found in PATH. check-ssu-taxonomy needs BLAST+ "
                  "(blastn) - run it in the echidna env.")
        return False
    return True


def run_blastn(query_fasta: Path, db: str, out_tsv: Path, threads: int = 1,
               task: str = 'megablast', max_target_seqs: int = 5,
               exe: str = 'blastn') -> int:
    """Run blastn (outfmt BLAST_OUTFMT) and write out_tsv. Returns exit code."""
    cmd = [exe, '-query', str(query_fasta), '-db', str(db),
           '-out', str(out_tsv), '-outfmt', BLAST_OUTFMT,
           '-max_target_seqs', str(max_target_seqs),
           '-num_threads', str(threads), '-task', task]
    try:
        subprocess.run(cmd, check=True, stdin=subprocess.DEVNULL)
        return 0
    except subprocess.CalledProcessError as e:
        log_error(f"blastn failed (exit {e.returncode})")
        return e.returncode


def parse_blast_best(out_tsv: Path) -> Dict[str, Dict]:
    """Best hit per query (highest bitscore) from a BLAST_OUTFMT TSV.

    Returns {qseqid: {pident, length, slen, bitscore, stitle}}.
    """
    best: Dict[str, Dict] = {}
    if not Path(out_tsv).exists():
        return best
    with open(out_tsv) as f:
        for line in f:
            if not line.strip():
                continue
            # stitle is the last requested field; keep it intact.
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 7:
                continue
            q, s, pident, length, slen, bitscore = parts[:6]
            stitle = '\t'.join(parts[6:]).strip()
            try:
                bs = float(bitscore)
            except ValueError:
                continue
            prev = best.get(q)
            if prev is None or bs > prev['bitscore']:
                best[q] = {
                    'sseqid': s,
                    'pident': _safe_float(pident),
                    'length': _safe_int(length),
                    'slen': _safe_int(slen),
                    'bitscore': bs,
                    'stitle': stitle,
                }
    return best


def ssu_taxon_at_rank(stitle: str, rank: str, scheme: str) -> Optional[str]:
    """Extract the hit's taxon name at `rank` from its subject title.

    scheme='gtdb' : lineage carries d__/p__/c__ prefixes (reuse the lineage
                    parser). scheme='silva' : positional, no prefixes, with a
                    leading accession token to drop.
    """
    if not stitle:
        return None
    if scheme == 'gtdb':
        return extract_rank_from_lineage(stitle, rank)
    # SILVA: "<accession> Domain;Phylum;Class;..." -> drop the accession token
    idx = _SILVA_RANK_INDEX.get(rank)
    if idx is None:
        return None
    body = stitle.split(None, 1)
    lineage = body[1] if len(body) == 2 else body[0]
    parts = [p.strip() for p in lineage.split(';') if p.strip()]
    if idx < len(parts):
        return parts[idx] or None
    return None


def normalize_taxon(name: Optional[str]) -> str:
    """Normalise a taxon name for comparison: lowercase, drop a GTDB
    polyphyly suffix (Firmicutes_A -> firmicutes). Empty for None/Unclassified."""
    if not name:
        return ''
    n = name.strip()
    if n.lower() in ('unclassified', 'na', 'n/a', ''):
        return ''
    n = re.sub(r'_[A-Z]+$', '', n)
    return n.lower()


def taxa_agree(ssu_taxon: Optional[str], genome_taxon: Optional[str]) -> Optional[bool]:
    """True/False agreement, or None when either side is unknown (cannot judge)."""
    a, b = normalize_taxon(ssu_taxon), normalize_taxon(genome_taxon)
    if not a or not b:
        return None
    return a == b


def _safe_float(v) -> Optional[float]:
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def _safe_int(v) -> Optional[int]:
    try:
        return int(v)
    except (TypeError, ValueError):
        return None
