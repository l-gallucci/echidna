"""
Canonical genome-name matching across tools.

Every per-MAG signal is joined by genome name, but CheckM2 (`Name`), GTDB-Tk
(`user_genome`), GUNC (`genome`, file-suffix stripped), MMseqs/classify-contigs
and the bin FASTA stems can spell the same genome differently (extensions,
`.contigs`, case). A mismatch makes a signal silently become N/A. These helpers
normalise both sides to one canonical key and let callers report what did NOT
match, so a naming problem is visible instead of silent.
"""

from typing import Dict, Iterable, List, Optional, Tuple


_FASTA_EXTS = ('.fa.gz', '.fasta.gz', '.fna.gz', '.fa', '.fasta', '.fna',
               '.fsa', '.gz')


def canonical_name(name: str, strip_suffixes: Iterable[str] = ()) -> str:
    """Normalise a genome name for cross-tool joins:
    strip a trailing FASTA extension, then any of the user-supplied
    `strip_suffixes` (e.g. '.contigs', '_genomic'), and lowercase.
    """
    n = (name or '').strip()
    low = n.lower()
    for ext in _FASTA_EXTS:
        if low.endswith(ext):
            n = n[: -len(ext)]
            low = n.lower()
            break
    for suf in strip_suffixes:
        if suf and low.endswith(suf.lower()):
            n = n[: -len(suf)]
            low = n.lower()
            break
    return low


def canonical_map(names: Iterable[str],
                  strip_suffixes: Iterable[str] = ()) -> Dict[str, str]:
    """{canonical_name: original_name} for a set of names."""
    return {canonical_name(n, strip_suffixes): n for n in names}


def match_report(bins: Iterable[str], source_keys: Iterable[str],
                 strip_suffixes: Iterable[str] = ()) -> Tuple[List[str], List[str]]:
    """Which bins have / lack a match in a source, by canonical name.
    Returns (matched_bins, unmatched_bins) using the ORIGINAL bin names."""
    src = set(canonical_name(k, strip_suffixes) for k in source_keys)
    matched, unmatched = [], []
    for b in bins:
        (matched if canonical_name(b, strip_suffixes) in src else unmatched).append(b)
    return matched, unmatched
