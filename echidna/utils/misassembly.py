"""
Misassembly data parsing utilities (anvi'o clipping files)
"""

import csv
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Set, Tuple

from .logging import log_info, log_warning


# A clipping record: (position, clipping_coverage)
ClipRecord = Tuple[int, int]


def find_clipping_files(misassembly_dir: Path) -> List[Path]:
    """Find all *-clipping.txt files in a directory."""
    files = sorted(misassembly_dir.glob('*-clipping.txt'))
    if not files:
        log_warning(f"No *-clipping.txt files found in {misassembly_dir}")
    return files


def parse_clipping_files(clipping_files: List[Path]) -> Dict[str, List[ClipRecord]]:
    """Read native anvi'o contig/pos/clipping or legacy Echidna columns.

    Header-only tables mean no events. Invalid schemas/values are errors,
    never a negative biological result. Positions are zero-based.
    """
    data = defaultdict(dict)
    for path in clipping_files:
        with open(path) as f:
            reader = csv.DictReader(f, delimiter='\t')
            fields = set(reader.fieldnames or [])
            native = {'contig', 'pos', 'clipping'}.issubset(fields)
            legacy = {'contig_name', 'position'}.issubset(fields)
            if not (native or legacy):
                raise ValueError(f'Unrecognised clipping columns in {path}')
            for lineno, row in enumerate(reader, 2):
                try:
                    name = row['contig' if native else 'contig_name'].strip()
                    pos = int(row['pos' if native else 'position'])
                    cov = int(row['clipping'] if native else row.get('clipping_coverage', '1'))
                    if not name or pos < 0 or cov < 0:
                        raise ValueError('negative value or empty contig')
                except (ValueError, TypeError, AttributeError) as e:
                    raise ValueError(f'Invalid clipping row {lineno} in {path}: {e}') from e
                # Repeated samples at the same position must not create repeated cuts.
                data[name][pos] = max(cov, data[name].get(pos, 0))
    return {name: sorted(records.items()) for name, records in data.items()}


def filter_by_coverage(
    clipping_data: Dict[str, List[ClipRecord]],
    min_coverage: int,
) -> Dict[str, List[ClipRecord]]:
    """
    Return a filtered copy keeping only clipping records where coverage >= min_coverage.
    Contigs with no records surviving the filter are dropped entirely.
    """
    if min_coverage <= 1:
        return clipping_data

    filtered: Dict[str, List[ClipRecord]] = {}
    for contig, records in clipping_data.items():
        passing = [r for r in records if r[1] >= min_coverage]
        if passing:
            filtered[contig] = passing

    n_removed = sum(len(v) for v in clipping_data.values()) - \
                sum(len(v) for v in filtered.values())
    log_info(
        f"Coverage filter (≥{min_coverage}): "
        f"kept {sum(len(v) for v in filtered.values()):,} positions, "
        f"dropped {n_removed:,} low-coverage positions"
    )
    return filtered


def load_clipping_data(
    misassembly_dir: Path,
    min_coverage: int = 0,
) -> Dict[str, List[ClipRecord]]:
    """
    Load all clipping data from a directory, optionally filtered by coverage.

    Returns dict mapping contig_name -> sorted list of (position, coverage).
    """
    files = find_clipping_files(misassembly_dir)
    if not files:
        raise ValueError(f"No *-clipping.txt files in {misassembly_dir}")

    data = parse_clipping_files(files)
    if min_coverage > 1:
        data = filter_by_coverage(data, min_coverage)

    n_contigs = len(data)
    n_positions = sum(len(v) for v in data.values())
    log_info(f"Clipping data: {n_contigs:,} contigs, {n_positions:,} positions")
    return data


def load_clipping_positions(
    misassembly_dir: Path,
    min_coverage: int = 0,
) -> Dict[str, List[int]]:
    """
    Load clipping positions (position-only, no coverage) from a directory.
    Backward-compatible wrapper around load_clipping_data.
    """
    data = load_clipping_data(misassembly_dir, min_coverage)
    return {contig: [pos for pos, _ in records] for contig, records in data.items()}


def get_contigs_with_misassemblies(
    misassembly_dir: Path,
    min_coverage: int = 0,
) -> Set[str]:
    """Return the set of contig names that have at least one clipping position."""
    return set(load_clipping_positions(misassembly_dir, min_coverage).keys())
