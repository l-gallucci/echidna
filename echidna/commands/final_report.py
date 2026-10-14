#!/usr/bin/env python3
"""
Echidna - final-report command

Join the per-MAG signals Echidna produces into one table: CheckM2
completeness/contamination (+tier), GTDB-Tk taxonomy, MIMAG rRNA/tRNA,
check-chimeras risk, classify-contigs per-contig split, and check-ssu-taxonomy
16S agreement. Every input is optional - whatever you pass is joined; missing
fields are 'NA'. A 'concerns' column flags the red signals for easy sorting.

This only reads TSVs already written by the other commands; it runs nothing and
needs no external tools.
"""

import sys
import argparse
from pathlib import Path
from typing import Dict, Set

from ..utils.checkm2 import parse_quality_report, find_quality_reports
from ..utils.chimera import parse_gtdbtk_summary, extract_taxonomy_level
from ..utils.rrna import parse_mimag_rna_summary
from ..utils.final_report import (
    norm_name, build_report, write_report,
    parse_chimera_report, parse_gunc_pass, parse_contig_taxonomy_bins, parse_ssu_report,
)
from ..utils.logging import log_info, log_warning, log_error, print_header, print_section


def _load_checkm2(report: Path, root: Path) -> Dict[str, Dict]:
    out: Dict[str, Dict] = {}
    paths = []
    if report:
        paths = [report]
    elif root:
        paths = find_quality_reports(root)
    for p in paths:
        for rec in parse_quality_report(p):
            out[norm_name(rec['name'])] = {
                'completeness': rec['completeness'],
                'contamination': rec['contamination'],
            }
    return out


def _load_taxonomy(gtdbtk_dir: Path, tax_level: str) -> Dict[str, Dict]:
    out: Dict[str, Dict] = {}
    for summary in sorted(gtdbtk_dir.rglob('gtdbtk.*.summary.tsv')):
        for name, rec in parse_gtdbtk_summary(summary).items():
            cls = rec.get('classification', '')
            out[norm_name(name)] = {
                'classification': cls,
                'rank_taxon': extract_taxonomy_level(cls, tax_level) if cls else 'NA',
            }
    return out


def _load_mimag(mimag_dir: Path) -> Dict[str, Dict]:
    out: Dict[str, Dict] = {}
    for summary in sorted(mimag_dir.rglob('mimag_rna_summary.tsv')):
        for name, rec in parse_mimag_rna_summary(summary).items():
            out[norm_name(name)] = rec
    return out


def register_parser(subparsers):
    parser = subparsers.add_parser(
        'final-report',
        help='Join per-MAG signals (CheckM2, GTDB-Tk, MIMAG, chimera, 16S) into one table',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Every input is optional; whatever you pass is joined by genome name (FASTA
extension stripped). Example with all of them:

  echidna final-report \\
    --checkm2-report checkm2/quality_report.tsv \\
    --gtdbtk-dir taxonomy/ \\
    --mimag-rna-dir mimag/ \\
    --chimera-report chimera_report.tsv \\
    --contig-taxonomy contig_taxonomy.tsv.bins.tsv \\
    --ssu-report ssu_taxonomy.tsv \\
    -o final_report.tsv

The 'concerns' column flags: contam>10, chimera_High/Medium, contig_chimera,
ssu_disagree, mimag_rna_incomplete. Empty concerns = no red flags (not a formal
pass; apply your own thresholds).
        """
    )
    parser.add_argument('--checkm2-report', type=Path,
                        help='A single CheckM2 quality_report.tsv')
    parser.add_argument('--checkm2-root', type=Path,
                        help='Directory to search for quality_report.tsv (recursive)')
    parser.add_argument('--gtdbtk-dir', type=Path, help='GTDB-Tk output directory')
    parser.add_argument('--mimag-rna-dir', type=Path,
                        help='run-mimag-rna output dir (mimag_rna_summary.tsv)')
    parser.add_argument('--chimera-report', type=Path,
                        help='check-chimeras chimera_report.tsv')
    parser.add_argument('--contig-taxonomy', type=Path,
                        help='classify-contigs <output>.bins.tsv')
    parser.add_argument('--ssu-report', type=Path,
                        help='check-ssu-taxonomy ssu_taxonomy.tsv')
    parser.add_argument('--genomes-dir', type=Path,
                        help='Optional: define the genome set from this directory '
                             '(else the union of names across the inputs)')
    parser.add_argument('--extensions', default='fa,fasta,fna',
                        help='Genome extensions for --genomes-dir (default: fa,fasta,fna)')
    parser.add_argument('--tax-level', default='phylum',
                        choices=['domain', 'phylum', 'class', 'order', 'family',
                                 'genus', 'species'],
                        help='GTDB rank shown in rank_taxon (default: phylum)')
    parser.add_argument('-o', '--output', type=Path, default=Path('final_report.tsv'),
                        help='Output TSV (default: final_report.tsv)')
    parser.set_defaults(func=run)
    return parser


def run(args):
    print_header("Echidna - Final Per-MAG Report")

    if args.checkm2_report and not args.checkm2_report.exists():
        log_error(f"CheckM2 report not found: {args.checkm2_report}")
        sys.exit(1)

    print_section("Loading sources")
    checkm2 = _load_checkm2(args.checkm2_report, args.checkm2_root) \
        if (args.checkm2_report or args.checkm2_root) else {}
    taxonomy = _load_taxonomy(args.gtdbtk_dir, args.tax_level) if args.gtdbtk_dir else {}
    mimag = _load_mimag(args.mimag_rna_dir) if args.mimag_rna_dir else {}
    chimera = parse_chimera_report(args.chimera_report) if args.chimera_report else {}
    gunc_pass = parse_gunc_pass(args.chimera_report) if args.chimera_report else {}
    contig = parse_contig_taxonomy_bins(args.contig_taxonomy) if args.contig_taxonomy else {}
    ssu = parse_ssu_report(args.ssu_report) if args.ssu_report else {}

    for label, m in (('CheckM2', checkm2), ('GTDB-Tk', taxonomy), ('MIMAG', mimag),
                     ('chimera', chimera), ('contig-taxonomy', contig), ('16S', ssu)):
        if m:
            log_info(f"  {label}: {len(m):,} genomes")

    # Universe of genomes
    universe: Set[str] = set()
    if args.genomes_dir:
        exts = [e.strip() for e in args.extensions.split(',')]
        for ext in exts:
            for p in args.genomes_dir.rglob(f"*.{ext}"):
                universe.add(norm_name(p.name))
    else:
        for m in (checkm2, taxonomy, mimag, chimera, contig, ssu):
            universe.update(m.keys())

    if not universe:
        log_error("No genomes found. Pass at least one source TSV or --genomes-dir.")
        sys.exit(1)

    print_section(f"Joining {len(universe):,} genomes")
    rows = build_report(universe, checkm2=checkm2, taxonomy=taxonomy, mimag=mimag,
                        chimera=chimera, contig=contig, ssu=ssu, gunc_pass=gunc_pass)
    write_report(rows, args.output)
    log_info(f"Final report written: {args.output}")

    print_section("CONCERNS SUMMARY")
    n_concern = sum(1 for r in rows if r['concerns'])
    print(f"  Genomes:              {len(rows):>6,}")
    print(f"  With concerns:        {n_concern:>6,}")
    if n_concern:
        from collections import Counter
        tally: Counter = Counter()
        for r in rows:
            for c in r['concerns'].split(';'):
                if c:
                    tally[c] += 1
        for flag, n in tally.most_common():
            print(f"    {flag:24s} {n:>6,}")

    log_info("\nDone.")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    register_parser(parser.add_subparsers(dest='command'))
    run(parser.parse_args())
