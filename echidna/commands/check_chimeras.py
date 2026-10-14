#!/usr/bin/env python3
"""
Echidna - check-chimeras command
Detect chimeric MAGs and circular contigs using GC heterogeneity and GTDB-Tk taxonomy signals.

Chimeric MAGs ("puzzle MAGs") are a common artefact in metagenomic binning where contigs
from multiple organisms are incorrectly grouped into a single bin. Circular contigs can
also be chimeric if they were assembled from reads belonging to two different organisms
that share similar coverage and k-mer profiles.

Detection strategy
------------------
1. **GC heterogeneity** (all bins):
   - Bins whose contigs show high GC variance (delta > 10 %) are flagged.

2. **Windowed GC analysis** (single-contig circular candidates ≥ 100 kbp):
   - A sliding window scans the contig for abrupt GC shifts that may indicate
     a chimeric junction (assembly artefact).

3. **GTDB-Tk taxonomy confidence** (optional):
   - Bins with low MSA recovery (< 10 %), classification warnings, or no
     FastANI reference match are flagged as unreliably placed.

4. **CheckM2 contamination** (optional):
   - Contamination > 5 % is a strong independent chimera indicator.
"""

import csv
import sys
import argparse
from pathlib import Path
from typing import Dict, List, Optional

from ..utils.chimera import (
    analyze_bin_gc,
    strong_gc_outlier,
    parse_gtdbtk_summary,
    extract_taxonomy_level,
    has_taxonomy_warning,
    assess_chimera_verdict,
)
from ..utils.checkm2 import parse_quality_report
from ..utils.names import canonical_name, match_report
from ..utils.validation import validate_dir_exists
from ..utils.logging import log_info, log_warning, log_error, print_header, print_section


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _find_gtdbtk_summaries(gtdbtk_dir: Path) -> List[Path]:
    """Locate GTDB-Tk summary TSV files under gtdbtk_dir."""
    candidates = [
        gtdbtk_dir / 'gtdbtk.bac120.summary.tsv',
        gtdbtk_dir / 'gtdbtk.ar53.summary.tsv',
        gtdbtk_dir / 'classify' / 'gtdbtk.bac120.summary.tsv',
        gtdbtk_dir / 'classify' / 'gtdbtk.ar53.summary.tsv',
    ]
    found = [p for p in candidates if p.exists()]
    if not found:
        found = sorted(gtdbtk_dir.rglob('gtdbtk.*.summary.tsv'))
    return found


def _load_taxonomy(gtdbtk_dir: Path) -> Dict[str, Dict]:
    summaries = _find_gtdbtk_summaries(gtdbtk_dir)
    if not summaries:
        log_warning(f"No GTDB-Tk summary files found under {gtdbtk_dir}")
        return {}
    taxonomy = {}
    for summary in summaries:
        log_info(f"  Loading taxonomy: {summary.name}")
        taxonomy.update(parse_gtdbtk_summary(summary))
    return taxonomy


def _load_contamination(checkm2_file: Path) -> Dict[str, float]:
    records = parse_quality_report(checkm2_file)
    return {rec['name']: rec['contamination'] for rec in records}


# ---------------------------------------------------------------------------
# Core logic
# ---------------------------------------------------------------------------

def _load_contig_foreign_frac(path: Path) -> Dict[str, float]:
    """classify-contigs <out>.bins.tsv -> {bin: foreign_frac}."""
    out: Dict[str, float] = {}
    with open(path, newline='') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            b = (row.get('bin') or '').strip()
            if not b:
                continue
            try:
                out[b] = float(row.get('foreign_frac') or 0.0)
            except ValueError:
                out[b] = 0.0
    return out


def _warn_unmatched(bin_names: List[str], source_keys, label: str,
                    strip: List[str]) -> None:
    """Warn (not fail) about bins with no match in a source, by canonical name."""
    _, unmatched = match_report(bin_names, source_keys, strip)
    if unmatched:
        shown = ', '.join(unmatched[:8]) + ('...' if len(unmatched) > 8 else '')
        log_warning(f"{len(unmatched):,}/{len(bin_names):,} bins have NO {label} "
                    f"match (signal will be N/A for them): {shown}")


def run_chimera_check(
    bins_dir: Path,
    output_file: Path,
    extensions: List[str],
    gtdbtk_dir: Optional[Path],
    checkm2_file: Optional[Path],
    gunc_dir: Optional[Path],
    contig_taxonomy_file: Optional[Path],
    min_foreign_frac: float,
    gc_min_length: int,
    gc_max_dev: float,
    taxonomy_level: str,
    name_strip: List[str],
    dry_run: bool,
) -> Dict[str, int]:
    """
    Decide a chimera verdict per bin from TWO signals: GUNC pass/fail and the
    per-contig taxonomy (MMseqs) foreign fraction. Both chimeric -> Chimeric;
    one -> Check; neither -> Clean. CheckM2 contamination, GTDB-Tk taxonomy and a
    strong GC outlier are REPORTED as columns but do not drive the verdict.
    Every source is joined by a canonical genome name (utils.names) and unmatched
    bins are warned about.
    """
    print_section("Loading auxiliary data")

    taxonomy: Dict[str, Dict] = _load_taxonomy(gtdbtk_dir) if gtdbtk_dir else {}
    if taxonomy:
        log_info(f"Taxonomy loaded for {len(taxonomy):,} genomes")

    contamination_map: Dict[str, float] = (
        _load_contamination(checkm2_file) if checkm2_file else {})
    if contamination_map:
        log_info(f"CheckM2 contamination loaded for {len(contamination_map):,} genomes")

    gunc_map: Dict[str, Dict] = {}
    if gunc_dir:
        from ..utils.chimera import parse_gunc_output
        gunc_map = parse_gunc_output(gunc_dir)
        if gunc_map:
            log_info(f"GUNC results loaded for {len(gunc_map):,} genomes")

    contig_frac_map: Dict[str, float] = (
        _load_contig_foreign_frac(contig_taxonomy_file) if contig_taxonomy_file else {})
    if contig_frac_map:
        log_info(f"Per-contig taxonomy (MMseqs) loaded for {len(contig_frac_map):,} bins")

    # Canonicalise all source keys once so the per-bin join is naming-robust.
    def _canon_map(m):
        return {canonical_name(k, name_strip): v for k, v in m.items()}
    tax_c = _canon_map(taxonomy)
    contam_c = _canon_map(contamination_map)
    gunc_c = _canon_map(gunc_map)
    frac_c = _canon_map(contig_frac_map)

    print_section("Scanning bin files")
    bin_files: List[Path] = []
    for ext in extensions:
        bin_files.extend(bins_dir.rglob(f"*.{ext}"))
    bin_files = sorted(set(bin_files))
    if not bin_files:
        log_error(f"No bin files found in {bins_dir} with extensions: {extensions}")
        return {'verdict': {}, 'gunc': {'pass': 0, 'fail': 0, 'na': 0}, 'gunc_used': bool(gunc_dir)}
    log_info(f"Found {len(bin_files):,} bin file(s)")

    bin_names = [b.stem for b in bin_files]
    if taxonomy:
        _warn_unmatched(bin_names, taxonomy.keys(), "GTDB-Tk", name_strip)
    if contamination_map:
        _warn_unmatched(bin_names, contamination_map.keys(), "CheckM2", name_strip)
    if gunc_map:
        _warn_unmatched(bin_names, gunc_map.keys(), "GUNC", name_strip)
    if contig_frac_map:
        _warn_unmatched(bin_names, contig_frac_map.keys(), "contig-taxonomy", name_strip)

    print_section("Deciding chimera verdict per bin")
    rows: List[Dict] = []
    verdict_counts: Dict[str, int] = {'Chimeric': 0, 'Check': 0, 'Clean': 0, 'Unknown': 0}
    gunc_counts: Dict[str, int] = {'pass': 0, 'fail': 0, 'na': 0}

    for bin_file in bin_files:
        bin_name = bin_file.stem
        key = canonical_name(bin_name, name_strip)
        log_info(f"\n  [{bin_name}]")

        try:
            gc_data = analyze_bin_gc(bin_file)
        except Exception as e:
            log_warning(f"  Could not analyse {bin_file.name}: {e}")
            continue

        bin_stats = gc_data['bin_gc_stats']
        # GC: informational only - a long contig far from the bin mean (short
        # contigs ignored). Does NOT affect the verdict.
        gc_out = strong_gc_outlier(gc_data['per_contig_gc'],
                                   gc_data['per_contig_len'],
                                   min_len=gc_min_length, max_dev=gc_max_dev)
        gc_out_str = 'N/A'
        if gc_out:
            gc_out_str = (f"{gc_out['contig']} {gc_out['gc'] * 100:.1f}% "
                          f"{gc_out['length']:,}bp vs bin_mean {gc_out['bin_mean'] * 100:.1f}%")

        tax_record = tax_c.get(key, {})
        classification = tax_record.get('classification', '')
        taxon = extract_taxonomy_level(classification, taxonomy_level) if classification else 'N/A'
        tax_warn = has_taxonomy_warning(tax_record) if tax_record else False

        contamination = contam_c.get(key)

        gunc_record = gunc_c.get(key, {})
        gunc_fail = (not gunc_record['pass_gunc']) if gunc_record else None
        gunc_css = gunc_record.get('css') if gunc_record else None
        if gunc_dir:
            if not gunc_record:
                gunc_counts['na'] += 1
            elif gunc_record['pass_gunc']:
                gunc_counts['pass'] += 1
            else:
                gunc_counts['fail'] += 1

        # Per-contig taxonomy signal: substantial foreign fraction = chimeric.
        contig_foreign_frac = frac_c.get(key)
        tax_chimeric = (None if contig_foreign_frac is None
                        else contig_foreign_frac >= min_foreign_frac)

        verdict, reasons = assess_chimera_verdict(
            gunc_fail=gunc_fail, tax_chimeric=tax_chimeric,
            gunc_css=gunc_css, contig_foreign_frac=contig_foreign_frac)
        verdict_counts[verdict] = verdict_counts.get(verdict, 0) + 1

        log_info(
            f"  Contigs={gc_data['n_contigs']}  "
            f"GUNC={'N/A' if gunc_fail is None else ('FAIL' if gunc_fail else 'pass')}  "
            f"foreign_frac={'N/A' if contig_foreign_frac is None else f'{contig_foreign_frac * 100:.1f}%'}  "
            f"Verdict={verdict}"
        )
        for reason in reasons:
            log_info(f"  ⚠  {reason}")
        if gc_out:
            log_info(f"  (info) strong GC outlier: {gc_out_str}")

        rows.append({
            'bin': bin_name,
            'n_contigs': gc_data['n_contigs'],
            'gc_mean_pct': f"{(bin_stats.get('mean') or float('nan')) * 100:.2f}",
            'gc_strong_outlier': gc_out_str,
            'checkm2_contamination': (
                f"{contamination:.2f}" if contamination is not None else 'N/A'
            ),
            'gtdbtk_taxonomy': taxon,
            'gtdbtk_warning': tax_warn,
            'gunc_pass': ('N/A' if not gunc_record else gunc_record['pass_gunc']),
            'gunc_css': ('N/A' if gunc_css is None else f"{gunc_css:.3f}"),
            'contig_foreign_frac': ('N/A' if contig_foreign_frac is None
                                    else f"{contig_foreign_frac:.4f}"),
            'chimera_verdict': verdict,
            'reasons': ' | '.join(reasons),
        })

    if not dry_run and rows:
        fieldnames = [
            'bin', 'n_contigs', 'gc_mean_pct', 'gc_strong_outlier',
            'checkm2_contamination', 'gtdbtk_taxonomy', 'gtdbtk_warning',
            'gunc_pass', 'gunc_css', 'contig_foreign_frac',
            'chimera_verdict', 'reasons',
        ]
        with open(output_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()
            writer.writerows(rows)
        log_info(f"\nChimera report written to: {output_file}")

    return {'verdict': verdict_counts, 'gunc': gunc_counts, 'gunc_used': bool(gunc_dir)}


# ---------------------------------------------------------------------------
# CLI registration
# ---------------------------------------------------------------------------

def register_parser(subparsers):
    """Register the check-chimeras subcommand."""
    parser = subparsers.add_parser(
        'check-chimeras',
        help='Detect chimeric MAGs and circular contigs using GC and taxonomy signals',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Verdict from TWO signals only (nothing else enters it)
---------------------
  GUNC gene-level     : pass.GUNC = False  -> GUNC says chimeric
  Per-contig taxonomy : foreign LENGTH fraction >= --min-foreign-frac (MMseqs via
                        --contig-taxonomy / classify-contigs) -> taxonomy chimeric

  both chimeric       -> Chimeric
  only one            -> Check   (review)
  neither             -> Clean
  (with a signal missing, the verdict uses the one available: one firing -> Check)

Reported but NOT in the verdict: CheckM2 contamination, GTDB-Tk taxonomy and a
strong GC outlier (a long contig far from the bin's mean GC; short contigs are
ignored) - all go to their own columns for inspection.

Examples
--------
  # Verdict from GUNC + per-contig taxonomy (the two that matter)
  echidna check-chimeras \\
    --bins-dir bins/ --output chimera_report.tsv \\
    --gunc-dir gunc_results/ \\
    --contig-taxonomy contig_taxonomy.tsv.bins.tsv

  # Add CheckM2 / GTDB-Tk as informational columns
  echidna check-chimeras --bins-dir bins/ --gunc-dir gunc_results/ \\
    --contig-taxonomy contig_taxonomy.tsv.bins.tsv \\
    --checkm2-report checkm2/quality_report.tsv --gtdbtk-dir gtdbtk_results/

  # Tolerate a naming suffix some tools add
  echidna check-chimeras --bins-dir bins/ --gunc-dir gunc/ --name-strip .contigs
        """
    )

    parser.add_argument('--bins-dir', type=Path, required=True,
                        help='Directory containing bin/MAG FASTA files')
    parser.add_argument('-o', '--output', type=Path, default=Path('chimera_report.tsv'),
                        help='Output TSV report (default: chimera_report.tsv)')

    # Verdict signals
    parser.add_argument('--gunc-dir', type=Path,
                        help='GUNC output directory (pass.GUNC drives the verdict)')
    parser.add_argument('--contig-taxonomy', type=Path,
                        help='classify-contigs <output>.bins.tsv (per-contig taxonomy; '
                             'foreign-fraction drives the verdict)')
    parser.add_argument('--min-foreign-frac', type=float, default=0.05,
                        help='Foreign LENGTH fraction at/above which the per-contig '
                             'taxonomy signal counts as chimeric (default: 0.05)')

    # Informational-only sources (reported, not in the verdict)
    parser.add_argument('--gtdbtk-dir', type=Path,
                        help='GTDB-Tk output directory (reported column only)')
    parser.add_argument('--checkm2-report', type=Path,
                        help='CheckM2 quality_report.tsv (reported column only)')
    parser.add_argument('--name-strip', default='',
                        help='Comma-separated suffixes to strip from genome names when '
                             'joining sources (e.g. ".contigs,_genomic"); FASTA extensions '
                             'are always stripped')

    # GC outlier (reported column only)
    parser.add_argument('--gc-min-length', type=int, default=5000,
                        help='Ignore contigs shorter than this for the GC-outlier column '
                             '(default: 5000 bp)')
    parser.add_argument('--gc-max-dev', type=float, default=0.10,
                        help='Report a GC outlier when a long contig deviates from the '
                             "bin's mean GC by more than this (default: 0.10 = 10 points)")

    parser.add_argument('--taxonomy-level', default='phylum',
                        choices=['domain', 'phylum', 'class', 'order', 'family',
                                 'genus', 'species'],
                        help='Taxonomic level to report in output (default: phylum)')
    parser.add_argument('--extensions', default='fa,fasta,fna',
                        help='Comma-separated bin file extensions (default: fa,fasta,fna)')
    parser.add_argument('--dry-run', action='store_true',
                        help='Print results to terminal without writing output file')

    parser.set_defaults(func=run)
    return parser


def run(args):
    """Execute the check-chimeras command."""
    validate_dir_exists(args.bins_dir, "Bins directory")

    if args.gtdbtk_dir:
        validate_dir_exists(args.gtdbtk_dir, "GTDB-Tk directory")

    if args.checkm2_report and not args.checkm2_report.exists():
        log_error(f"CheckM2 report not found: {args.checkm2_report}")
        sys.exit(1)

    if args.gunc_dir:
        validate_dir_exists(args.gunc_dir, "GUNC directory")

    if args.contig_taxonomy and not args.contig_taxonomy.exists():
        log_error(f"contig-taxonomy file not found: {args.contig_taxonomy}")
        sys.exit(1)

    if not args.gunc_dir and not args.contig_taxonomy:
        log_warning("Neither --gunc-dir nor --contig-taxonomy given; the verdict "
                    "needs at least one of them (all bins will be 'Unknown').")

    extensions = [e.strip() for e in args.extensions.split(',')]
    name_strip = [s.strip() for s in args.name_strip.split(',') if s.strip()]

    print_header("Echidna - Chimera Check")

    result = run_chimera_check(
        bins_dir=args.bins_dir,
        output_file=args.output,
        extensions=extensions,
        gtdbtk_dir=args.gtdbtk_dir,
        checkm2_file=args.checkm2_report,
        gunc_dir=args.gunc_dir,
        contig_taxonomy_file=args.contig_taxonomy,
        min_foreign_frac=args.min_foreign_frac,
        gc_min_length=args.gc_min_length,
        gc_max_dev=args.gc_max_dev,
        taxonomy_level=args.taxonomy_level,
        name_strip=name_strip,
        dry_run=args.dry_run,
    )

    verdict_counts = result['verdict']
    if not verdict_counts:
        log_error("No bins were analysed — check --bins-dir and --extensions")
        sys.exit(1)

    print_section("CHIMERA VERDICT SUMMARY")
    total = sum(verdict_counts.values())
    for level in ('Chimeric', 'Check', 'Clean', 'Unknown'):
        count = verdict_counts.get(level, 0)
        marker = '⚠ ' if level in ('Chimeric', 'Check') else '  '
        print(f"  {marker}{level:8s}: {count:>6,}")

    flagged = verdict_counts.get('Chimeric', 0) + verdict_counts.get('Check', 0)
    if total > 0:
        pct = 100 * flagged / total
        print(f"\n  Flagged (Chimeric + Check): {flagged:,} / {total:,} ({pct:.1f}%)")

    # GUNC pass/fail tally (only when --gunc-dir was supplied)
    if result['gunc_used']:
        g = result['gunc']
        evaluated = g['pass'] + g['fail']
        print(f"\n  GUNC pass.GUNC=True:     {g['pass']:>6,}")
        print(f"  GUNC chimeras (False):   {g['fail']:>6,}")
        if g['na']:
            print(f"  GUNC no result:          {g['na']:>6,}")
        if evaluated:
            pct = 100 * g['fail'] / evaluated
            print(f"  GUNC chimera rate:       {pct:.1f}% of {evaluated:,} evaluated")

    if args.dry_run:
        log_warning("\nDRY-RUN mode: no output file written")

    log_info("\n✓ Done!")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    register_parser(parser.add_subparsers(dest='command'))
    run(parser.parse_args())
