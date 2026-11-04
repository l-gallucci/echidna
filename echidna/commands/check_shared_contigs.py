"""Report or resolve exact whole-contig sharing, preserving original MAG FASTAs."""
import argparse
from pathlib import Path

from ..utils.bin_overlap import analyze_bins, discover_bins, write_reports
from ..utils.logging import log_info, log_warning, print_header
from ..utils.bin_resolution import Inventory, load_quality, load_clusters, plan_resolution, write_resolution
from ..utils.bin_taxonomy import BinTaxonomy, TaxonomyEvidence, validate_thresholds
from ..utils.contig_taxonomy import VALID_RANKS


def add_taxonomy_options(parser):
    parser.add_argument('--contig-taxonomy', type=Path,
                        help='Sequence-linked classify-contigs per-contig TSV; gate MAGs and replacements '
                             'by taxonomy, not the per-bin .bins.tsv')
    parser.add_argument('--taxonomy-rank', choices=VALID_RANKS,
                        help='Require this rank in the reused table (default: use table rank)')
    parser.add_argument('--min-taxonomy-resolved-frac', type=float, default=0.5,
                        help='Minimum fraction of all MAG bases resolved at the rank (default: 0.5)')
    parser.add_argument('--min-foreign-frac', type=float, default=0.05,
                        help='Reject when minority taxa reach this fraction of resolved bases (default: 0.05)')


def register_parser(subparsers):
    parser = subparsers.add_parser(
        'check-shared-contigs',
        help='Report or resolve shared MAG contigs and conflicting FASTA IDs',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''Examples:
  echidna check-shared-contigs --bins-dir bins/ -o overlap/
  echidna check-shared-contigs --bins-dir Selected/ --recursive -o overlap/
  echidna check-shared-contigs --bins-dir metabat/ semibin/ -o overlap/

One FASTA per MAG; .fa/.fna/.fasta and their .gz/.gzip forms are recognized.
Matching ignores case, sequence whitespace and strand orientation, but requires
identical full-length sequences. Partial overlaps and ANI are not assessed.
Shared sequences and same-ID/different-sequence conflicts are reported separately.
Default: report only. Optional --exclude-mags or --remove-shared-contigs exports
resolved copies, ranked by --quality-report; originals are never changed.
Cluster replacements require --candidates-dir plus --drep-dir or --cluster-table.
--contig-taxonomy reuses sequence-linked classify-contigs assignments to gate all
retained MAGs and candidates by taxonomic coverage and foreign-base fraction.
''')
    parser.add_argument('--bins-dir', type=Path, nargs='+', required=True,
                        help='One or more directories containing genome FASTAs')
    parser.add_argument('-o', '--output-dir', type=Path, required=True,
                        help='New or empty output directory for TSV reports and summary.json')
    parser.add_argument('--recursive', action='store_true',
                        help='Include nested folders such as Selected/HQ and Selected/MQ')
    parser.add_argument('--extensions', default='fa,fna,fasta',
                        help='Comma-separated FASTA extensions (default: fa,fna,fasta; gzip auto-detected)')
    resolution = parser.add_mutually_exclusive_group()
    resolution.add_argument('--exclude-mags', action='store_true',
                            help='Keep priority MAGs; export conflicting MAGs separately')
    resolution.add_argument('--remove-shared-contigs', action='store_true',
                            help='Keep shared contigs in priority MAGs and trim lower-priority copies')
    parser.add_argument('--quality-report', type=Path,
                        help='CheckM2 TSV or CSV/TSV with genome, completeness and contamination')
    parser.add_argument('--rank-by', choices=['completeness', 'quality-score'], default=None,
                        help='Priority: completeness then contamination (default), or completeness - 5*contamination')
    clusters = parser.add_mutually_exclusive_group()
    clusters.add_argument('--drep-dir', type=Path,
                          help='dRep work directory; read data_tables/Cdb.csv secondary clusters')
    clusters.add_argument('--cluster-table', type=Path,
                          help='CSV/TSV with genome and cluster (or secondary_cluster)')
    parser.add_argument('--candidates-dir', type=Path, nargs='+',
                        help='Original candidate MAG FASTAs for same-cluster replacements')
    parser.add_argument('--require-retained', action='store_true',
                        help='After writing resolution reports, return an error if no MAG passed selection')
    add_taxonomy_options(parser)
    parser.set_defaults(func=run)
    return parser


def run(args):
    print_header('Echidna - Shared contigs between MAGs')
    if args.output_dir.exists() and (not args.output_dir.is_dir() or any(args.output_dir.iterdir())):
        raise ValueError(f'Output must be a new or empty directory: {args.output_dir}')
    mode = 'exclude-mags' if args.exclude_mags else ('remove-contigs' if args.remove_shared_contigs else None)
    cluster_path = args.cluster_table or (args.drep_dir / 'data_tables/Cdb.csv' if args.drep_dir else None)
    taxonomy_path = getattr(args, 'contig_taxonomy', None)
    if not mode and any((args.quality_report, args.rank_by, cluster_path, args.candidates_dir, taxonomy_path)):
        raise ValueError('Quality and replacement options require --exclude-mags or --remove-shared-contigs')
    if mode and not args.quality_report:
        raise ValueError('Overlap resolution requires --quality-report; quality is not inferred from FASTA length')
    if getattr(args, 'require_retained', False) and not mode:
        raise ValueError('--require-retained requires a resolution mode')
    if bool(args.candidates_dir) != bool(cluster_path):
        raise ValueError('Replacements require both --candidates-dir and --drep-dir or --cluster-table')
    validate_thresholds(getattr(args, 'min_taxonomy_resolved_frac', 0.5), getattr(args, 'min_foreign_frac', 0.05))
    if not taxonomy_path and (getattr(args, 'taxonomy_rank', None) or
                             getattr(args, 'min_taxonomy_resolved_frac', 0.5) != 0.5 or
                             getattr(args, 'min_foreign_frac', 0.05) != 0.05):
        raise ValueError('Taxonomy rank/threshold options require --contig-taxonomy')
    files = discover_bins(args.bins_dir, args.recursive, args.extensions)
    log_info(f'Checking {len(files):,} MAG FASTAs')
    result = analyze_bins(files)
    if mode:
        inventory = Inventory(files, result, args.candidates_dir, args.recursive, args.extensions)
        clusters, eligible, unavailable = load_clusters(cluster_path, inventory) if cluster_path else ({}, set(), [])
        quality = load_quality(args.quality_report, inventory)
        taxonomy = BinTaxonomy(TaxonomyEvidence(
            taxonomy_path, getattr(args, 'taxonomy_rank', None),
            getattr(args, 'min_taxonomy_resolved_frac', 0.5), getattr(args, 'min_foreign_frac', 0.05)), inventory
        ) if taxonomy_path else None
        resolution = plan_resolution(inventory, quality, mode, args.rank_by or 'completeness', clusters, eligible, taxonomy)
        metadata = {'quality_report': str(args.quality_report.resolve()),
                    'cluster_table': str(cluster_path.resolve()) if cluster_path else None,
                    'candidate_directories': [str(path.resolve()) for path in (args.candidates_dir or [])]}
        metadata['contig_taxonomy'] = str(taxonomy_path.resolve()) if taxonomy_path else None
        summary = write_resolution(result, inventory, quality, resolution, args.output_dir, clusters, unavailable, metadata)
        selected = summary['resolution']
        log_info(f"Final retained MAGs: {selected['n_retained_mags']:,}; accepted replacements: {selected['n_replacements']:,}")
        log_info('Global overlap recheck on exported FASTAs: passed (0 shared full-length sequences)')
        if taxonomy and selected['n_retained_mags']:
            log_info(f'Final taxonomy recheck: passed for {selected["taxonomy_retained_mags_checked"]} MAGs '
                     f'at rank {selected["taxonomy_rank"]}; see taxonomy_checks.tsv')
        if not selected['n_retained_mags']:
            log_warning('No MAGs passed selection; retained/ is empty. Inspect resolution.tsv and taxonomy_checks.tsv')
            if getattr(args, 'require_retained', False):
                raise ValueError(f'No MAGs passed selection; complete reports saved in {args.output_dir}')
        if selected['quality_reassessment_required']:
            log_warning('Trimmed MAG quality must be reassessed with CheckM2; final quality fields are blank in resolution.tsv')
        if selected['n_replacements']:
            log_warning('Replacements need their own downstream taxonomy/chimera/RNA checks; do not reuse the old representative reports')
        if selected['n_unresolved_clusters']:
            log_warning('Some clusters have no compatible replacement; see unresolved_clusters.tsv')
        if selected['n_unavailable_cluster_members']:
            log_warning('Some cluster members have no supplied FASTA; see unavailable_cluster_members.tsv')
    else:
        summary = write_reports(result, args.output_dir)
    prefix = 'Input ' if mode else ''
    log_info(f"{prefix}shared full-length sequences: {summary['n_shared_sequences']:,}")
    log_info(f"{prefix}MAGs with shared contigs: {summary['n_bins_with_shared_contigs']:,}")
    log_info(f"{prefix}MAG pairs with shared contigs: {summary['n_pairs_with_shared_contigs']:,}")
    if summary['n_conflicting_ids']:
        log_warning(f"Same-ID/different-sequence groups: {summary['n_conflicting_ids']:,}; see id_conflicts.tsv")
    log_info(f'Reports: {args.output_dir}')
