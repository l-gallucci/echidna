"""Report exact full-length contigs shared by MAGs, without refining the bins."""
import argparse
from pathlib import Path

from ..utils.bin_overlap import analyze_bins, discover_bins, write_reports
from ..utils.logging import log_info, log_warning, print_header


def register_parser(subparsers):
    parser = subparsers.add_parser(
        'check-shared-contigs',
        help='Report exact contigs shared by MAGs and conflicting FASTA IDs',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''Examples:
  echidna check-shared-contigs --bins-dir bins/ -o overlap/
  echidna check-shared-contigs --bins-dir Selected/ --recursive -o overlap/
  echidna check-shared-contigs --bins-dir metabat/ semibin/ -o overlap/

One FASTA per MAG; .fa/.fna/.fasta and their .gz/.gzip forms are recognized.
Matching ignores case, sequence whitespace and strand orientation, but requires
identical full-length sequences. Partial overlaps and ANI are not assessed.
Shared sequences and same-ID/different-sequence conflicts are reported separately.
No sequences are removed or reassigned; Binette refinement is a separate task.
''')
    parser.add_argument('--bins-dir', type=Path, nargs='+', required=True,
                        help='One or more directories containing genome FASTAs')
    parser.add_argument('-o', '--output-dir', type=Path, required=True,
                        help='New or empty output directory for TSV reports and summary.json')
    parser.add_argument('--recursive', action='store_true',
                        help='Include nested folders such as Selected/HQ and Selected/MQ')
    parser.add_argument('--extensions', default='fa,fna,fasta',
                        help='Comma-separated FASTA extensions (default: fa,fna,fasta; gzip auto-detected)')
    parser.set_defaults(func=run)
    return parser


def run(args):
    print_header('Echidna - Shared contigs between MAGs')
    if args.output_dir.exists() and (not args.output_dir.is_dir() or any(args.output_dir.iterdir())):
        raise ValueError(f'Output must be a new or empty directory: {args.output_dir}')
    files = discover_bins(args.bins_dir, args.recursive, args.extensions)
    log_info(f'Checking {len(files):,} MAG FASTAs')
    result = analyze_bins(files)
    summary = write_reports(result, args.output_dir)
    log_info(f"Shared full-length sequences: {summary['n_shared_sequences']:,}")
    log_info(f"MAGs with shared contigs: {summary['n_bins_with_shared_contigs']:,}")
    log_info(f"MAG pairs with shared contigs: {summary['n_pairs_with_shared_contigs']:,}")
    if summary['n_conflicting_ids']:
        log_warning(f"Same-ID/different-sequence groups: {summary['n_conflicting_ids']:,}; see id_conflicts.tsv")
    log_info(f'Reports: {args.output_dir}')
