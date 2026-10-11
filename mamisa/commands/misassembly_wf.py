#!/usr/bin/env python3
"""
MaMISA - misassembly-wf command

Chains the misassembly-cleaning branch end to end (like `gtdbtk classify_wf`
chains its steps), calling each sub-command in the right conda env:

  anvi      run-anvi-misassemblies   (anvi'o env)  -> *-clipping.txt
  process   process-large-contigs    (CheckM2 via --checkm2-env)
  readchim  check-read-chimeras      (optional; needs Kraken2/Kaiju taxonomy)
  classify  classify-clipping
  filter    filter-misassemblies     -> assembly_clean.fa

PacBio HiFi only (the soft-clipping method needs long reads, Trigodet 2025).

Resume/partial runs: --from <step> starts at a step, --skip <step[,step]> drops
steps, --dry-run prints the commands without running them.
"""

import argparse
from pathlib import Path
from typing import List

from ..utils.workflow_run import Step, execute_plan
from ..utils.logging import log_error, print_header

ALL_KEYS = ['anvi', 'process', 'readchim', 'classify', 'filter']


def build_plan(args) -> List[Step]:
    """Build the ordered Step list from parsed args. Pure (no side effects)."""
    out = args.output_dir
    misasm_dir = out / '01_misassemblies'
    anvi_prefix = misasm_dir / args.anvi_prefix
    plc_dir = out / '02_large_contigs'
    chimera_dir = out / '03_chimera'
    classified = out / '04_classified' / 'clipping_classified.tsv'
    clean = out / 'assembly_clean.fa'

    use_readchim = bool(args.kraken2_output or args.kaiju_output)

    steps: List[Step] = []

    steps.append(Step(
        'anvi', 'anvi-script-find-misassemblies', args.anvio_env,
        ['mamisa', 'run-anvi-misassemblies',
         '-b', str(args.bam), '-o', str(anvi_prefix), '-T', str(args.threads)],
    ))

    steps.append(Step(
        'process', 'process-large-contigs', None,
        ['mamisa', 'process-large-contigs',
         '-a', str(args.assembly), '-m', str(misasm_dir), '-o', str(plc_dir),
         '--threads', str(args.threads), '--checkm2-env', args.checkm2_env,
         '--max-length', str(args.max_length),
         '--min-completeness', str(args.min_completeness),
         '--max-contamination', str(args.max_contamination)],
    ))

    if use_readchim:
        rc = ['mamisa', 'check-read-chimeras',
              '--bam', str(args.bam), '-o', str(chimera_dir)]
        if args.kaiju_output:
            rc += ['--kaiju-output', str(args.kaiju_output)]
            if args.kaiju_names:
                rc += ['--kaiju-names', str(args.kaiju_names)]
        else:
            rc += ['--kraken2-output', str(args.kraken2_output)]
            if args.kraken2_report:
                rc += ['--kraken2-report', str(args.kraken2_report)]
        steps.append(Step('readchim', 'check-read-chimeras', None, rc))

    cl = ['mamisa', 'classify-clipping',
          '--bam', str(args.bam), '-m', str(misasm_dir), '-o', str(classified)]
    if use_readchim:
        cl += ['--taxonomy-windows', str(chimera_dir / 'chimera_read_windows.tsv')]
    steps.append(Step('classify', 'classify-clipping', None, cl))

    fm = ['mamisa', 'filter-misassemblies',
          '-a', str(plc_dir / 'assembly_for_filtering.fa'),
          '-m', str(misasm_dir),
          '-g', str(plc_dir / 'HQ_extracted'),
          '-o', str(clean)]
    if use_readchim:
        fm += ['--chimera-report', str(chimera_dir / 'chimera_read_report.tsv')]
    steps.append(Step('filter', 'filter-misassemblies', None, fm))

    return steps


def register_parser(subparsers):
    parser = subparsers.add_parser(
        'misassembly-wf',
        help='Run the whole misassembly-cleaning branch (PacBio HiFi only)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
PacBio HiFi only. Map HiFi reads to their own assembly (minimap2 -ax map-hifi,
keep supplementary), sort and index, then:

  mamisa misassembly-wf \\
    --assembly assembly.fa --bam hifi_mapping.bam -o misasm_wf/ --threads 20 \\
    --anvio-env anvio-9 --checkm2-env checkm2

Add read-level taxonomy (enables the check-read-chimeras step):
    --kraken2-output reads.kraken2 --kraken2-report report.txt
  or
    --kaiju-output reads.kaiju --kaiju-names names.dmp

Steps: anvi -> process -> readchim (if taxonomy) -> classify -> filter.
Resume after a fix with --from <step>; drop steps with --skip a,b; preview with
--dry-run.
        """
    )
    parser.add_argument('-a', '--assembly', type=Path, required=True,
                        help='Assembly FASTA (the reads were mapped to this)')
    parser.add_argument('-b', '--bam', type=Path, required=True,
                        help='Sorted, indexed HiFi BAM mapped to the assembly')
    parser.add_argument('-o', '--output-dir', type=Path, required=True,
                        help='Base output directory for all steps')
    parser.add_argument('--threads', type=int, default=1, help='Threads (default: 1)')

    parser.add_argument('--anvio-env', default='anvio-9',
                        help='Conda env with anvi\'o + MaMISA (default: anvio-9)')
    parser.add_argument('--checkm2-env', default='checkm2',
                        help='Conda env with CheckM2 (default: checkm2)')
    parser.add_argument('--anvi-prefix', default='MisAsm',
                        help='Output prefix for anvi clipping files (default: MisAsm)')

    parser.add_argument('--max-length', type=int, default=300000,
                        help='process-large-contigs --max-length (default: 300000)')
    parser.add_argument('--min-completeness', type=float, default=50.0,
                        help='process-large-contigs --min-completeness (default: 50)')
    parser.add_argument('--max-contamination', type=float, default=10.0,
                        help='process-large-contigs --max-contamination (default: 10)')

    tax = parser.add_mutually_exclusive_group()
    tax.add_argument('--kraken2-output', type=Path,
                     help='Kraken2 per-read output (enables check-read-chimeras)')
    tax.add_argument('--kaiju-output', type=Path,
                     help='Kaiju per-read output (enables check-read-chimeras)')
    parser.add_argument('--kraken2-report', type=Path, help='Kraken2 report (for names)')
    parser.add_argument('--kaiju-names', type=Path, help='NCBI names.dmp (for Kaiju names)')

    parser.add_argument('--from', dest='from_step', metavar='STEP',
                        help=f'Start at this step ({", ".join(ALL_KEYS)})')
    parser.add_argument('--skip', dest='skip_steps', metavar='STEP[,STEP]',
                        help='Comma-separated steps to skip')
    parser.add_argument('--dry-run', action='store_true',
                        help='Print the commands without running them')
    parser.set_defaults(func=run)
    return parser


def run(args):
    print_header("MaMISA - Misassembly workflow")

    if not args.dry_run:
        for sub in ('01_misassemblies', '02_large_contigs', '03_chimera',
                    '04_classified'):
            (args.output_dir / sub).mkdir(parents=True, exist_ok=True)

    skip = [s.strip() for s in args.skip_steps.split(',')] if args.skip_steps else None
    steps = build_plan(args)
    code = execute_plan(steps, args.from_step, skip, args.dry_run)
    if code != 0:
        raise SystemExit(code)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    register_parser(parser.add_subparsers(dest='command'))
    run(parser.parse_args())
