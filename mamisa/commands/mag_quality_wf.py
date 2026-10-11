#!/usr/bin/env python3
"""
MaMISA - mag-quality-wf command

Chains the MAG quality + taxonomy branch end to end, each sub-command in its env:

  checkm2    run-checkm2          (CheckM2 env)
  organize1  organize-mags        tier the bins (original names, symlinks)
  gtdbtk     run-gtdbtk           (GTDB-Tk env; GTDBTK_DATA_PATH)
  gunc       run-gunc             (GUNC env; --contig-taxonomy-output)
  mimag      run-mimag-rna        (MIMAG env; optional, full-MIMAG HQ gate)
  chimeras   check-chimeras       GC + GUNC + taxonomy chimera report
  organize2  organize-mags        final rename + quality gate

The genome names are frozen after organize1 so GUNC, check-chimeras and mimag all
see the same names (no rename until organize2). Start binning yourself; this
workflow takes the bins.

Resume/partial: --from <step>, --skip <step[,step]>, --dry-run.
"""

import argparse
from pathlib import Path
from typing import List

from ..utils.workflow_run import Step, execute_plan
from ..utils.logging import print_header

ALL_KEYS = ['checkm2', 'organize1', 'gtdbtk', 'gunc', 'mimag', 'chimeras', 'organize2']


def build_plan(args) -> List[Step]:
    """Build the ordered Step list from parsed args. Pure (no side effects)."""
    out = args.output_dir
    ext = args.genome_ext.lstrip('.')
    suffix = '.' + ext

    ckm2_dir = out / '01_checkm2'
    filtered_dir = out / '02_filtered'
    selected = filtered_dir / 'Selected'
    tax_dir = out / '03_taxonomy'
    gunc_dir = out / '04_gunc'
    mimag_dir = out / '05_mimag'
    chimera_report = out / '06_chimera' / 'chimera_report.tsv'
    final_dir = out / '07_final'
    quality_report = ckm2_dir / 'quality_report.tsv'

    use_mimag = not args.skip_mimag

    steps: List[Step] = []

    ck = ['mamisa', 'run-checkm2', '--genome-dir', str(args.bins_dir),
          '-o', str(ckm2_dir), '--threads', str(args.threads), '--extension', ext]
    if args.checkm2_db:
        ck += ['--database', str(args.checkm2_db)]
    steps.append(Step('checkm2', 'run-checkm2', args.checkm2_env, ck))

    steps.append(Step(
        'organize1', 'organize-mags (tier)', args.mamisa_env,
        ['mamisa', 'organize-mags', '--checkm2-root', str(ckm2_dir),
         '--genomes-dir', str(args.bins_dir), '-o', str(filtered_dir),
         '--tiers', args.tiers, '--symlink', '--extensions', ext],
    ))

    gt_env_vars = {'GTDBTK_DATA_PATH': str(args.gtdbtk_data)} if args.gtdbtk_data else None
    steps.append(Step(
        'gtdbtk', 'run-gtdbtk', args.gtdbtk_env,
        ['mamisa', 'run-gtdbtk', '--selected-dir', str(selected),
         '-o', str(tax_dir), '--cpus', str(args.threads),
         '--tiers', args.tiers, '--extension', ext],
        gt_env_vars,
    ))

    gu = ['mamisa', 'run-gunc', '--selected-dir', str(selected),
          '-o', str(gunc_dir), '--file-suffix', suffix,
          '--threads', str(args.threads), '--contig-taxonomy-output']
    if args.gunc_db:
        gu += ['--db-file', str(args.gunc_db)]
    steps.append(Step('gunc', 'run-gunc', args.gunc_env, gu))

    if use_mimag:
        mi = ['mamisa', 'run-mimag-rna', '--selected-dir', str(selected),
              '-o', str(mimag_dir), '--threads', str(args.threads),
              '--tiers', args.tiers]
        if args.split_by_domain:
            mi += ['--split-by-domain', '--gtdbtk-dir', str(tax_dir)]
        steps.append(Step('mimag', 'run-mimag-rna', args.mimag_env, mi))

    steps.append(Step(
        'chimeras', 'check-chimeras', args.mamisa_env,
        ['mamisa', 'check-chimeras', '--bins-dir', str(selected),
         '--gunc-dir', str(gunc_dir), '--checkm2-report', str(quality_report),
         '--gtdbtk-dir', str(tax_dir), '-o', str(chimera_report)],
    ))

    o2 = ['mamisa', 'organize-mags', '--checkm2-root', str(ckm2_dir),
          '--genomes-dir', str(args.bins_dir), '-o', str(final_dir),
          '--tiers', args.tiers, '--copy', '--rename',
          '--gtdbtk-dir', str(tax_dir), '--extensions', ext]
    if use_mimag:
        o2 += ['--mimag-rna-dir', str(mimag_dir)]
    steps.append(Step('organize2', 'organize-mags (final rename + gate)',
                      args.mamisa_env, o2))

    return steps


def register_parser(subparsers):
    parser = subparsers.add_parser(
        'mag-quality-wf',
        help='Run the whole MAG quality + taxonomy branch from a bin set',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Binning is done upstream; this takes the bins.

  mamisa mag-quality-wf \\
    --bins-dir bins/ -o mag_wf/ --threads 40 --tiers HQ,MQ \\
    --gtdbtk-data /data/gtdbtk_r220 \\
    --gunc-db /data/gunc/gunc_db_gtdb214.dmnd \\
    --checkm2-env checkm2 --gtdbtk-env gtdbtk --gunc-env gunc --mimag-env mimag \\
    --split-by-domain

Steps: checkm2 -> organize1 -> gtdbtk -> gunc -> mimag -> chimeras -> organize2.
Names are frozen after organize1 so gunc/chimeras/mimag match. Skip the rRNA/tRNA
gate with --skip-mimag. Resume with --from <step>, drop with --skip a,b, preview
with --dry-run.
        """
    )
    parser.add_argument('--bins-dir', type=Path, required=True,
                        help='Directory of bin FASTA files (from binning)')
    parser.add_argument('-o', '--output-dir', type=Path, required=True,
                        help='Base output directory for all steps')
    parser.add_argument('--threads', type=int, default=1, help='Threads (default: 1)')
    parser.add_argument('--tiers', default='HQ,MQ',
                        help='Quality tiers to carry downstream (default: HQ,MQ)')
    parser.add_argument('--genome-ext', default='fa',
                        help='Bin file extension without dot (default: fa)')

    parser.add_argument('--mamisa-env', default=None,
                        help='Env for the pure-MaMISA steps (default: current env)')
    parser.add_argument('--checkm2-env', default='checkm2',
                        help='Conda env with CheckM2 (default: checkm2)')
    parser.add_argument('--gtdbtk-env', default='gtdbtk',
                        help='Conda env with GTDB-Tk (default: gtdbtk)')
    parser.add_argument('--gunc-env', default='gunc',
                        help='Conda env with GUNC (default: gunc)')
    parser.add_argument('--mimag-env', default='mimag',
                        help='Conda env with barrnap + tRNAscan-SE (default: mimag)')

    parser.add_argument('--gtdbtk-data', type=Path,
                        help='GTDB-Tk reference data dir (sets GTDBTK_DATA_PATH)')
    parser.add_argument('--gunc-db', type=Path,
                        help='GUNC diamond database (.dmnd); else uses $GUNC_DB')
    parser.add_argument('--checkm2-db', type=Path,
                        help='CheckM2 database (.dmnd); else CheckM2 default')

    parser.add_argument('--split-by-domain', action='store_true',
                        help='run-mimag-rna picks bac/arc per genome from GTDB-Tk')
    parser.add_argument('--skip-mimag', action='store_true',
                        help='Skip the rRNA/tRNA MIMAG HQ gate entirely')

    parser.add_argument('--from', dest='from_step', metavar='STEP',
                        help=f'Start at this step ({", ".join(ALL_KEYS)})')
    parser.add_argument('--skip', dest='skip_steps', metavar='STEP[,STEP]',
                        help='Comma-separated steps to skip')
    parser.add_argument('--dry-run', action='store_true',
                        help='Print the commands without running them')
    parser.set_defaults(func=run)
    return parser


def run(args):
    print_header("MaMISA - MAG quality workflow")

    if not args.dry_run:
        for sub in ('01_checkm2', '06_chimera'):
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
