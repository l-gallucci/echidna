#!/usr/bin/env python3
"""
Echidna - mag-quality-wf command

Chains the MAG quality + taxonomy branch end to end, each sub-command in its env:

  checkm2    run-checkm2          (CheckM2 env)
  organize1  organize-mags        tier the bins (original names, symlinks)
  gtdbtk     run-gtdbtk           (GTDB-Tk env; GTDBTK_DATA_PATH)
  gunc       run-gunc             (GUNC env)
  mmseqs     run-mmseqs-contigs   (MMseqs2 env; only with --mmseqs-db)
  contigtax  classify-contigs     per-contig taxonomy -> foreign-fraction TSV
  mimag      run-mimag-rna        (MIMAG env; optional, full-MIMAG HQ gate)
  chimeras   check-chimeras       verdict from GUNC (+ per-contig taxonomy)
  organize2  organize-mags        final rename + quality gate

The optional overlap step reuses prior contig taxonomy to gate resolved MAGs
before final QC. The chimera verdict uses GUNC alone unless --mmseqs-db or
--contig-taxonomy adds the per-contig signal. A supplied sequence-linked table
is reaggregated, without rerunning the classifier.

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
from .check_shared_contigs import add_taxonomy_options

ALL_KEYS = ['overlap', 'checkm2', 'organize1', 'gtdbtk', 'gunc', 'mmseqs', 'contigtax',
            'mimag', 'chimeras', 'organize2']


def build_plan(args) -> List[Step]:
    """Build the ordered Step list from parsed args. Pure (no side effects)."""
    out = args.output_dir
    resolution_mode = getattr(args, 'resolve_overlaps', None)
    taxonomy_table = getattr(args, 'contig_taxonomy', None)
    taxonomy_rank = getattr(args, 'taxonomy_rank', None) or 'order'
    foreign_threshold = getattr(args, 'min_foreign_frac', 0.05)
    if taxonomy_table and args.mmseqs_db:
        raise ValueError('Use --contig-taxonomy to reuse classifications OR --mmseqs-db to generate them')
    if resolution_mode and (not taxonomy_table or not getattr(args, 'quality_report', None)):
        raise ValueError('Workflow overlap resolution requires --quality-report and --contig-taxonomy '
                         '(sequence-linked per-contig TSV covering representatives and candidate MAGs)')
    if not resolution_mode and any(getattr(args, flag, None) for flag in
                                   ('quality_report', 'rank_by', 'drep_dir', 'cluster_table', 'candidates_dir')):
        raise ValueError('Quality ranking and cluster replacement options require --resolve-overlaps')
    bins = out / '00_overlap' / 'retained' if resolution_mode else args.bins_dir
    ext = 'fa' if resolution_mode else args.genome_ext.lstrip('.')
    suffix = '.' + ext

    ckm2_dir = out / '01_checkm2'
    filtered_dir = out / '02_filtered'
    selected = filtered_dir / 'Selected'
    tax_dir = out / '03_taxonomy'
    gunc_dir = out / '04_gunc'
    ctax_dir = out / '04b_contig_taxonomy'
    mmseqs_tsv = ctax_dir / 'contig_taxonomy.tsv'
    classify_out = ctax_dir / 'classify_contigs.tsv'
    contig_bins_tsv = Path(str(classify_out) + '.bins.tsv')
    mimag_dir = out / '05_mimag'
    chimera_report = out / '06_chimera' / 'chimera_report.tsv'
    final_dir = out / '07_final'
    quality_report = ckm2_dir / 'quality_report.tsv'

    use_mimag = not args.skip_mimag

    steps: List[Step] = []
    if resolution_mode:
        overlap = ['echidna', 'check-shared-contigs', '--bins-dir', str(args.bins_dir),
                   '--recursive', '--extensions', args.genome_ext.lstrip('.'),
                   '--' + ('exclude-mags' if resolution_mode == 'exclude-mags' else 'remove-shared-contigs'),
                   '--quality-report', str(args.quality_report), '--contig-taxonomy', str(taxonomy_table),
                   '--require-retained',
                   '--taxonomy-rank', taxonomy_rank,
                   '--min-taxonomy-resolved-frac', str(getattr(args, 'min_taxonomy_resolved_frac', 0.5)),
                   '--min-foreign-frac', str(foreign_threshold), '-o', str(out / '00_overlap')]
        for attr in ('rank_by', 'drep_dir', 'cluster_table'):
            value = getattr(args, attr, None)
            if value:
                overlap += ['--' + attr.replace('_', '-'), str(value)]
        if getattr(args, 'candidates_dir', None):
            overlap += ['--candidates-dir'] + [str(path) for path in args.candidates_dir]
        steps.append(Step('overlap', 'check-shared-contigs (overlap + taxonomy gate)', args.echidna_env, overlap))

    ck = ['echidna', 'run-checkm2', '--genome-dir', str(bins),
          '-o', str(ckm2_dir), '--threads', str(args.threads), '--extension', ext]
    if args.checkm2_db:
        ck += ['--database', str(args.checkm2_db)]
    steps.append(Step('checkm2', 'run-checkm2', args.checkm2_env, ck))

    steps.append(Step(
        'organize1', 'organize-mags (tier)', args.echidna_env,
        ['echidna', 'organize-mags', '--checkm2-root', str(ckm2_dir),
         '--genomes-dir', str(bins), '-o', str(filtered_dir),
         '--tiers', args.tiers, '--symlink', '--extensions', ext],
    ))

    gt_env_vars = {'GTDBTK_DATA_PATH': str(args.gtdbtk_data)} if args.gtdbtk_data else None
    steps.append(Step(
        'gtdbtk', 'run-gtdbtk', args.gtdbtk_env,
        ['echidna', 'run-gtdbtk', '--selected-dir', str(selected),
         '-o', str(tax_dir), '--cpus', str(args.threads),
         '--tiers', args.tiers, '--extension', ext],
        gt_env_vars,
    ))

    gu = ['echidna', 'run-gunc', '--selected-dir', str(selected),
          '-o', str(gunc_dir), '--file-suffix', suffix,
          '--threads', str(args.threads)]
    if args.gunc_db:
        gu += ['--db-file', str(args.gunc_db)]
    steps.append(Step('gunc', 'run-gunc', args.gunc_env, gu))

    # Per-contig taxonomy can be generated or reused, always reaggregated for the
    # selected genomes rather than copying the old representative's bin summary.
    use_contig_tax = args.mmseqs_db is not None or taxonomy_table is not None
    if args.mmseqs_db is not None:
        steps.append(Step(
            'mmseqs', 'run-mmseqs-contigs', args.mmseqs_env,
            ['echidna', 'run-mmseqs-contigs', '--bins-dir', str(selected),
             '--db', str(args.mmseqs_db), '-o', str(mmseqs_tsv),
             '--threads', str(args.threads), '--extensions', ext],
        ))
    if use_contig_tax:
        source = ['--taxonomy-table', str(taxonomy_table)] if taxonomy_table else ['--mmseqs-output', str(mmseqs_tsv)]
        steps.append(Step(
            'contigtax', 'classify-contigs', args.echidna_env,
            ['echidna', 'classify-contigs', '--bins-dir', str(selected),
             '-o', str(classify_out), '--rank', taxonomy_rank,
             '--min-foreign-frac', str(foreign_threshold), '--extensions', ext] + source,
        ))

    if use_mimag:
        mi = ['echidna', 'run-mimag-rna', '--selected-dir', str(selected),
              '-o', str(mimag_dir), '--threads', str(args.threads),
              '--tiers', args.tiers, '--extension', ext]
        if args.split_by_domain:
            mi += ['--split-by-domain', '--gtdbtk-dir', str(tax_dir)]
        steps.append(Step('mimag', 'run-mimag-rna', args.mimag_env, mi))

    chim = ['echidna', 'check-chimeras', '--bins-dir', str(selected),
            '--gunc-dir', str(gunc_dir), '--checkm2-report', str(quality_report),
            '--gtdbtk-dir', str(tax_dir), '--extensions', ext, '-o', str(chimera_report)]
    if use_contig_tax:
        chim += ['--contig-taxonomy', str(contig_bins_tsv), '--min-foreign-frac', str(foreign_threshold)]
    steps.append(Step('chimeras', 'check-chimeras', args.echidna_env, chim))

    o2 = ['echidna', 'organize-mags', '--checkm2-root', str(ckm2_dir),
          '--genomes-dir', str(bins), '-o', str(final_dir),
          '--tiers', args.tiers, '--copy', '--rename',
          '--gtdbtk-dir', str(tax_dir), '--extensions', ext]
    if use_mimag:
        o2 += ['--mimag-rna-dir', str(mimag_dir)]
    steps.append(Step('organize2', 'organize-mags (final rename + gate)',
                      args.echidna_env, o2))

    return steps


def register_parser(subparsers):
    parser = subparsers.add_parser(
        'mag-quality-wf',
        help='Run the whole MAG quality + taxonomy branch from a bin set',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Binning is done upstream; this takes the bins.

  echidna mag-quality-wf \\
    --bins-dir bins/ -o mag_wf/ --threads 40 --tiers HQ,MQ \\
    --gtdbtk-data /data/gtdbtk_r220 \\
    --gunc-db /data/gunc/gunc_db_gtdb214.dmnd \\
    --mmseqs-db /data/mmseqs_gtdb/gtdb \\
    --checkm2-env checkm2 --gtdbtk-env gtdbtk --gunc-env gunc --mimag-env mimag \\
    --split-by-domain

Steps: [overlap] -> checkm2 -> organize1 -> gtdbtk -> gunc -> [mmseqs -> contigtax]
-> mimag -> chimeras -> organize2. --mmseqs-db generates contig classifications;
--contig-taxonomy reuses a sequence-linked classify-contigs TSV instead.
--resolve-overlaps requires that table plus --quality-report for prior ranking.
Optional --drep-dir/--cluster-table with --candidates-dir enables replacements.
All final QC follows the resolved genomes; the original bin reports are not reused.
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
    parser.add_argument('--resolve-overlaps', choices=['exclude-mags', 'remove-contigs'],
                        help='Before QC, resolve overlaps and gate all retained MAGs by reused contig taxonomy')
    parser.add_argument('--quality-report', type=Path,
                        help='Existing quality table for overlap ranking, including candidate MAGs; '
                             'CheckM2 is run again on the resolved collection')
    parser.add_argument('--rank-by', choices=['completeness', 'quality-score'])
    cluster = parser.add_mutually_exclusive_group()
    cluster.add_argument('--drep-dir', type=Path)
    cluster.add_argument('--cluster-table', type=Path)
    parser.add_argument('--candidates-dir', type=Path, nargs='+')
    add_taxonomy_options(parser)

    parser.add_argument('--echidna-env', default=None,
                        help='Env for the pure-Echidna steps (default: current env)')
    parser.add_argument('--checkm2-env', default='checkm2',
                        help='Conda env with CheckM2 (default: checkm2)')
    parser.add_argument('--gtdbtk-env', default='gtdbtk',
                        help='Conda env with GTDB-Tk (default: gtdbtk)')
    parser.add_argument('--gunc-env', default='gunc',
                        help='Conda env with GUNC (default: gunc)')
    parser.add_argument('--mimag-env', default='mimag',
                        help='Conda env with barrnap + tRNAscan-SE (default: mimag)')
    parser.add_argument('--mmseqs-env', default='mmseqs2',
                        help='Conda env with MMseqs2 (default: mmseqs2); used only '
                             'when --mmseqs-db is given')

    parser.add_argument('--gtdbtk-data', type=Path,
                        help='GTDB-Tk reference data dir (sets GTDBTK_DATA_PATH)')
    parser.add_argument('--gunc-db', type=Path,
                        help='GUNC diamond database (.dmnd); else uses $GUNC_DB')
    parser.add_argument('--checkm2-db', type=Path,
                        help='CheckM2 database (.dmnd); else CheckM2 default')
    parser.add_argument('--mmseqs-db', type=Path,
                        help='MMseqs2 taxonomy DB prefix (e.g. /data/mmseqs_gtdb/gtdb). '
                             'When given, adds the mmseqs + contigtax steps so the '
                             'chimera verdict gets the per-contig taxonomy signal; '
                             'without it the verdict uses GUNC alone.')

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
    print_header("Echidna - MAG quality workflow")

    if not args.dry_run:
        subdirs = ['01_checkm2', '06_chimera']
        if args.mmseqs_db is not None or args.contig_taxonomy is not None:
            subdirs.append('04b_contig_taxonomy')
        for sub in subdirs:
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
