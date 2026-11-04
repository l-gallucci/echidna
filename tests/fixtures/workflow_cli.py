"""Fixture producers for expensive tools; downstream Echidna commands run normally.

This checks DAG execution, paths and reruns. It does not validate the external
scientific tools or replace a benchmark on real genomes and databases.
"""
import csv
import os
from pathlib import Path
import sys

command = sys.argv[1]
args = sys.argv[2:]
if command not in ('run-checkm2', 'run-gtdbtk', 'run-gunc', 'run-mimag-rna'):
    os.execv(sys.executable, [sys.executable, '-m', 'echidna'] + sys.argv[1:])


def arg(name):
    if name == '--output' and name not in args:
        name = '-o'
    return args[args.index(name) + 1]


out = Path(arg('--output'))
out.mkdir(parents=True, exist_ok=True)
ext = arg('--extension') if '--extension' in args else arg('--file-suffix').lstrip('.')
selected = '--selected-dir' in args
genomes = Path(arg('--selected-dir' if selected else '--genome-dir'))
names = sorted(p.stem for p in (genomes.rglob('*.' + ext) if selected else genomes.glob('*.' + ext)))


def table(name, fields, rows):
    with (out / name).open('w', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(fields)
        w.writerows(rows)


if command == 'run-checkm2':
    if '--database' in args:
        assert Path(arg('--database')).is_file()
    table('quality_report.tsv', ['Name', 'Completeness', 'Contamination'],
          [[n, 95, 1] for n in names])
elif command == 'run-gtdbtk':
    assert Path(os.environ['GTDBTK_DATA_PATH']).is_dir()
    table('gtdbtk.bac120.summary.tsv', ['user_genome', 'classification', 'msa_percent', 'warnings'],
          [[n, 'd__Bacteria;p__Firmicutes;g__Test_A;s__', 95, ''] for n in names])
elif command == 'run-gunc':
    assert Path(arg('--db-file')).is_file()
    table('GUNC.gtdb_214.maxCSS_level.tsv',
          ['genome', 'clade_separation_score', 'contamination_portion', 'pass.GUNC'],
          [[n, 0.9 if 'bad' in n else 0.01, 0.1, 'bad' not in n] for n in names])
else:
    table('mimag_rna_summary.tsv',
          ['genome', 'n_5S', 'n_16S', 'n_23S', 'n_trna_aa',
           'rrna_complete', 'trna_complete', 'mimag_rna_ok'],
          [[n, 1, 1, 1, 20, True, True, True] for n in names])
