"""Real Snakemake execution using fixture producers for the heavy tools."""
import csv
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

SNAKEMAKE = shutil.which('snakemake') or str(Path(sys.executable).parent / 'snakemake')
ROOT = Path(__file__).resolve().parents[1]
pytestmark = pytest.mark.skipif(not Path(SNAKEMAKE).is_file(), reason='Snakemake not installed')


@pytest.mark.parametrize('rna', [True, False])
def test_dag_execution_and_rerun_from_external_directory(tmp_path, rna):
    work = tmp_path / 'work with spaces'
    bins = work / 'bins with spaces'
    bins.mkdir(parents=True)
    for name in ('S1.good', 'S2.bad'):
        (bins / (name + '.fna')).write_text('>ctg circular=true\n' + 'ACGT' * 1000 + '\n')
    db = work / 'db with spaces'
    db.mkdir()
    (db / 'gunc.dmnd').write_text('fixture')
    (db / 'checkm2.dmnd').write_text('fixture')
    out = work / 'results with spaces'
    config = work / 'config.json'
    config.write_text(json.dumps({
        'genomes_dir': str(bins), 'genome_ext': 'fna', 'outdir': str(out), 'threads': 1,
        'gtdbtk_data': str(db), 'gunc_db': str(db / 'gunc.dmnd'),
        'checkm2_db': str(db / 'checkm2.dmnd'), 'mimag_rna': rna,
        'sample_regex': r"^(S\d+)[.]", 'use_named_envs': False,
    }))
    bindir = work / 'bin'
    bindir.mkdir()
    cli = bindir / 'echidna'
    # Use the same Python as pytest, with the actual CLI for downstream steps.
    cli.write_text('#!' + sys.executable + '\n' + (ROOT / 'tests/fixtures/workflow_cli.py').read_text())
    cli.chmod(0o755)
    env = dict(os.environ, PATH=str(bindir) + os.pathsep + os.environ['PATH'],
               PYTHONPATH=str(ROOT), XDG_CACHE_HOME=str(work / 'cache'))
    cmd = [SNAKEMAKE, '--cores', '1', '-s', str(ROOT / 'workflow/Snakefile'),
           '--configfile', str(config)]

    def run(*extra):
        proc = subprocess.run(cmd + list(extra), cwd=work, env=env,
                              capture_output=True, text=True, timeout=90)
        assert proc.returncode == 0, proc.stdout + proc.stderr
        return proc.stdout + proc.stderr

    run()
    with (out / 'chimera_report.tsv').open() as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    assert {row['chimera_verdict'] for row in rows} == {'Clean', 'Check'}
    assert len(list((out / 'organized/Selected/HQ').glob('*.fna'))) == 2
    assert (out / 'mimag_rna/mimag_rna_summary.tsv').exists() is rna
    assert 'Nothing to be done' in run('-n')
    # A changed quality assessment must move files between tiers on a forced rerun.
    report = out / 'checkm2/quality_report.tsv'
    report.write_text(report.read_text().replace('\t95\t', '\t80\t'))
    run('--forcerun', 'organize_mags')
    assert not list((out / 'organized/Selected/HQ').glob('*.fna'))
    assert len(list((out / 'organized/Selected/MQ').glob('*.fna'))) == 2
    # Editing a FASTA (parent-directory mtime unchanged) reruns upstream QC too.
    fasta = bins / 'S1.good.fna'
    fasta.write_text(fasta.read_text() + 'ACGT\n')
    plan = run('-n')
    assert 'rule run_checkm2:' in plan and 'rule run_gtdbtk:' in plan
    run()
    assert len(list((out / 'organized/Selected/HQ').glob('*.fna'))) == 2
    assert not list((out / 'organized/Selected/MQ').glob('*.fna'))
    assert 'Nothing to be done' in run('-n')
