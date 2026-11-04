"""Real workflow execution; fixture producers stand in only for heavy tools."""
import csv
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

from echidna.utils.bin_overlap import sequence_digest

ROOT = Path(__file__).resolve().parents[1]
SNAKEMAKE = shutil.which('snakemake') or str(Path(sys.executable).parent / 'snakemake')


def setup_case(tmp_path):
    work = tmp_path / 'work with spaces'
    selected, candidates = work / 'selected', work / 'candidates'
    selected.mkdir(parents=True)
    candidates.mkdir()
    for root, name, seqs in [(selected, 'high', ['AACCGT']), (selected, 'original', ['AACCGT', 'CCCCAA']),
                             (candidates, 'bad', ['CCCCAA', 'GGGATA']), (candidates, 'good', ['ACACCC'])]:
        (root / (name + '.fna')).write_text(''.join(f'>{name}_{i}\n{seq}\n' for i, seq in enumerate(seqs)))
    quality = work / 'quality.tsv'
    quality.write_text('Name\tCompleteness\tContamination\nhigh\t100\t0\noriginal\t95\t0\nbad\t98\t0\ngood\t90\t0\n')
    clusters = work / 'clusters.tsv'
    clusters.write_text('genome\tcluster\nhigh.fna\t1\noriginal.fna\t2\nbad.fna\t2\ngood.fna\t2\n')
    taxonomy = work / 'taxonomy.tsv'
    with taxonomy.open('w', newline='') as stream:
        writer = csv.writer(stream, delimiter='\t')
        writer.writerow(['sequence_sha256', 'contig_length', 'taxonomy_rank', 'rank_taxon'])
        for seq, taxon in [('AACCGT', 'A'), ('CCCCAA', 'A'), ('GGGATA', 'B'), ('ACACCC', 'A')]:
            writer.writerow([sequence_digest(seq), len(seq), 'order', taxon])
    db = work / 'db'
    db.mkdir()
    (db / 'gunc.dmnd').write_text('fixture')
    bindir = work / 'bin'
    bindir.mkdir()
    cli = bindir / 'echidna'
    cli.write_text('#!' + sys.executable + '\n' + (ROOT / 'tests/fixtures/workflow_cli.py').read_text())
    cli.chmod(0o755)
    env = dict(os.environ, PATH=str(bindir) + os.pathsep + os.environ['PATH'],
               PYTHONPATH=str(ROOT), XDG_CACHE_HOME=str(work / 'cache'))
    return work, selected, candidates, quality, clusters, taxonomy, db, cli, env


def read_rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


@pytest.mark.parametrize('mode', ['exclude-mags', 'remove-contigs'])
@pytest.mark.skipif(not Path(SNAKEMAKE).is_file(), reason='Snakemake not installed')
def test_snakemake_replacement_gets_fresh_qc_and_reaggregated_taxonomy(tmp_path, mode):
    work, selected, candidates, quality, clusters, taxonomy, db, cli, env = setup_case(tmp_path)
    out = work / 'results'
    config = work / 'config.json'
    settings = {'genomes_dir': str(selected), 'genome_ext': 'fna', 'outdir': str(out), 'threads': 1,
                'gtdbtk_data': str(db), 'gunc_db': str(db / 'gunc.dmnd'), 'mimag_rna': mode == 'remove-contigs',
                'overlap_resolution': mode, 'overlap_quality_report': str(quality),
                'overlap_cluster_table': str(clusters), 'overlap_candidates_dirs': [str(candidates)],
                'contig_taxonomy': str(taxonomy), 'use_named_envs': False}
    config.write_text(json.dumps(settings))
    cmd = [SNAKEMAKE, '--cores', '1', '-s', str(ROOT / 'workflow/Snakefile'), '--configfile', str(config)]
    def run(*extra):
        result = subprocess.run(cmd + list(extra), cwd=work, env=env, capture_output=True, text=True, timeout=90)
        assert result.returncode == 0, result.stdout + result.stderr
        return result.stdout + result.stderr
    run()
    assert {row['Name'] for row in read_rows(out / 'checkm2/quality_report.tsv')} == {'high', 'good'}
    assert {row['bin'] for row in read_rows(out / 'contig_taxonomy.tsv.bins.tsv')} == {'high', 'good'}
    assert {row['chimera_verdict'] for row in read_rows(out / 'chimera_report.tsv')} == {'Clean'}
    assert {path.stem for path in (out / 'overlap_resolution/retained').glob('*.fa')} == {'high', 'good'}
    if mode == 'remove-contigs':
        assert {row['genome'] for row in read_rows(out / 'mimag_rna/mimag_rna_summary.tsv')} == {'high', 'good'}
    assert 'Nothing to be done' in run('-n')
    # Candidate changes must rerun selection and all downstream QC, without stale copies.
    fasta = candidates / 'good.fna'
    fasta.write_text(fasta.read_text().replace('good_0', 'good_renamed'))
    assert 'resolve_overlaps' in run('-n')
    run()
    assert {row['Name'] for row in read_rows(out / 'checkm2/quality_report.tsv')} == {'high', 'good'}
    assert 'Nothing to be done' in run('-n')
    # A changed candidate sequence with no taxonomy must remove the old winner
    # and propagate a changed genome list through the checkpoint and tier outputs.
    fasta.write_text('>new_unknown_sequence\nAGAGTC\n')
    run()
    expected = {'high'} if mode == 'exclude-mags' else {'high', 'original'}
    assert {row['Name'] for row in read_rows(out / 'checkm2/quality_report.tsv')} == expected
    assert {row['bin'] for row in read_rows(out / 'chimera_report.tsv')} == expected
    assert not (out / 'overlap_resolution/retained/good.fa').exists()
    assert len(list((out / 'organized/Selected/HQ').glob('*.fa'))) == len(expected)
    assert 'Nothing to be done' in run('-n')


@pytest.mark.parametrize('mode', ['exclude-mags', 'remove-contigs'])
def test_chained_workflow_checks_the_replacement_instead_of_old_representative(tmp_path, mode):
    work, selected, candidates, quality, clusters, taxonomy, db, cli, env = setup_case(tmp_path)
    out = work / 'chained'
    cmd = [str(cli), 'mag-quality-wf', '--bins-dir', str(selected), '--genome-ext', 'fna', '-o', str(out),
           '--resolve-overlaps', mode, '--quality-report', str(quality), '--contig-taxonomy', str(taxonomy),
           '--cluster-table', str(clusters), '--candidates-dir', str(candidates), '--gtdbtk-data', str(db),
           '--gunc-db', str(db / 'gunc.dmnd'), '--checkm2-env', '', '--gtdbtk-env', '', '--gunc-env', '']
    cmd += ['--mimag-env', ''] if mode == 'remove-contigs' else ['--skip-mimag']
    result = subprocess.run(cmd, cwd=work, env=env, capture_output=True, text=True, timeout=90)
    assert result.returncode == 0, result.stdout + result.stderr
    assert {row['Name'] for row in read_rows(out / '01_checkm2/quality_report.tsv')} == {'high', 'good'}
    rows = read_rows(out / '04b_contig_taxonomy/classify_contigs.tsv.bins.tsv')
    assert {row['bin'] for row in rows} == {'high', 'good'}
    assert {row['foreign_frac'] for row in rows} == {'0.0000'}
    assert {row['bin'] for row in read_rows(out / '06_chimera/chimera_report.tsv')} == {'high', 'good'}
    assert len(list((out / '07_final/Selected/HQ').glob('*.fa'))) == 2
    if mode == 'remove-contigs':
        assert {row['genome'] for row in read_rows(out / '05_mimag/mimag_rna_summary.tsv')} == {'high', 'good'}
