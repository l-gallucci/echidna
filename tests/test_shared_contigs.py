"""Shared full-contig reporting: identity, strand, collisions and CLI contracts."""
import argparse
import csv
import gzip
import json
from pathlib import Path
import subprocess
import sys

import pytest

from echidna.commands import check_shared_contigs as cli
from echidna.utils.bin_overlap import analyze_bins, discover_bins, read_bin, write_reports


def fasta(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path


def table(path):
    with path.open() as f:
        return list(csv.DictReader(f, delimiter='\t'))


def test_exact_sharing_descriptions_and_renamed_reverse_strand(tmp_path):
    bins = tmp_path / 'bins'
    fasta(bins / 'A.fa', '>ctg_1 circular=true\naac cgt\n>unique_A\nATAT\n')
    fasta(bins / 'B.fa', '>renamed\nACGGTT\n>unique_B\nGGGA\n')
    result = analyze_bins(discover_bins([bins]))
    assert len(result['shared']) == 1
    assert not result['conflicts']
    pair = result['pairs'][('A.fa', 'B.fa')]
    assert pair['n_shared_sequences'] == 1
    assert pair['shared_unique_bp'] == 6
    assert pair['shared_fraction_a'] == pytest.approx(0.6)
    assert pair['shared_fraction_b'] == pytest.approx(0.6)
    assert pair['sequence_jaccard'] == pytest.approx(1 / 3)
    out = tmp_path / 'reports'
    summary = write_reports(result, out)
    rows = table(out / 'shared_contigs.tsv')
    assert {r['contig_id'] for r in rows} == {'ctg_1', 'renamed'}
    assert {r['n_bins'] for r in rows} == {'2'}
    assert len({r['sequence_sha256'] for r in rows}) == 1
    assert summary['n_bins_with_shared_contigs'] == 2
    assert json.loads((out / 'summary.json').read_text()) == summary


def test_shared_ids_use_the_entire_header_token(tmp_path):
    a = fasta(tmp_path / 'A.fa', '>NODE_1_length_6_cov_1\nAACCGT\n')
    b = fasta(tmp_path / 'B.fa', '>NODE_2_length_6_cov_1\nAACCAA\n')
    result = analyze_bins({'A': a, 'B': b})
    assert not result['shared'] and not result['conflicts'] and not result['pairs']


def test_same_id_different_sequence_is_a_conflict_not_shared(tmp_path):
    a = fasta(tmp_path / 'A.fa', '>contig_1\nAACCGT\n')
    b = fasta(tmp_path / 'B.fa', '>contig_1\nTTAACC\n')
    result = analyze_bins({'A': a, 'B': b})
    summary = write_reports(result, tmp_path / 'reports')
    assert summary['n_shared_sequences'] == 0
    assert summary['n_conflicting_ids'] == 1
    assert not table(tmp_path / 'reports/shared_contigs.tsv')
    assert not table(tmp_path / 'reports/bin_pairs.tsv')
    assert len(table(tmp_path / 'reports/id_conflicts.tsv')) == 2
    assert all(b['n_conflicting_ids'] == 1 for b in result['bins'].values())


def test_shared_sequence_can_coexist_with_conflicting_id(tmp_path):
    files = {name: fasta(tmp_path / (name + '.fa'), '>x\n' + seq + '\n')
             for name, seq in [('A', 'AACCGT'), ('B', 'ACGGTT'), ('C', 'TTAACC')]}
    result = analyze_bins(files)
    assert len(result['shared']) == 1 and len(result['conflicts']) == 1
    assert list(result['pairs']) == [('A', 'B')]
    assert result['bins']['C']['n_shared_contigs'] == 0


def test_duplicate_sequence_within_bin_does_not_inflate_pairs(tmp_path):
    a = fasta(tmp_path / 'A.fa', '>x\nAACCGT\n>copy\nACGGTT\n')
    b = fasta(tmp_path / 'B.fa', '>y\nAACCGT\n')
    result = analyze_bins({'A': a, 'B': b})
    pair = result['pairs'][('A', 'B')]
    assert pair['n_shared_sequences'] == 1 and pair['shared_unique_bp'] == 6
    assert pair['unique_bp_a'] == pair['unique_bp_b'] == 6
    assert pair['shared_fraction_a'] == pair['shared_fraction_b'] == 1
    assert result['bins']['A']['n_contigs'] == 2
    assert result['bins']['A']['n_unique_sequences'] == 1
    assert result['bins']['A']['shared_bp'] == 12
    assert result['bins']['A']['shared_fraction_bp'] == 1
    alone = analyze_bins({'A': a})
    assert not alone['shared'] and alone['bins']['A']['n_shared_contigs'] == 0


def test_three_bins_pair_counts_and_per_bin_bases(tmp_path):
    files = {name: fasta(tmp_path / (name + '.fa'), '>x\nAACCGT\n') for name in 'ABC'}
    result = analyze_bins(files)
    summary = write_reports(result, tmp_path / 'reports')
    assert len(result['pairs']) == 3
    assert summary['n_shared_sequences'] == 1
    assert summary['n_pairs_with_shared_contigs'] == 3
    assert all(row['n_shared_contigs'] == 1 and row['shared_bp'] == 6 for row in result['bins'].values())


def test_clean_bins_produce_header_only_shared_reports(tmp_path):
    a = fasta(tmp_path / 'A.fa', '>a\nAAAA\n')
    b = fasta(tmp_path / 'B.fa', '>b\nCCCC\n')
    out = tmp_path / 'reports'
    summary = write_reports(analyze_bins({'A': a, 'B': b}), out)
    assert summary['n_shared_sequences'] == summary['n_conflicting_ids'] == 0
    for name in ('shared_contigs.tsv', 'id_conflicts.tsv', 'bin_pairs.tsv'):
        assert not table(out / name)
        assert len((out / name).read_text().splitlines()) == 1
    assert len(table(out / 'bin_summary.tsv')) == 2


def test_recursive_gzip_and_multiple_roots_preserve_bin_paths(tmp_path):
    first = tmp_path / 'first'
    second = tmp_path / 'second'
    fasta(first / 'HQ/bin.fa', '>a\nAACCGT\n')
    fasta(first / 'MQ/bin.fa', '>b\nTTAACC\n')
    second.mkdir()
    with gzip.open(second / 'bin.FNA.GZ', 'wt') as f:
        f.write('>renamed\nACGGTT\n')
    files = discover_bins([first, second], recursive=True)
    assert set(files) == {'set1/HQ/bin.fa', 'set1/MQ/bin.fa', 'set2/bin.FNA.GZ'}
    result = analyze_bins(files)
    assert list(result['pairs']) == [('set1/HQ/bin.fa', 'set2/bin.FNA.GZ')]
    assert not (first / 'HQ/shared_contigs.tsv').exists()


def test_iupac_reverse_complement(tmp_path):
    a = fasta(tmp_path / 'A.fa', '>a\nACGTRYSWKMBDHVN\n')
    b = fasta(tmp_path / 'B.fa', '>b\nNBDHVKMWSRYACGT\n')
    assert len(analyze_bins({'A': a, 'B': b})['shared']) == 1


@pytest.mark.parametrize('text,problem', [
    ('', 'No FASTA records'),
    ('ACGT\n>a\nACGT\n', 'before first FASTA header'),
    ('>\nACGT\n', 'Empty FASTA header'),
    ('>a\n', 'Empty sequence'),
    ('>a\n>b\nACGT\n', 'Empty sequence'),
    ('>a\nACGTZ\n', 'Invalid DNA characters'),
    ('>a\nACGT\n>a duplicate\nACGT\n', 'Duplicate contig ID'),
])
def test_malformed_inputs_are_errors(tmp_path, text, problem):
    path = fasta(tmp_path / 'bad.fa', text)
    with pytest.raises(ValueError, match=problem):
        list(read_bin(path))


def test_missing_empty_and_duplicate_inputs(tmp_path):
    empty = tmp_path / 'empty'
    empty.mkdir()
    with pytest.raises(ValueError, match='directory not found'):
        discover_bins([tmp_path / 'absent'])
    with pytest.raises(ValueError, match='No matching FASTA'):
        discover_bins([empty])
    fasta(empty / 'bin.fa', '>x\nACGT\n')
    with pytest.raises(ValueError, match='included more than once'):
        discover_bins([empty, empty])
    (empty / 'alias.fa').symlink_to(empty / 'bin.fa')
    aliases = analyze_bins(discover_bins([empty]))
    assert list(aliases['pairs']) == [('alias.fa', 'bin.fa')]


def test_existing_report_is_not_overwritten(tmp_path):
    path = fasta(tmp_path / 'bins/A.fa', '>a\nACGT\n')
    out = tmp_path / 'out'
    out.mkdir()
    marker = out / 'summary.json'
    marker.write_text('previous results')
    with pytest.raises(ValueError, match='new or empty'):
        write_reports(analyze_bins({'A': path}), out)
    assert marker.read_text() == 'previous results'


def test_cli_and_invalid_input_without_partial_reports(tmp_path):
    root = Path(__file__).resolve().parents[1]
    bins = tmp_path / 'bins with spaces'
    a = fasta(bins / 'A.fa', '>shared\nAACCGT\n')
    b = fasta(bins / 'B.fa', '>shared\nAACCGT\n')
    before = (a.read_bytes(), b.read_bytes())
    out = tmp_path / 'output with spaces'
    cmd = [sys.executable, '-m', 'echidna', 'check-shared-contigs', '--bins-dir', str(bins), '-o', str(out)]
    proc = subprocess.run(cmd, cwd=root, capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    assert len(table(out / 'shared_contigs.tsv')) == 2
    assert (a.read_bytes(), b.read_bytes()) == before
    assert subprocess.run(cmd, cwd=root, capture_output=True).returncode == 1
    fasta(b, '>bad\nXYZ\n')
    bad_out = tmp_path / 'invalid'
    cmd[-1] = str(bad_out)
    failed = subprocess.run(cmd, cwd=root, capture_output=True, text=True)
    assert failed.returncode == 1 and 'Invalid DNA characters' in failed.stderr
    assert not bad_out.exists()


def test_parser_accepts_multiple_dirs_and_recursive():
    parser = argparse.ArgumentParser()
    cli.register_parser(parser.add_subparsers())
    args = parser.parse_args(['check-shared-contigs', '--bins-dir', 'a', 'b', '--recursive', '-o', 'out'])
    assert args.bins_dir == [Path('a'), Path('b')] and args.recursive
