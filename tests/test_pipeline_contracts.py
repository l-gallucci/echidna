"""Regression tests across real producer/consumer schemas, including negative controls.

External programs are mocked; these tests do not claim end-to-end biological
accuracy. Native anvi'o 9 fixtures follow its documented column names.
"""
import argparse
import csv
import io
import json
import subprocess
import shutil
from collections import Counter
from pathlib import Path

import pytest

from echidna.commands import (process_large_contigs as plc, filter_checkm2 as organize,
    filter_misassemblies as filtering, check_read_chimeras as reads,
    classify_clipping as clipping, misassembly_wf, mag_quality_wf, run_mimag_rna,
    check_chimeras, classify_contigs, remove_hq_contigs, final_report)
from echidna.commands.setup_workflow import set_yaml_scalar
from echidna.utils import bam_stats, zero_coverage, misassembly, contig_taxonomy
from echidna.utils.final_report import build_report
from echidna.utils.fasta import read_fasta_streaming
from echidna.utils.read_taxonomy import compute_taxonomy_stats, assess_read_chimera_risk
from echidna.utils.ssu_taxonomy import taxa_agree


def args_for(module, argv):
    parser = argparse.ArgumentParser()
    module.register_parser(parser.add_subparsers(dest='command'))
    return parser.parse_args(argv)


@pytest.fixture
def native_clips(tmp_path):
    path = tmp_path / 'MisAsm-clipping.txt'
    path.write_text('contig\tlength\tpos\trelative_pos\tcov\tclipping\tclipping_ratio\n'
                    'contig_7\t400\t200\t0.5\t120\t120\t1.0\n')
    return tmp_path


def test_native_clipping_counts_and_selective_filter(native_clips):
    assert misassembly.load_clipping_data(native_clips, 50) == {'contig_7': [(200, 120)]}
    assert misassembly.load_clipping_data(native_clips, 121) == {}


def test_legacy_clipping_and_duplicate_positions(tmp_path):
    (tmp_path / 'legacy-clipping.txt').write_text(
        'contig_name\tsample_name\tposition\tclipping_coverage\n'
        'c\ts1\t10\t80\nc\ts2\t10\t100\n')
    assert misassembly.load_clipping_positions(tmp_path) == {'c': [10]}
    assert misassembly.load_clipping_data(tmp_path) == {'c': [(10, 100)]}


def test_no_clipping_file_is_not_a_clean_result(tmp_path):
    with pytest.raises(ValueError, match='No .*clipping'):
        misassembly.load_clipping_data(tmp_path)


@pytest.mark.parametrize('body', ['garbage\n', 'contig\tpos\tclipping\nc\t10\tbad\n'])
def test_invalid_clipping_fails(tmp_path, body):
    (tmp_path / 'x-clipping.txt').write_text(body)
    with pytest.raises(ValueError):
        misassembly.load_clipping_data(tmp_path)


def test_native_zero_cov_coordinates(tmp_path):
    (tmp_path / 'x-zero_cov.txt').write_text('contig\tlength\trange\trange_size\n'
                                          'c\t1200\t0-500\t500\nc\t1200\t700-1200\t500\n')
    assert zero_coverage.parse_zero_cov_files(tmp_path) == [('c', 0, 500), ('c', 700, 1200)]


def test_legacy_zero_cov_coordinates(tmp_path):
    (tmp_path / 'x-zero_cov.txt').write_text('contig_name\tsample_name\tstart\tend\nc\ts\t1\t500\n')
    assert zero_coverage.parse_zero_cov_files(tmp_path) == [('c', 0, 500)]


def test_invalid_zero_cov_fails(tmp_path):
    (tmp_path / 'x-zero_cov.txt').write_text('contig\trange\nc\tinvalid\n')
    with pytest.raises(ValueError):
        zero_coverage.parse_zero_cov_files(tmp_path)


def test_large_contig_keeps_misassembly_and_unassessed_contigs(tmp_path, native_clips):
    assembly = tmp_path / 'a.fa'
    assembly.write_text('>contig_7\n' + 'A' * 400 + '\n>c2\n' + 'C' * 400 + '\n')
    qc = tmp_path / 'qc'
    qc.mkdir()
    (qc / 'quality_report.tsv').write_text('Name\tCompleteness\tContamination\n'
        f'{plc.safe_filename("contig_7")}\t95\t1\n')
    out = tmp_path / 'out'
    plc.run(args_for(plc, ['process-large-contigs', '-a', str(assembly), '-m', str(native_clips),
        '-o', str(out), '--max-length', '100', '--skip-checkm2', '--checkm2-results', str(qc)]))
    assert set(n for n, _ in read_fasta_streaming(out / 'assembly_for_filtering.fa')) == {'contig_7', 'c2'}
    assert not list((out / 'HQ_extracted').glob('*.fa'))
    assert 'contig_7' in (out / 'contig_name_map.tsv').read_text()


def test_no_large_contigs_produces_contractual_output(tmp_path):
    assembly = tmp_path / 'a.fa'
    assembly.write_text('>c\nACGT\n')
    out = tmp_path / 'out'
    plc.run(args_for(plc, ['process-large-contigs', '-a', str(assembly), '-m', str(tmp_path), '-o', str(out)]))
    assert (out / 'assembly_for_filtering.fa').read_text() == assembly.read_text()
    assert (out / 'HQ_extracted').is_dir()


def test_read_window_producer_consumer_roundtrip(tmp_path):
    path = tmp_path / 'windows.tsv'
    reads.write_tsv([{'contig': 'c', 'window_start': 1000, 'window_end': 2000,
                      'is_taxon_shift': True}], reads.WINDOW_FIELDS, path)
    shifts = clipping.load_taxonomy_shifts(path, flank=0)
    assert clipping.has_taxonomy_shift('c', 1500, shifts)
    assert not clipping.has_taxonomy_shift('c', 3000, shifts)


def test_workflow_propagates_assembly_and_extension():
    ma = args_for(misassembly_wf, ['misassembly-wf', '-a', 'a.fa', '-b', 'a.bam', '-o', 'out', '--kraken2-output', 'r.tsv'])
    rc = next(s for s in misassembly_wf.build_plan(ma) if s.key == 'readchim')
    assert rc.cmd[rc.cmd.index('--assembly') + 1] == 'a.fa'
    qa = args_for(mag_quality_wf, ['mag-quality-wf', '--bins-dir', 'bins', '-o', 'out', '--genome-ext', 'fna'])
    plan = mag_quality_wf.build_plan(qa)
    mi = next(s for s in plan if s.key == 'mimag')
    assert args_for(run_mimag_rna, mi.cmd[1:]).extension == 'fna'
    chim = next(s for s in plan if s.key == 'chimeras')
    assert chim.cmd[chim.cmd.index('--extensions') + 1] == 'fna'


def sam_process(monkeypatch, text, exit_code=0):
    class FakeProc:
        returncode = exit_code
        def __init__(self, *a, **kw):
            self.stdout = io.StringIO(text)
        def wait(self):
            return self.returncode
    monkeypatch.setattr(bam_stats.subprocess, 'Popen', FakeProc)


def test_long_reads_overlap_window_and_retain_trailing_sa(monkeypatch):
    # 0-based start 1000, 9kb matched, soft clipping at 10000.
    sam = 'r\t0\tc\t1001\t60\t9000M100S\t*\t0\t0\t' + 'ACGT' * 2275 + '\t*\tSA:Z:other,1,+,100M,60,0;\n'
    sam_process(monkeypatch, sam)
    result = bam_stats.collect_position_stats(Path('mock.bam'), {'c': [(10000, 1)]},
        500, 20, {'c': 1}, {'c': 20000}, 0, 0)[('c', 10000)]
    assert result['n_clipped_reads'] == 1
    assert result['sa_partner_contig'] == 'other'
    assert result['local_depth'] == pytest.approx(500 / 1001)
    assert result['discordant_fraction'] == 'N/A'


def test_uniform_aligned_depth_is_independent_of_read_length(monkeypatch):
    sam = ''.join(f'r{i}\t0\tc\t1\t60\t20000M\t*\t0\t0\t*\t*\n' for i in range(10))
    sam_process(monkeypatch, sam)
    result = bam_stats.collect_position_stats(Path('mock.bam'), {'c': [(10000, 1)]},
        500, 20, {'c': 10}, {'c': 20000}, 0, 0)[('c', 10000)]
    assert result['local_depth'] == 10
    assert result['depth_ratio'] == 1
    assert clipping.classify_position(result)[0] == 'low_confidence'


def test_samtools_failure_does_not_produce_statistics(monkeypatch):
    sam_process(monkeypatch, '', exit_code=2)
    with pytest.raises(RuntimeError, match='samtools view failed'):
        bam_stats.collect_position_stats(Path('mock.bam'), {'c': [(10000, 1)]}, 500, 20, {}, {}, 0, 0)


def test_unpaired_is_not_discordant():
    acc = bam_stats._PosAccumulator(10000, 'c')
    for flag in [0, 16] * 5:
        acc.add_read(flag, 0, 'ACGT' * 25, '100M', 10000, 0, 0)
    result = acc.to_dict('c', 20000, 10)
    assert result['discordant_fraction'] == 'N/A'
    assert clipping.classify_position(result)[0] == 'low_confidence'


def test_blast_failure_is_not_zero_hits(monkeypatch):
    def fail(*a, **kw):
        raise subprocess.CalledProcessError(2, ['blastn'], stderr='synthetic DB error')
    monkeypatch.setattr(zero_coverage.subprocess, 'run', fail)
    with pytest.raises(RuntimeError, match='blastn failed'):
        zero_coverage.blast_regions_vs_reads(Path('q'), Path('db'))


def test_successful_blast_without_hits_retains_method(monkeypatch):
    monkeypatch.setattr(zero_coverage.subprocess, 'run', lambda *a, **k: subprocess.CompletedProcess(a, 0, stdout='', stderr=''))
    assert zero_coverage.blast_regions_vs_reads(Path('q'), Path('db')) == {}
    assert zero_coverage.assign_verdict(0, float('nan'), False) == 'unsupported'
    assert zero_coverage.assign_verdict(1, float('nan'), False) == 'supported'


def test_taxonomy_without_resolved_contigs_is_unknown(tmp_path):
    bins = tmp_path / 'bins'
    bins.mkdir()
    (bins / 'x.fa').write_text('>c\nACGT\n')
    summary = contig_taxonomy.summarize_bin_lineage({'x': ['c']}, {}, 'order', {'c': 4})
    _, rows = classify_contigs.build_rows(summary, {})
    source = tmp_path / 'ct.tsv'
    classify_contigs.write_tsv(rows, classify_contigs.BIN_FIELDS, source)
    out = tmp_path / 'result.tsv'
    check_chimeras.run_chimera_check(bins, out, ['fa'], None, None, None, source, .05, 5000, .1, 'phylum', [], False)
    with out.open() as f:
        assert next(csv.DictReader(f, delimiter='\t'))['chimera_verdict'] == 'Unknown'


def test_measured_clean_signal_stays_clean(tmp_path):
    table = tmp_path / 'ct.tsv'
    table.write_text('bin\tforeign_frac\tn_contigs_resolved\nx\t0\t3\n')
    assert check_chimeras._load_contig_foreign_frac(table) == {'x': 0}
    from echidna.utils.chimera import assess_chimera_verdict
    assert assess_chimera_verdict(False, False)[0] == 'Clean'
    assert assess_chimera_verdict(True, True)[0] == 'Chimeric'


@pytest.mark.parametrize('value', ['bad', 'nan', '-1', '1.1'])
def test_invalid_taxonomy_fraction_fails(tmp_path, value):
    p = tmp_path / 'ct.tsv'
    p.write_text(f'bin\tforeign_frac\nx\t{value}\n')
    with pytest.raises(ValueError):
        check_chimeras._load_contig_foreign_frac(p)


def test_unmatched_reads_are_not_classified():
    stats = compute_taxonomy_stats(Counter({-1: 100}), {}, True)
    assert stats['n_classified'] == stats['n_taxa'] == 0
    assert stats['classified_fraction'] == 0
    assert assess_read_chimera_risk(stats['dominant_fraction'], stats['n_taxa'], 100, 0, 0)[0] == 'Insufficient'


def test_hq_contigs_read_from_bin_contents(tmp_path):
    hq = tmp_path / 'hq'
    hq.mkdir()
    (hq / 'bin.42.fa').write_text('>contig_7\nAAAA\n>contig_9\nCCCC\n')
    assembly = tmp_path / 'a.fa'
    assembly.write_text('>contig_7\nAAAA\n>contig_9\nCCCC\n>contig_42\nGGGG\n')
    out = tmp_path / 'out.fa'
    ids = remove_hq_contigs.load_hq_ids_from_dir(hq)
    assert ids == filtering.get_hq_contig_ids(hq) == {'contig_7', 'contig_9'}
    remove_hq_contigs.process_assembly(assembly, ids, out)
    assert list(read_fasta_streaming(out)) == [('contig_42', 'GGGG')]


def test_reused_tier_directory_cannot_leave_obsolete_hq(tmp_path):
    bins = tmp_path / 'bins'
    bins.mkdir()
    (bins / 'x.fa').write_text('>c\nACGT\n')
    out = tmp_path / 'out'
    organize.link_or_copy_genomes([{'name': 'x', 'tier': 'HQ'}], bins, out, ['HQ', 'MQ'], ['fa'], {}, {}, mode='copy')
    with pytest.raises(ValueError, match='fresh output'):
        organize.link_or_copy_genomes([{'name': 'x', 'tier': 'MQ'}], bins, out, ['HQ', 'MQ'], ['fa'], {}, {}, mode='copy')


def test_empty_clipping_is_valid_and_preserves_contigs(tmp_path):
    (tmp_path / 'x-clipping.txt').write_text('contig\tpos\tclipping\n')
    assembly = tmp_path / 'a.fa'
    assembly.write_text('>c\nACGTACGT\n')
    out = tmp_path / 'out.fa'
    filtering.process_assembly(assembly, None, tmp_path, out, 'split', 1, False)
    assert out.read_text() == assembly.read_text()


def test_mag_verdict_cannot_silently_be_used_as_contig_risk(tmp_path):
    p = tmp_path / 'chim.tsv'
    p.write_text('bin\tchimera_verdict\nx\tChimeric\n')
    with pytest.raises(ValueError, match='per-contig'):
        filtering.load_chimera_flagged(p)


def test_yaml_regex_round_trip():
    value = r'^(S\d+)_"x"'
    out = set_yaml_scalar('sample_regex: "x"\n', 'sample_regex', value, quote=True)
    assert json.loads(out.split(': ', 1)[1]) == value


def test_gtbd_clade_suffixes_remain_distinct():
    assert taxa_agree('Firmicutes_A', 'Firmicutes_B') is False
    assert taxa_agree('Firmicutes_A', 'firmicutes_a') is True


def test_circular_header_metadata_and_independent_quality_flag(tmp_path):
    p = tmp_path / 'a.fa'
    p.write_text('>c length=100 circular=true\nACGT\n')
    assert list(read_fasta_streaming(p, full_header=True))[0][0].endswith('circular=true')
    assert filtering.is_circular_candidate('c circular=true', 4, 200000)
    assert not filtering.is_circular_candidate('c circular=false', 500000, 200000)
    categories, _ = filtering.classify_contigs(p, {'c'}, {'c': [2]}, set(), True, 200000)
    assert categories['hq_circular_clip'] == {'c'}
    rows = build_report({'x'}, checkm2={'x': {'completeness': 99., 'contamination': 1.}},
                        organized_tiers={'x': 'MQ'}, structures={'x': {'n_contigs': 1, 'circular_annotated': False}})
    assert rows[0]['tier'] == 'HQ' and rows[0]['organized_tier'] == 'MQ'
    assert rows[0]['near_complete_candidate'] is True
    assert rows[0]['circular_annotated'] is False


def test_read_rank_rollup_clean_control_and_window_output(tmp_path, monkeypatch):
    bam = tmp_path / 'mock.bam'
    bam.write_text('BAM stream is mocked')
    assembly = tmp_path / 'a.fa'
    assembly.write_text('>c\n' + 'ACGT' * 12500 + '\n')
    tax = tmp_path / 'reads.tsv'
    tax.write_text(''.join(f'C\tr{i}\t{11 if i < 10 else 12}\n' for i in range(20)))
    nodes = tmp_path / 'nodes.dmp'
    nodes.write_text('1\t|\t1\t|\tno rank\t|\n10\t|\t1\t|\torder\t|\n'
                     '11\t|\t10\t|\tspecies\t|\n12\t|\t10\t|\tgenus\t|\n')
    monkeypatch.setattr(reads, 'check_samtools', lambda: True)
    monkeypatch.setattr(reads, 'stream_bam', lambda *a, **k: iter((f'r{i}', 'c', 1000 + i) for i in range(20)))
    out = tmp_path / 'out'
    reads.run(args_for(reads, ['check-read-chimeras', '--bam', str(bam), '--assembly', str(assembly),
        '--kraken2-output', str(tax), '--taxonomy-rank', 'order', '--nodes-dmp', str(nodes), '-o', str(out)]))
    with (out / 'chimera_read_report.tsv').open() as f:
        row = next(csv.DictReader(f, delimiter='\t'))
        assert row['n_taxa'] == '1' and row['chimera_risk'] == 'Clean'
    assert (out / 'chimera_read_windows.tsv').exists()


def test_requested_empty_rna_gate_demotes_hq_end_to_end(tmp_path):
    bins = tmp_path / 'bins'
    bins.mkdir()
    (bins / 'x.fa').write_text('>c\nACGT\n')
    qc = tmp_path / 'qc'
    qc.mkdir()
    (qc / 'quality_report.tsv').write_text('Name\tCompleteness\tContamination\nx\t99\t1\n')
    rna = tmp_path / 'rna'
    rna.mkdir()
    out = tmp_path / 'out'
    organize.run(args_for(organize, ['organize-mags', '--checkm2-root', str(qc), '--genomes-dir', str(bins),
        '--mimag-rna-dir', str(rna), '--copy', '-o', str(out)]))
    assert not (out / 'Selected/HQ/x.fa').exists()
    assert (out / 'Selected/MQ/x.fa').exists()
    with (out / 'merged_quality.tsv').open() as f:
        assert next(csv.DictReader(f, delimiter='\t'))['tier'] == 'MQ'


def test_final_report_reads_actual_tier_and_circular_metadata(tmp_path):
    bins = tmp_path / 'bins'
    bins.mkdir()
    (bins / 'x.fa').write_text('>c circular=true\nACGT\n')
    quality = tmp_path / 'quality_report.tsv'
    quality.write_text('Name\tCompleteness\tContamination\nx\t99\t1\n')
    organized = tmp_path / 'merged_quality.tsv'
    organized.write_text('name\ttier\nx\tMQ\n')
    output = tmp_path / 'final.tsv'
    final_report.run(args_for(final_report, ['final-report', '--genomes-dir', str(bins),
        '--checkm2-report', str(quality), '--organized-quality', str(organized), '-o', str(output)]))
    with output.open() as f:
        row = next(csv.DictReader(f, delimiter='\t'))
    assert row['tier'] == 'HQ' and row['organized_tier'] == 'MQ'
    assert row['circular_annotated'] == 'True' and row['near_complete_candidate'] == 'True'


def test_long_contig_processing_keeps_circular_description(tmp_path):
    assembly = tmp_path / 'a.fa'
    assembly.write_text('>c circular=true\n' + 'A' * 400 + '\n')
    (tmp_path / 'x-clipping.txt').write_text('contig\tpos\tclipping\n')
    qc = tmp_path / 'qc'
    qc.mkdir()
    (qc / 'quality_report.tsv').write_text('Name\tCompleteness\tContamination\n'
                                           f'{plc.safe_filename("c")}\t99\t1\n')
    out = tmp_path / 'out'
    plc.run(args_for(plc, ['process-large-contigs', '-a', str(assembly), '-m', str(tmp_path),
        '-o', str(out), '--max-length', '100', '--skip-checkm2', '--checkm2-results', str(qc)]))
    extracted = next((out / 'HQ_extracted').glob('*.fa'))
    assert extracted.read_text().startswith('>c circular=true\n')


def test_cigar_deletion_does_not_invent_aligned_depth(monkeypatch):
    sam_process(monkeypatch, 'r\t0\tc\t1\t60\t100M1000D100M\t*\t0\t0\t*\t*\n')
    result = bam_stats.collect_position_stats(Path('mock.bam'), {'c': [(500, 1)]},
        20, 20, {'c': 1}, {'c': 1200}, 0, 0)[('c', 500)]
    assert result['local_depth'] == 0


@pytest.mark.skipif(shutil.which('samtools') is None, reason='optional real samtools smoke test')
def test_real_bam_depth_and_trailing_clips(tmp_path):
    """Check actual samtools output against a BAM with known base coverage."""
    sam = tmp_path / 'input.sam'
    lines = ['@HD\tVN:1.6\tSO:coordinate\n', '@SQ\tSN:uniform\tLN:20000\n',
             '@SQ\tSN:clipped\tLN:20000\n', '@SQ\tSN:partner\tLN:20000\n']
    for i in range(10):
        lines.append(f'u{i}\t0\tuniform\t1\t60\t20000M\t*\t0\t0\t' + 'ACGT' * 5000 + '\t*\n')
    for i in range(10):
        lines.append(f'c{i}\t0\tclipped\t1001\t60\t9000M100S\t*\t0\t0\t' + 'ACGT' * 2275 + '\t*\tSA:Z:partner,1,+,100M,60,0;\n')
    sam.write_text(''.join(lines))
    bam = tmp_path / 'input.bam'
    subprocess.run(['samtools', 'view', '-b', '-o', str(bam), str(sam)], check=True, capture_output=True)
    mean = bam_stats.get_contig_mean_depths(bam)
    stats = bam_stats.collect_position_stats(bam, {'uniform': [(10000, 1)], 'clipped': [(10000, 1)]},
        500, 20, mean, {'uniform': 20000, 'clipped': 20000}, 0, 0)
    uniform = stats[('uniform', 10000)]
    clipped = stats[('clipped', 10000)]
    assert mean['uniform'] == uniform['local_depth'] == 10
    assert uniform['depth_ratio'] == 1
    assert uniform['discordant_fraction'] == 'N/A'
    assert clipping.classify_position(uniform)[0] == 'low_confidence'
    assert clipped['n_clipped_reads'] == clipped['sa_partner_support'] == 10
    assert clipped['sa_partner_fraction'] == 1
