"""Resolution policies, whole-set replacement validation and metadata failures."""
import argparse
from collections import Counter
import csv
import json
from pathlib import Path

import pytest

from echidna.commands import check_shared_contigs as cli
from echidna.utils.bin_overlap import analyze_bins, discover_bins, read_bin
from echidna.utils.bin_resolution import Inventory, load_quality, load_clusters, plan_resolution, write_resolution

S = 'AACCGT'
A = 'ATAT'
B = 'CCCCAA'
C = 'GGGATA'
D = 'ACACCC'
E = 'AGAGTC'


def fasta(folder, name, sequences):
    folder.mkdir(parents=True, exist_ok=True)
    p = folder / (name + '.fa')
    p.write_text(''.join(f'>{name}_{i} circular=true\n{seq}\n' for i, seq in enumerate(sequences)))
    return p


def table(path, fields, rows, sep='\t'):
    with path.open('w', newline='') as f:
        w = csv.writer(f, delimiter=sep)
        w.writerow(fields)
        w.writerows(rows)
    return path


def read_table(path):
    with path.open() as f:
        return list(csv.DictReader(f, delimiter='\t'))


def setup(tmp_path, genomes, qualities, candidates=None, clusters=None):
    selected = tmp_path / 'selected'
    for name, sequences in genomes.items():
        fasta(selected, name, sequences)
    roots = None
    if candidates:
        folder = tmp_path / 'candidates'
        for name, sequences in candidates.items():
            fasta(folder, name, sequences)
        roots = [folder]
    files = discover_bins([selected])
    before = analyze_bins(files)
    inventory = Inventory(files, before, roots)
    qfile = table(tmp_path / 'quality.tsv', ['Name', 'Completeness', 'Contamination'],
                  [[name, comp, cont] for name, (comp, cont) in qualities.items()])
    quality = load_quality(qfile, inventory)
    membership, eligible, unavailable = {}, set(), []
    if clusters:
        cfile = table(tmp_path / 'Cdb.csv', ['genome', 'secondary_cluster', 'primary_cluster'],
                      [[name + '.fa', cluster, 'same_primary'] for name, cluster in clusters.items()], sep=',')
        membership, eligible, unavailable = load_clusters(cfile, inventory)
    return before, inventory, quality, membership, eligible, unavailable


def plan(data, mode='exclude-mags', rank='completeness'):
    before, inventory, quality, clusters, eligible, _ = data
    return plan_resolution(inventory, quality, mode, rank, clusters, eligible)


def emit(data, resolution, output):
    before, inventory, quality, clusters, _, unavailable = data
    return write_resolution(before, inventory, quality, resolution, output, clusters, unavailable)


def test_exclude_keeps_best_and_saves_originals(tmp_path):
    data = setup(tmp_path, {'high': [S, A], 'low': [S, B], 'clean': [C]},
                 {'high': (99, 1), 'low': (90, 0), 'clean': (80, 0)})
    original = {p: p.read_bytes() for p in data[1].files.values()}
    resolution = plan(data)
    assert set(resolution['retained']) == {'high.fa', 'clean.fa'}
    summary = emit(data, resolution, tmp_path / 'out')
    assert (tmp_path / 'out/excluded_mags/low.fa').is_file()
    assert len(list(read_bin(tmp_path / 'out/retained/high.fa'))) == 2
    assert summary['resolution']['final_n_shared_sequences'] == 0
    assert all(p.read_bytes() == content for p, content in original.items())
    assert all('.echidna-resolution-' not in row['file_path'] for row in read_table(tmp_path / 'out/after/bin_summary.tsv'))


def test_conflict_chain_does_not_collapse_component_to_one_mag(tmp_path):
    data = setup(tmp_path, {'a': [S], 'b': [S, B], 'c': [B]},
                 {'a': (99, 0), 'b': (98, 0), 'c': (97, 0)})
    assert set(plan(data)['retained']) == {'a.fa', 'c.fa'}


def test_trim_saves_removed_contigs_preserves_headers_and_invalidates_quality(tmp_path):
    data = setup(tmp_path, {'high': [S, A], 'low': [S, B]}, {'high': (99, 1), 'low': (90, 1)})
    summary = emit(data, plan(data, 'remove-contigs'), tmp_path / 'out')
    low = list(read_bin(tmp_path / 'out/retained/low.fa', full_header=True))
    assert low == [('low_1 circular=true', B)]
    assert list(read_bin(tmp_path / 'out/removed_contigs/low.fa'))[0][1] == S
    rows = {r['genome']: r for r in read_table(tmp_path / 'out/resolution.tsv')}
    assert rows['low.fa']['quality_status'] == 'requires_reassessment'
    assert rows['low.fa']['final_completeness'] == rows['low.fa']['final_contamination'] == ''
    assert rows['high.fa']['final_completeness'] == '99.000000'
    assert summary['resolution']['quality_reassessment_required']
    assert not read_table(tmp_path / 'out/after/shared_contigs.tsv')


def test_entirely_shared_lower_mag_is_excluded_not_exported_empty(tmp_path):
    data = setup(tmp_path, {'high': [S, A], 'low': [S]}, {'high': (99, 0), 'low': (90, 0)})
    resolution = plan(data, 'remove-contigs')
    assert resolution['plans']['low.fa']['action'] == 'excluded'
    emit(data, resolution, tmp_path / 'out')
    assert not (tmp_path / 'out/retained/low.fa').exists()
    assert (tmp_path / 'out/excluded_mags/low.fa').exists()


def test_completeness_ties_and_optional_contamination_score(tmp_path):
    data = setup(tmp_path, {'a': [S], 'b': [S]}, {'a': (99, 10), 'b': (98, 0)})
    assert set(plan(data)['retained']) == {'a.fa'}
    assert set(plan(data, rank='quality-score')['retained']) == {'b.fa'}
    data[2]['b.fa'] = (99, 0)
    assert set(plan(data)['retained']) == {'b.fa'}
    data[2]['a.fa'] = (99, 0)
    assert set(plan(data)['retained']) == {'a.fa'}  # stable exact-label tie break


def test_best_compatible_cluster_alternative_checked_against_all_retained(tmp_path):
    data = setup(tmp_path, {'high': [S], 'low': [S, B], 'other': [C]},
                 {'high': (100, 0), 'low': (95, 0), 'other': (90, 0), 'best_alt': (99, 0), 'next_alt': (98, 0)},
                 {'best_alt': [C, D], 'next_alt': [B, E]},
                 {'high': '1', 'low': '2', 'other': '3', 'best_alt': '2', 'next_alt': '2'})
    resolution = plan(data)
    assert set(resolution['retained']) == {'high.fa', 'other.fa', 'candidate:next_alt.fa'}
    assert resolution['attempts'][0]['decision'] == 'shared_contigs'
    assert json.loads(resolution['attempts'][0]['blocking_mags']) == ['other.fa']
    assert resolution['attempts'][1]['decision'] == 'accepted'
    summary = emit(data, resolution, tmp_path / 'out')
    assert summary['resolution']['n_replacements'] == 1
    assert summary['resolution']['final_n_shared_sequences'] == 0


def test_replacements_cannot_share_with_previous_replacements(tmp_path):
    data = setup(tmp_path, {'high': [S], 'low1': [S, A], 'low2': [S, B]},
                 {'high': (100, 0), 'low1': (95, 0), 'low2': (90, 0), 'alt1': (96, 0), 'alt2': (94, 0)},
                 {'alt1': [D], 'alt2': [D]},
                 {'high': '1', 'low1': '2', 'low2': '3', 'alt1': '2', 'alt2': '3'})
    resolution = plan(data)
    assert set(resolution['retained']) == {'high.fa', 'candidate:alt1.fa'}
    assert resolution['unresolved'][0]['original'] == 'low2.fa'
    assert json.loads(resolution['attempts'][1]['blocking_mags']) == ['candidate:alt1.fa']


def test_trim_can_be_replaced_by_whole_compatible_alternative(tmp_path):
    data = setup(tmp_path, {'high': [S], 'low': [S, B]},
                 {'high': (100, 0), 'low': (95, 0), 'alt': (94, 0)},
                 {'alt': [B, C]}, {'high': '1', 'low': '2', 'alt': '2'})
    resolution = plan(data, 'remove-contigs')
    assert resolution['plans']['low.fa']['action'] == 'replaced'
    summary = emit(data, resolution, tmp_path / 'out')
    assert not summary['resolution']['quality_reassessment_required']
    assert not list((tmp_path / 'out/removed_contigs').iterdir())
    assert len(list(read_bin(tmp_path / 'out/retained/alt.fa'))) == 2


def test_no_compatible_alternative_leaves_trimmed_mag_marked_for_qc(tmp_path):
    data = setup(tmp_path, {'high': [S], 'low': [S, B]},
                 {'high': (100, 0), 'low': (95, 0), 'alt': (94, 0)},
                 {'alt': [S, C]}, {'high': '1', 'low': '2', 'alt': '2'})
    resolution = plan(data, 'remove-contigs')
    assert resolution['plans']['low.fa']['action'] == 'trimmed'
    assert resolution['unresolved'][0]['reason'] == 'no_compatible_candidate'


def test_selected_copy_in_candidate_pool_is_not_a_second_candidate(tmp_path):
    data = setup(tmp_path, {'high': [S], 'low': [S, B]},
                 {'high': (100, 0), 'low': (95, 0)},
                 {'high': [S], 'low': [S, B]}, {'high': '1', 'low': '2'})
    assert data[4] == set()
    resolution = plan(data)
    assert resolution['unresolved'][0]['reason'] == 'no_available_candidate'


def test_missing_fasta_members_are_reported(tmp_path):
    data = setup(tmp_path, {'a': [S]}, {'a': (90, 0)}, clusters={'a': '1', 'absent': '1'})
    assert data[5] == [{'genome': 'absent.fa', 'cluster': '1'}]


@pytest.mark.parametrize('values', [('nan', '1'), ('95', 'inf'), ('101', '1'), ('95', '-1'), ('?', '0')])
def test_bad_quality_is_rejected(tmp_path, values):
    data = setup(tmp_path, {'a': [S]}, {'a': (90, 0)})
    path = table(tmp_path / 'bad.tsv', ['genome', 'completeness', 'contamination'], [['a', *values]])
    with pytest.raises(ValueError, match='Invalid quality'):
        load_quality(path, data[1])


def test_missing_and_duplicate_quality_is_rejected(tmp_path):
    data = setup(tmp_path, {'a': [S], 'b': [B]}, {'a': (90, 0)})
    with pytest.raises(ValueError, match='Quality missing'):
        plan(data)
    path = table(tmp_path / 'dup.tsv', ['genome', 'completeness', 'contamination'],
                 [['a', 90, 0], ['a.fa', 90, 0]])
    with pytest.raises(ValueError, match='Duplicate quality'):
        load_quality(path, data[1])


def test_primary_clusters_are_not_used_and_conflicting_memberships_fail(tmp_path):
    data = setup(tmp_path, {'a': [S]}, {'a': (90, 0)})
    path = table(tmp_path / 'bad.csv', ['genome', 'primary_cluster'], [['a.fa', '1']], sep=',')
    with pytest.raises(ValueError, match='primary_cluster alone'):
        load_clusters(path, data[1])
    path = table(path, ['genome', 'cluster'], [['a.fa', '1'], ['a.fa', '2']], sep=',')
    with pytest.raises(ValueError, match='multiple clusters'):
        load_clusters(path, data[1])


def test_cluster_replacement_expects_one_selected_rep_per_cluster(tmp_path):
    data = setup(tmp_path, {'a': [S], 'b': [B]}, {'a': (90, 0), 'b': (80, 0)})
    path = table(tmp_path / 'clusters.tsv', ['genome', 'cluster'], [['a', '1'], ['b', '1']])
    with pytest.raises(ValueError, match='Multiple selected MAGs'):
        load_clusters(path, data[1])


def test_input_change_aborts_before_publishing_output(tmp_path):
    data = setup(tmp_path, {'a': [S]}, {'a': (90, 0)})
    resolution = plan(data)
    data[1].files['a.fa'].write_text('>changed\nCCCCCC\n')
    with pytest.raises(ValueError, match='changed during resolution'):
        emit(data, resolution, tmp_path / 'out')
    assert not (tmp_path / 'out').exists()
    assert not list(tmp_path.glob('.echidna-resolution-*'))


def test_cli_drep_and_required_resolution_flags(tmp_path):
    data = setup(tmp_path, {'high': [S], 'low': [S, B]},
                 {'high': (100, 0), 'low': (95, 0), 'alt': (94, 0)},
                 {'alt': [C]}, {'high': '1', 'low': '2', 'alt': '2'})
    drep = tmp_path / 'drep/data_tables'
    drep.mkdir(parents=True)
    (drep / 'Cdb.csv').write_text((tmp_path / 'Cdb.csv').read_text())
    parser = argparse.ArgumentParser()
    cli.register_parser(parser.add_subparsers())
    base = ['check-shared-contigs', '--bins-dir', str(tmp_path / 'selected'), '-o', str(tmp_path / 'out')]
    args = parser.parse_args(base + ['--exclude-mags', '--quality-report', str(tmp_path / 'quality.tsv'),
                                    '--drep-dir', str(drep.parent), '--candidates-dir', str(tmp_path / 'candidates')])
    cli.run(args)
    assert (tmp_path / 'out/retained/alt.fa').is_file()
    with pytest.raises(ValueError, match='requires --quality-report'):
        cli.run(parser.parse_args(base[:-1] + [str(tmp_path / 'out2'), '--exclude-mags']))
    with pytest.raises(ValueError, match='require --exclude-mags'):
        cli.run(parser.parse_args(base[:-1] + [str(tmp_path / 'out3'), '--quality-report', str(tmp_path / 'quality.tsv')]))
    with pytest.raises(SystemExit):
        parser.parse_args(base + ['--exclude-mags', '--remove-shared-contigs'])


def test_pending_cluster_is_retried_after_another_replacement(tmp_path):
    data = setup(tmp_path, {'high': [S], 'mid': [S, B], 'low': [S, C]},
                 {'high': (100, 0), 'mid': (95, 0), 'low': (90, 0), 'alt_mid': (94, 0), 'alt_low': (89, 0)},
                 {'alt_mid': [C], 'alt_low': [D]},
                 {'high': '1', 'mid': '2', 'low': '3', 'alt_mid': '2', 'alt_low': '3'})
    resolution = plan(data, 'remove-contigs')
    assert set(resolution['retained']) == {'high.fa', 'candidate:alt_mid.fa', 'candidate:alt_low.fa'}
    assert any(a['candidate'] == 'candidate:alt_mid.fa' and a['round'] == 2 and a['decision'] == 'accepted'
               for a in resolution['attempts'])
    assert not resolution['unresolved']


def test_original_mag_is_restored_when_replacement_frees_its_contigs(tmp_path):
    data = setup(tmp_path, {'high': [S], 'mid': [S, B], 'low': [B, C]},
                 {'high': (100, 0), 'mid': (95, 0), 'low': (90, 0), 'alt_mid': (94, 0)},
                 {'alt_mid': [D]}, {'high': '1', 'mid': '2', 'low': '3', 'alt_mid': '2'})
    resolution = plan(data, 'remove-contigs')
    assert resolution['plans']['low.fa']['action'] == 'restored'
    summary = emit(data, resolution, tmp_path / 'out')
    assert summary['resolution']['n_restored_originals'] == 1
    assert not summary['resolution']['quality_reassessment_required']
    assert len(list(read_bin(tmp_path / 'out/retained/low.fa'))) == 2


def test_clusters_without_overlap_are_unchanged(tmp_path):
    data = setup(tmp_path, {'a': [S], 'b': [B]}, {'a': (90, 0), 'b': (80, 0)},
                 clusters={'a': '1', 'b': '2'})
    resolution = plan(data)
    assert not resolution['attempts'] and not resolution['unresolved']
    assert set(resolution['retained']) == {'a.fa', 'b.fa'}
    emit(data, resolution, tmp_path / 'out')


def test_missing_quality_for_available_alternative_is_an_error(tmp_path):
    data = setup(tmp_path, {'high': [S], 'low': [S, B]}, {'high': (100, 0), 'low': (95, 0)},
                 {'alt': [C]}, {'high': '1', 'low': '2', 'alt': '2'})
    with pytest.raises(ValueError, match='Quality missing'):
        plan(data)


def test_ambiguous_names_require_exact_labels(tmp_path):
    first, second = tmp_path / 'first', tmp_path / 'second'
    fasta(first, 'a', [S])
    fasta(second, 'a', [B])
    files = discover_bins([first, second])
    inventory = Inventory(files, analyze_bins(files))
    path = table(tmp_path / 'quality.tsv', ['genome', 'completeness', 'contamination'], [['a', 90, 0]])
    with pytest.raises(ValueError, match='Ambiguous genome'):
        load_quality(path, inventory)
    table(path, ['genome', 'completeness', 'contamination'], [['set1/a.fa', 90, 0], ['set2/a.fa', 80, 0]])
    assert len(load_quality(path, inventory)) == 2


def test_export_recheck_rejects_a_corrupt_plan(tmp_path):
    data = setup(tmp_path, {'a': [S], 'b': [S, B]}, {'a': (90, 0), 'b': (80, 0)})
    resolution = plan(data)
    resolution['retained']['b.fa'] = set(data[1].counts('b.fa'))
    resolution['plans']['b.fa']['action'] = 'retained'
    with pytest.raises(RuntimeError, match='Exported MAGs still share'):
        emit(data, resolution, tmp_path / 'out')
    assert not (tmp_path / 'out').exists()


def test_output_names_cannot_collide_with_generated_prefixes(tmp_path):
    first, second = tmp_path / 'first', tmp_path / 'second'
    fasta(first, 'a', [S])
    fasta(first, '000002__a', [D])
    fasta(second, 'a', [B])
    files = discover_bins([first, second])
    before = analyze_bins(files)
    inventory = Inventory(files, before)
    quality = {key: (90, 0) for key in files}
    resolution = plan_resolution(inventory, quality, 'exclude-mags')
    write_resolution(before, inventory, quality, resolution, tmp_path / 'out')
    assert len(list((tmp_path / 'out/retained').glob('*.fa'))) == 3
    assert len({row['output_file'] for row in read_table(tmp_path / 'out/resolution.tsv')}) == 3


def test_generic_cluster_table_and_empty_existing_output(tmp_path):
    data = setup(tmp_path, {'a': [S]}, {'a': (90, 0)})
    path = table(tmp_path / 'clusters.tsv', ['genome', 'cluster'], [['a.fa', 'group1']])
    clusters, _, _ = load_clusters(path, data[1])
    assert clusters == {'a.fa': 'group1'}
    out = tmp_path / 'out'
    out.mkdir()
    emit(data, plan(data), out)
    assert (out / 'retained/a.fa').is_file()
    with pytest.raises(ValueError, match='new or empty'):
        emit(data, plan(data), out)


def test_contigs_freed_by_replacement_are_restored_to_trimmed_mag(tmp_path):
    data = setup(tmp_path, {'high': [S], 'mid': [S, B], 'low': [S, B, C]},
                 {'high': (100, 0), 'mid': (95, 0), 'low': (90, 0), 'alt_mid': (94, 0)},
                 {'alt_mid': [D]}, {'high': '1', 'mid': '2', 'low': '3', 'alt_mid': '2'})
    resolution = plan(data, 'remove-contigs')
    assert resolution['retained']['low.fa'] == set(data[1].read_counts(
        fasta(tmp_path / 'expected', 'low', [B, C])))
    assert resolution['plans']['low.fa']['removed'] == set(data[1].read_counts(
        fasta(tmp_path / 'expected', 'removed', [S])))
    emit(data, resolution, tmp_path / 'out')
    assert [seq for _, seq in read_bin(tmp_path / 'out/retained/low.fa')] == [B, C]
    assert [seq for _, seq in read_bin(tmp_path / 'out/removed_contigs/low.fa')] == [S]


def test_initially_empty_remnant_can_recover_freed_contigs(tmp_path):
    data = setup(tmp_path, {'high': [S], 'mid': [S, B], 'low': [S, B]},
                 {'high': (100, 0), 'mid': (95, 0), 'low': (90, 0), 'alt_mid': (94, 0)},
                 {'alt_mid': [D]}, {'high': '1', 'mid': '2', 'low': '3', 'alt_mid': '2'})
    resolution = plan(data, 'remove-contigs')
    assert resolution['plans']['low.fa']['action'] == 'trimmed'
    emit(data, resolution, tmp_path / 'out')
    assert [seq for _, seq in read_bin(tmp_path / 'out/retained/low.fa')] == [B]
    assert resolution['unresolved'][0]['final_action'] == 'trimmed'
