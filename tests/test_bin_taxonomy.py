"""Sequence-linked taxonomy reuse, selection gates and final FASTA checks."""
import argparse
import csv
import json

import pytest

from echidna.commands import classify_contigs, check_shared_contigs, run_mmseqs_contigs
from echidna.utils.bin_overlap import read_bin, sequence_digest
from echidna.utils.bin_taxonomy import TaxonomyEvidence, BinTaxonomy
from echidna.utils.bin_resolution import Inventory, plan_resolution, write_resolution
from echidna.utils.bin_overlap import analyze_bins, discover_bins

S, A, B, C, D, E = 'AACCGT', 'ATAT', 'CCCCAA', 'GGGATA', 'ACACCC', 'AGAGTC'


def taxonomy_table(path, entries, rank='order'):
    with path.open('w', newline='') as stream:
        writer = csv.writer(stream, delimiter='\t')
        writer.writerow(['contig', 'rank_taxon', 'sequence_sha256', 'contig_length', 'taxonomy_rank'])
        for i, (seq, taxon) in enumerate(entries):
            writer.writerow([f'old_{i}', taxon, sequence_digest(seq), len(seq), rank])
    return path


def fasta(root, name, sequences):
    root.mkdir(parents=True, exist_ok=True)
    path = root / (name + '.fa')
    path.write_text(''.join(f'>{name}_{i}\n{seq}\n' for i, seq in enumerate(sequences)))
    return path


def parser_for(command):
    parser = argparse.ArgumentParser()
    command.register_parser(parser.add_subparsers())
    return parser


def data(tmp_path, genomes, entries, candidates=None, **thresholds):
    for name, seqs in genomes.items():
        fasta(tmp_path / 'selected', name, seqs)
    for name, seqs in (candidates or {}).items():
        fasta(tmp_path / 'candidates', name, seqs)
    files = discover_bins([tmp_path / 'selected'])
    before = analyze_bins(files)
    inventory = Inventory(files, before, [tmp_path / 'candidates'] if candidates else None)
    evidence = TaxonomyEvidence(taxonomy_table(tmp_path / 'tax.tsv', entries), **thresholds)
    return before, inventory, BinTaxonomy(evidence, inventory)


def test_sequence_match_survives_renaming_and_reverse_complement(tmp_path):
    evidence = TaxonomyEvidence(taxonomy_table(tmp_path / 'tax.tsv', [(S, 'OrderA')]))
    path = tmp_path / 'renamed.fa'
    path.write_text('>entirely_new_name\nACGGTT\n')
    assert evidence.assess_file(path)['taxonomy_status'] == 'coherent'


def test_zero_and_sparse_taxonomy_are_unassessed_not_coherent(tmp_path):
    evidence = TaxonomyEvidence(taxonomy_table(tmp_path / 'tax.tsv', [(A, 'OrderA')]))
    missing = evidence.assess_file(fasta(tmp_path, 'missing', [S]))
    sparse = evidence.assess_file(fasta(tmp_path, 'sparse', [A, S]))
    assert missing['taxonomy_status'] == sparse['taxonomy_status'] == 'insufficient_taxonomy'
    assert missing['taxonomy_foreign_fraction'] is None
    assert sparse['taxonomy_resolved_fraction_bp'] == 0.4


def test_foreign_fraction_uses_resolved_bp_but_coverage_uses_all_bp(tmp_path):
    evidence = TaxonomyEvidence(taxonomy_table(tmp_path / 'tax.tsv', [(A, 'OrderA'), (S, 'OrderB')]))
    result = evidence.assess_file(fasta(tmp_path, 'mixed', [A, S, B]))
    assert result['taxonomy_status'] == 'incoherent_taxonomy'
    assert result['taxonomy_foreign_fraction'] == 0.4
    assert result['taxonomy_resolved_fraction_bp'] == 10 / 16


def test_conflicting_same_sequence_assignments_are_unresolved(tmp_path):
    evidence = TaxonomyEvidence(taxonomy_table(tmp_path / 'tax.tsv', [(S, 'OrderA'), (S, 'OrderB')]))
    result = evidence.assess_file(fasta(tmp_path, 'conflict', [S]))
    assert result['taxonomy_status'] == 'insufficient_taxonomy'
    assert result['taxonomy_conflicting_bp'] == len(S)


@pytest.mark.parametrize('option,value', [('min_resolved_frac', 0), ('min_resolved_frac', float('nan')),
                                         ('min_foreign_frac', -1), ('min_foreign_frac', float('inf'))])
def test_invalid_thresholds_fail(tmp_path, option, value):
    with pytest.raises(ValueError):
        TaxonomyEvidence(taxonomy_table(tmp_path / 'tax.tsv', [(S, 'OrderA')]), **{option: value})


def test_reject_legacy_bin_summary_and_rank_mismatch(tmp_path):
    old = tmp_path / 'old.tsv'
    old.write_text('contig\tbin\trank_taxon\na\tb\tOrderA\n')
    with pytest.raises(ValueError, match='sequence-linked'):
        TaxonomyEvidence(old)
    rich = taxonomy_table(tmp_path / 'rich.tsv', [(S, 'OrderA')])
    with pytest.raises(ValueError, match='rank'):
        TaxonomyEvidence(rich, rank='genus')


@pytest.mark.parametrize('mode', ['exclude-mags', 'remove-contigs'])
def test_skip_complete_but_incoherent_or_unclassified_replacements(tmp_path, mode):
    before, inv, tax = data(tmp_path, {'high': [S, A], 'original': [S, B]},
                           [(seq, 'OrderA') for seq in (S, A, B, D)] + [(C, 'OrderB')],
                           {'overlap': [S, D], 'bad': [B, C], 'unknown': [E], 'good': [D]})
    quality = {'high.fa': (100, 0), 'original.fa': (90, 0), 'candidate:bad.fa': (99, 0),
               'candidate:unknown.fa': (98, 0), 'candidate:good.fa': (95, 0), 'candidate:overlap.fa': (99.5, 0)}
    clusters = {'high.fa': '1', 'original.fa': '2', **{k: '2' for k in quality if k.startswith('candidate:')}}
    plan = plan_resolution(inv, quality, mode, clusters=clusters,
                           eligible=set(quality) - inv.primary, taxonomy=tax)
    assert plan['plans']['original.fa']['replaced_by'] == 'candidate:good.fa'
    decisions = {row['candidate']: row['decision'] for row in plan['attempts']}
    assert decisions['candidate:overlap.fa'] == 'shared_contigs'
    assert decisions['candidate:bad.fa'] == 'incoherent_taxonomy'
    assert decisions['candidate:unknown.fa'] == 'insufficient_taxonomy'
    output = tmp_path / 'out'
    summary = write_resolution(before, inv, quality, plan, output, clusters)
    assert summary['resolution']['final_taxonomy_check'] == 'passed'
    assert summary['resolution']['taxonomy_retained_mags_checked'] == 2
    assert (output / 'taxonomy_checks.tsv').exists()
    assert not (output / 'retained/bad.fa').exists()


def test_reject_incoherent_primary_without_overlap_and_recover_cluster(tmp_path):
    before, inv, tax = data(tmp_path, {'original': [B, C]}, [(B, 'A'), (C, 'B'), (D, 'A')], {'alt': [D]})
    quality = {'original.fa': (99, 0), 'candidate:alt.fa': (95, 0)}
    plan = plan_resolution(inv, quality, 'exclude-mags', clusters={key: '1' for key in quality},
                           eligible={'candidate:alt.fa'}, taxonomy=tax)
    assert plan['plans']['original.fa']['action'] == 'replaced'
    assert plan['plans']['original.fa']['exclusion_reason'] == 'incoherent_taxonomy'


def test_no_passing_mags_still_writes_truthful_reports(tmp_path):
    before, inv, tax = data(tmp_path, {'original': [B, C]}, [(B, 'A'), (C, 'B')])
    quality = {'original.fa': (99, 0)}
    plan = plan_resolution(inv, quality, 'exclude-mags', taxonomy=tax)
    summary = write_resolution(before, inv, quality, plan, tmp_path / 'out')
    assert summary['resolution']['n_retained_mags'] == 0
    assert summary['resolution']['taxonomy_retained_mags_checked'] == 0


def test_contig_recovery_cannot_reintroduce_taxonomic_incoherence(tmp_path):
    before, inv, tax = data(tmp_path, {'high': [S], 'mid': [S, B], 'low': [S, B, C]},
                           [(S, 'A'), (B, 'A'), (C, 'B'), (D, 'A')], {'alt': [D]})
    quality = {'high.fa': (100, 0), 'mid.fa': (95, 0), 'low.fa': (90, 0), 'candidate:alt.fa': (94, 0)}
    clusters = {'high.fa': '1', 'mid.fa': '2', 'low.fa': '3', 'candidate:alt.fa': '2'}
    plan = plan_resolution(inv, quality, 'remove-contigs', clusters=clusters,
                           eligible={'candidate:alt.fa'}, taxonomy=tax)
    assert plan['retained']['low.fa'] == {sequence_digest(C)}
    assert any(row['stage'] == 'contig_recovery' and row['taxonomy_status'] == 'incoherent_taxonomy'
               for row in tax.checks)
    write_resolution(before, inv, quality, plan, tmp_path / 'out', clusters)


def test_actual_export_taxonomy_checked_even_if_plan_corrupted(tmp_path):
    before, inv, tax = data(tmp_path, {'original': [B, C]}, [(B, 'A'), (C, 'B')])
    quality = {'original.fa': (99, 0)}
    plan = plan_resolution(inv, quality, 'exclude-mags', taxonomy=tax)
    plan['retained']['original.fa'] = set(inv.counts('original.fa'))
    plan['plans']['original.fa']['action'] = 'retained'
    with pytest.raises(RuntimeError, match='taxonomy recheck'):
        write_resolution(before, inv, quality, plan, tmp_path / 'out')
    assert not (tmp_path / 'out').exists()


def test_classify_output_reused_on_renamed_contigs_and_new_bin_membership(tmp_path):
    bins = tmp_path / 'bins'
    fasta(bins, 'source', [A, S])
    mm = tmp_path / 'mm.tsv'
    mm.write_text('source_0\t1\torder\tA\td_Bacteria;o_OrderA;g_X\n'
                  'source_1\t2\torder\tB\td_Bacteria;o_OrderB;g_Y\n')
    output = tmp_path / 'classified.tsv'
    parser = parser_for(classify_contigs)
    classify_contigs.run(parser.parse_args(['classify-contigs', '--bins-dir', str(bins),
                         '--mmseqs-output', str(mm), '-o', str(output)]))
    assert TaxonomyEvidence(output).rank == 'order'
    final = tmp_path / 'final'
    fasta(final, 'new', [A])
    reused = tmp_path / 'reused.tsv'
    classify_contigs.run(parser.parse_args(['classify-contigs', '--bins-dir', str(final),
                         '--taxonomy-table', str(output), '-o', str(reused)]))
    with reused.open() as stream:
        row = next(csv.DictReader(stream, delimiter='\t'))
    assert row['contig'] == 'new_0' and row['rank_taxon'] == 'OrderA' and row['is_foreign'] == 'False'
    with (tmp_path / 'reused.tsv.bins.tsv').open() as stream:
        assert next(csv.DictReader(stream, delimiter='\t'))['foreign_frac'] == '0.0000'


def test_id_collisions_rejected_before_classification_and_shared_ids_deduplicated(tmp_path):
    bins = tmp_path / 'bins'
    bins.mkdir()
    (bins / 'a.fa').write_text('>same\nAACCGT\n')
    (bins / 'b.fa').write_text('>same\nACGGTT\n')
    assert run_mmseqs_contigs.concat_bins(bins, ['fa'], tmp_path / 'concat.fa') == 1
    (bins / 'b.fa').write_text('>same\nATAT\n')
    with pytest.raises(ValueError, match='different sequences'):
        run_mmseqs_contigs.concat_bins(bins, ['fa'], tmp_path / 'concat.fa')


def test_reuse_handles_same_id_with_different_sequences_across_bins(tmp_path):
    bins = tmp_path / 'bins'
    bins.mkdir()
    (bins / 'a.fa').write_text('>same\nAACCGT\n')
    (bins / 'b.fa').write_text('>same\nATAT\n')
    table = taxonomy_table(tmp_path / 'table.tsv', [(S, 'A'), (A, 'B')])
    out = tmp_path / 'out.tsv'
    parser = parser_for(classify_contigs)
    classify_contigs.run(parser.parse_args(['classify-contigs', '--bins-dir', str(bins),
                         '--taxonomy-table', str(table), '-o', str(out)]))
    with out.open() as stream:
        assert {row['bin']: row['rank_taxon'] for row in csv.DictReader(stream, delimiter='\t')} == {'a': 'A', 'b': 'B'}


def test_require_retained_stops_automation_but_preserves_diagnostics(tmp_path):
    bins = tmp_path / 'bins'
    fasta(bins, 'original', [B, C])
    taxonomy = taxonomy_table(tmp_path / 'tax.tsv', [(B, 'A'), (C, 'B')])
    quality = tmp_path / 'quality.tsv'
    quality.write_text('Name\tCompleteness\tContamination\noriginal\t99\t0\n')
    out = tmp_path / 'out'
    args = parser_for(check_shared_contigs).parse_args(['check-shared-contigs', '--bins-dir', str(bins),
               '--exclude-mags', '--require-retained', '--quality-report', str(quality),
               '--contig-taxonomy', str(taxonomy), '-o', str(out)])
    with pytest.raises(ValueError, match='No MAGs passed'):
        check_shared_contigs.run(args)
    assert (out / 'taxonomy_checks.tsv').exists()
    assert json.loads((out / 'summary.json').read_text())['resolution']['final_taxonomy_check'] == 'no_retained_mags'
