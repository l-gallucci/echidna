"""Tests for echidna.utils.contig_taxonomy and the classify-contigs command.
All pure-function tests; no classifier is executed."""

import argparse
from pathlib import Path

import pytest

from echidna.utils.contig_taxonomy import (
    parse_nodes_dmp,
    parse_mmseqs_taxonomy,
    parse_mmseqs_lineage,
    extract_rank_from_lineage,
    resolve_at_rank,
    assign_contigs_to_bins,
    summarize_bin_taxa,
    summarize_bin_lineage,
)
from echidna.commands import classify_contigs, run_mmseqs_contigs


# A tiny taxonomy tree:
#   root(1) -> Bacteria(2, superkingdom) -> Proteobacteria(10, phylum)
#     -> class(20) -> OrderA(30) -> family(40) -> genus(50) -> speciesA(60)
#                  -> OrderB(31) -> speciesB(61)
NODES_DMP = """\
1\t|\t1\t|\tno rank\t|
2\t|\t1\t|\tsuperkingdom\t|
10\t|\t2\t|\tphylum\t|
20\t|\t10\t|\tclass\t|
30\t|\t20\t|\torder\t|
40\t|\t30\t|\tfamily\t|
50\t|\t40\t|\tgenus\t|
60\t|\t50\t|\tspecies\t|
31\t|\t20\t|\torder\t|
61\t|\t31\t|\tspecies\t|
"""

NAMES_DMP = """\
30\t|\tOrderA\t|\t\t|\tscientific name\t|
31\t|\tOrderB\t|\t\t|\tscientific name\t|
60\t|\tSpeciesA\t|\t\t|\tscientific name\t|
61\t|\tSpeciesB\t|\t\t|\tscientific name\t|
"""


@pytest.fixture
def nodes(tmp_path):
    f = tmp_path / "nodes.dmp"
    f.write_text(NODES_DMP)
    return parse_nodes_dmp(f)


def test_parse_nodes_dmp(nodes):
    assert nodes[60] == (50, 'species')
    assert nodes[30] == (20, 'order')
    assert nodes[2] == (1, 'superkingdom')


def test_resolve_at_rank_walks_up_to_order(nodes):
    assert resolve_at_rank(60, 'order', nodes) == 30     # speciesA -> OrderA
    assert resolve_at_rank(61, 'order', nodes) == 31     # speciesB -> OrderB


def test_resolve_at_rank_above_rank_is_none(nodes):
    # a phylum-level taxid cannot be resolved DOWN to order
    assert resolve_at_rank(10, 'order', nodes) is None


def test_resolve_at_rank_domain_alias(nodes):
    # 'domain' (GTDB) maps to NCBI 'superkingdom'
    assert resolve_at_rank(60, 'domain', nodes) == 2


def test_resolve_at_rank_unclassified_and_unknown(nodes):
    assert resolve_at_rank(0, 'order', nodes) is None
    assert resolve_at_rank(99999, 'order', nodes) is None


def test_assign_contigs_to_bins(tmp_path):
    (tmp_path / "binA.fa").write_text(">c1 foo\nACGT\n>c2\nTTTT\n")
    (tmp_path / "binB.fa").write_text(">c3\nGGGG\n")
    bins = assign_contigs_to_bins(tmp_path, ['fa'])
    assert bins['binA'] == ['c1', 'c2']      # header first word only
    assert bins['binB'] == ['c3']


def test_summarize_bin_taxa_flags_minority_contig(nodes):
    bin_contigs = {
        'chimera_bin': ['c1', 'c2', 'c3'],   # 2x OrderA + 1x OrderB -> foreign c3
        'clean_bin': ['d1', 'd2'],           # both OrderA
    }
    contig_taxid = {'c1': 60, 'c2': 60, 'c3': 61, 'd1': 60, 'd2': 60}
    out = summarize_bin_taxa(bin_contigs, contig_taxid, 'order', nodes)

    chim = out['chimera_bin']
    assert chim['chimera'] is True
    assert chim['n_taxa'] == 2
    assert chim['majority_taxon'] == 30                  # OrderA is majority (2 contigs)
    assert [c for c, _ in chim['foreign_contigs']] == ['c3']

    clean = out['clean_bin']
    assert clean['chimera'] is False
    assert clean['n_taxa'] == 1
    assert clean['foreign_contigs'] == []


def test_summarize_bin_taxa_unresolved_not_counted(nodes):
    # a contig classified only to phylum cannot resolve at order -> not a split
    bin_contigs = {'b': ['c1', 'c2']}
    contig_taxid = {'c1': 60, 'c2': 10}      # c2 is phylum-level
    out = summarize_bin_taxa(bin_contigs, contig_taxid, 'order', nodes)
    assert out['b']['n_taxa'] == 1
    assert out['b']['n_contigs_resolved'] == 1
    assert out['b']['chimera'] is False


def test_summarize_with_names(tmp_path, nodes):
    names_file = tmp_path / "names.dmp"
    names_file.write_text(NAMES_DMP)
    from echidna.utils.read_taxonomy import parse_names_dmp
    names = parse_names_dmp(names_file)
    out = summarize_bin_taxa({'b': ['c1', 'c2']}, {'c1': 60, 'c2': 61},
                             'order', nodes, names)
    assert out['b']['majority_name'] in ('OrderA', 'OrderB')
    assert set(out['b']['taxa']) == {'OrderA', 'OrderB'}


def test_build_rows_marks_foreign(nodes):
    summary = summarize_bin_taxa(
        {'b': ['c1', 'c2', 'c3']}, {'c1': 60, 'c2': 60, 'c3': 61}, 'order', nodes)
    contig_rows, bin_rows = classify_contigs.build_rows(summary, {'c1': 60, 'c2': 60, 'c3': 61})
    foreign = {r['contig']: r['is_foreign'] for r in contig_rows}
    assert foreign == {'c1': False, 'c2': False, 'c3': True}
    assert bin_rows[0]['chimera_flag'] is True
    assert bin_rows[0]['n_foreign_contigs'] == 1


def test_classify_contigs_rank_default_and_mutually_exclusive():
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers()
    classify_contigs.register_parser(sub)
    args = parser.parse_args(['classify-contigs', '--bins-dir', 'bins',
                              '--kraken2-output', 'c.k2',
                              '--nodes-dmp', 'nodes.dmp'])
    assert args.rank == 'order'
    assert args.kaiju_output is None and args.mmseqs_output is None
    # MMseqs2 output is accepted as a classifier
    a2 = parser.parse_args(['classify-contigs', '--bins-dir', 'bins',
                            '--mmseqs-output', 'c.mmseqs.tsv',
                            '--nodes-dmp', 'nodes.dmp'])
    assert a2.mmseqs_output is not None
    # two classifiers at once is rejected
    with pytest.raises(SystemExit):
        parser.parse_args(['classify-contigs', '--bins-dir', 'bins',
                           '--kraken2-output', 'c.k2', '--kaiju-output', 'c.kj',
                           '--nodes-dmp', 'nodes.dmp'])


def test_parse_mmseqs_taxonomy(tmp_path):
    f = tmp_path / "tax.tsv"
    f.write_text("c1\t60\tspecies\tSpeciesA\n"
                 "c2\t61\tspecies\tSpeciesB\n"
                 "c3\t0\tno rank\tunclassified\n")
    out = parse_mmseqs_taxonomy(f)
    assert out == {'c1': 60, 'c2': 61, 'c3': 0}


def test_mmseqs_build_commands_order():
    steps = run_mmseqs_contigs.build_commands(
        'mmseqs', Path('w/contigs.fa'), Path('/db/gtdb'),
        Path('out.tsv'), Path('w'), 8, [])
    labels = [s[0] for s in steps]
    assert labels == ['createdb', 'taxonomy', 'createtsv']
    # taxonomy step references the target DB, threads and emits the lineage column
    tax = steps[1][1]
    assert '/db/gtdb' in tax and '--threads' in tax and '8' in tax
    assert '--tax-lineage' in tax and '1' in tax
    # createtsv writes the final output
    assert str(Path('out.tsv')) in steps[2][1]


def test_mmseqs_build_commands_user_lineage_not_duplicated():
    steps = run_mmseqs_contigs.build_commands(
        'mmseqs', Path('w/c.fa'), Path('/db'), Path('o.tsv'), Path('w'), 4,
        ['--tax-lineage', '2'])
    tax = steps[1][1]
    assert tax.count('--tax-lineage') == 1       # user's value kept, not doubled


# ---- lineage mode (MMseqs2 --tax-lineage / GTDB; no nodes.dmp) -------------

def test_extract_rank_from_lineage_single_and_double_underscore():
    mm = "d_Bacteria;p_Pseudomonadota;c_Gammaproteobacteria;o_Enterobacterales"
    gtdb = "d__Bacteria;p__Pseudomonadota;c__Alphaproteobacteria"
    assert extract_rank_from_lineage(mm, 'class') == 'Gammaproteobacteria'
    assert extract_rank_from_lineage(mm, 'order') == 'Enterobacterales'
    assert extract_rank_from_lineage(gtdb, 'class') == 'Alphaproteobacteria'
    # rank absent from lineage -> None (contig only classified higher)
    assert extract_rank_from_lineage(mm, 'genus') is None
    assert extract_rank_from_lineage("d_Bacteria", 'class') is None


def test_parse_mmseqs_lineage(tmp_path):
    f = tmp_path / "tax.tsv"
    f.write_text(
        "c1\t1224\tclass\tGammaproteobacteria\td_Bacteria;p_Pseudomonadota;c_Gammaproteobacteria\n"
        "c2\t28211\tclass\tAlphaproteobacteria\td_Bacteria;p_Pseudomonadota;c_Alphaproteobacteria\n"
        "c3\t0\tno rank\tunclassified\n")          # no lineage column -> omitted
    out = parse_mmseqs_lineage(f)
    assert set(out) == {'c1', 'c2'}
    assert out['c1'].endswith('c_Gammaproteobacteria')


def test_summarize_bin_lineage_flags_foreign_class():
    bin_contigs = {'bin1': ['c1', 'c2', 'c3']}
    lineage = {
        'c1': "d_Bacteria;p_Pseudomonadota;c_Gammaproteobacteria",
        'c2': "d_Bacteria;p_Pseudomonadota;c_Gammaproteobacteria",
        'c3': "d_Bacteria;p_Pseudomonadota;c_Alphaproteobacteria",  # foreign
    }
    out = summarize_bin_lineage(bin_contigs, lineage, 'class')
    rec = out['bin1']
    assert rec['chimera'] is True and rec['n_taxa'] == 2
    assert rec['majority_taxon'] == 'Gammaproteobacteria'   # key is the name string
    assert [c for c, _ in rec['foreign_contigs']] == ['c3']


def test_summarize_bin_lineage_unresolved_above_rank():
    # a contig classified only to phylum is unresolved at class -> no split
    bin_contigs = {'b': ['c1', 'c2']}
    lineage = {'c1': "d_Bacteria;p_Pseudomonadota;c_Gammaproteobacteria",
               'c2': "d_Bacteria;p_Pseudomonadota"}
    out = summarize_bin_lineage(bin_contigs, lineage, 'class')
    assert out['b']['n_taxa'] == 1 and out['b']['chimera'] is False


def test_classify_contigs_nodes_optional_for_mmseqs():
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers()
    classify_contigs.register_parser(sub)
    # MMseqs2 path without nodes.dmp parses fine (lineage mode)
    args = parser.parse_args(['classify-contigs', '--bins-dir', 'bins',
                              '--mmseqs-output', 'c.tsv', '--rank', 'class'])
    assert args.nodes_dmp is None and args.rank == 'class'


def test_mmseqs_run_requires_db_and_output():
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers()
    run_mmseqs_contigs.register_parser(sub)
    args = parser.parse_args(['run-mmseqs-contigs', '--bins-dir', 'bins',
                              '--db', '/db/gtdb', '-o', 'out.tsv'])
    assert args.threads == 1 and args.keep_tmp is False
