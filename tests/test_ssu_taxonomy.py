"""Tests for echidna.utils.ssu_taxonomy and the check-ssu-taxonomy command.
Pure-function tests; blastn/barrnap are not executed."""

import argparse
from pathlib import Path

import pytest

from echidna.utils.ssu_taxonomy import (
    extract_ssu_records, longest_ssu, parse_blast_best,
    ssu_taxon_at_rank, silva_lineage_parts, normalize_taxon, taxa_agree,
)
from echidna.commands import check_ssu_taxonomy, run_mimag_rna


def test_extract_ssu_records_filters_16s(tmp_path):
    f = tmp_path / "binA.rrna.fna"
    f.write_text(
        ">16S_rRNA::c1:1-1500(+)\nAAAACCCC\n"
        ">23S_rRNA::c1:2000-4800(+)\nGGGGTTTT\n"
        ">18S_rRNA::c2:1-1800(+)\nACGTACGT\n"
        ">5S_rRNA::c1:5000-5100(+)\nTTTT\n")
    recs = extract_ssu_records(f)
    names = [n for n, _ in recs]
    assert any(n.startswith('16S_rRNA') for n in names)
    assert not any(n.startswith(('18S_rRNA', '23S_rRNA', '5S_rRNA')) for n in names)


def test_longest_ssu():
    recs = [('16S_rRNA::a', 'AAA'), ('16S_rRNA::b', 'AAAAAAAA')]
    assert longest_ssu(recs)[0] == '16S_rRNA::b'
    assert longest_ssu([]) is None


def test_collect_ssu_records_keeps_multiple_16s_only(tmp_path):
    d = tmp_path / "rrna"
    d.mkdir()
    (d / "binA.rrna.fna").write_text(
        ">16S_rRNA::c1:1-1500(+)\n" + "A" * 1500 + "\n"
        ">16S_rRNA::c2:1-1400(+)\n" + "C" * 1400 + "\n"
        ">18S_rRNA::c3:1-1800(+)\n" + "G" * 1800 + "\n"
        ">23S_rRNA::c4:1-2900(+)\n" + "T" * 2900 + "\n")
    out = check_ssu_taxonomy.collect_ssu_records(d)
    assert [r['type'] for r in out['binA']] == ['16S', '16S']
    assert [r['len'] for r in out['binA']] == [1500, 1400]


def test_parse_blast_best_keeps_highest_bitscore_and_stitle(tmp_path):
    f = tmp_path / "b.tsv"
    # qseqid sseqid pident length slen bitscore stitle(may contain spaces/;)
    f.write_text(
        "binA\tacc1\t88.0\t1400\t1500\t900.0\tacc1 d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria\n"
        "binA\tacc2\t99.0\t1480\t1500\t1400.0\tacc2 d__Bacteria;p__Pseudomonadota;c__Alphaproteobacteria\n")
    best = parse_blast_best(f)
    assert best['binA']['bitscore'] == 1400.0        # higher bitscore wins
    assert 'Alphaproteobacteria' in best['binA']['stitle']
    assert best['binA']['pident'] == 99.0


def test_ssu_taxon_at_rank_gtdb():
    stitle = "acc1 d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales"
    assert ssu_taxon_at_rank(stitle, 'class', 'gtdb') == 'Gammaproteobacteria'
    assert ssu_taxon_at_rank(stitle, 'order', 'gtdb') == 'Enterobacterales'
    assert ssu_taxon_at_rank(stitle, 'genus', 'gtdb') is None


def test_ssu_taxon_at_rank_silva_positional():
    # SILVA: "<accession> Domain;Phylum;Class;Order;Family;Genus;Species"
    stitle = "AY123.1.1500 Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;Ent;Esch;coli"
    assert ssu_taxon_at_rank(stitle, 'class', 'silva') == 'Gammaproteobacteria'
    assert ssu_taxon_at_rank(stitle, 'phylum', 'silva') == 'Proteobacteria'
    assert ssu_taxon_at_rank(stitle, 'order', 'silva') == 'Enterobacterales'


def test_ssu_taxon_at_rank_silva_handles_species_spaces_and_no_accession():
    stitle = ("Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;"
              "Enterobacteriaceae;Escherichia-Shigella;Escherichia coli")
    assert ssu_taxon_at_rank(stitle, 'domain', 'silva') == 'Bacteria'
    assert ssu_taxon_at_rank(stitle, 'genus', 'silva') == 'Escherichia-Shigella'
    assert ssu_taxon_at_rank(stitle, 'species', 'silva') == 'Escherichia coli'


def test_silva_lineage_parts_drops_accession_without_breaking_taxa():
    stitle = ("AB123456.1.1490 Archaea;Halobacterota;Methanosarcinia;"
              "Methanosarcinales;Methanosarcinaceae;Methanosarcina;"
              "Methanosarcina mazei")
    parts = silva_lineage_parts(stitle)
    assert parts[0] == 'Archaea'
    assert parts[5] == 'Methanosarcina'
    assert parts[6] == 'Methanosarcina mazei'


def test_normalize_taxon_strips_gtdb_suffix():
    assert normalize_taxon('Firmicutes_A') == 'firmicutes'
    assert normalize_taxon('Gammaproteobacteria') == 'gammaproteobacteria'
    assert normalize_taxon('Unclassified') == ''
    assert normalize_taxon(None) == ''


def test_taxa_agree():
    assert taxa_agree('Gammaproteobacteria', 'Gammaproteobacteria') is True
    assert taxa_agree('Alphaproteobacteria', 'Gammaproteobacteria') is False
    assert taxa_agree('Firmicutes_A', 'Firmicutes') is True       # suffix-normalised
    assert taxa_agree(None, 'Gammaproteobacteria') is None         # unknown -> no judgement
    assert taxa_agree('Gammaproteobacteria', 'Unclassified') is None


def test_check_ssu_taxonomy_parser_defaults():
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers()
    check_ssu_taxonomy.register_parser(sub)
    args = parser.parse_args(['check-ssu-taxonomy', '--ssu-fasta-dir', 'rrna',
                              '--ssu-db', '/db/ssu', '--gtdbtk-dir', 'tax'])
    assert args.rank == 'order'
    assert args.ssu_taxonomy == 'gtdb'
    assert args.min_length == 1000


def test_run_mimag_rna_outseq_flag_default_off():
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers()
    run_mimag_rna.register_parser(sub)
    args = parser.parse_args(['run-mimag-rna', '--genome-dir', 'g', '-o', 'out'])
    assert args.outseq is False
    args2 = parser.parse_args(['run-mimag-rna', '--genome-dir', 'g', '-o', 'out', '--outseq'])
    assert args2.outseq is True
