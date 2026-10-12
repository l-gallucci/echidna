"""Tests for echidna.utils.final_report and the final-report command.
Pure-function tests."""

import argparse
from pathlib import Path

import pytest

from echidna.utils.final_report import (
    norm_name, build_report, parse_chimera_report,
    parse_contig_taxonomy_bins, parse_ssu_report,
)
from echidna.commands import final_report


def test_norm_name_strips_extension():
    assert norm_name('binA.fa') == 'binA'
    assert norm_name('binA.fasta') == 'binA'
    assert norm_name('binA.fna.gz') == 'binA'
    assert norm_name('bin.1') == 'bin.1'        # not a known ext -> untouched


def test_parse_chimera_report(tmp_path):
    f = tmp_path / "chim.tsv"
    f.write_text("bin\tchimera_risk\nbinA.fa\tHigh\nbinB\tClean\n")
    out = parse_chimera_report(f)
    assert out == {'binA': 'High', 'binB': 'Clean'}


def test_parse_contig_taxonomy_bins(tmp_path):
    f = tmp_path / "ct.bins.tsv"
    f.write_text("bin\tn_taxa\tchimera_flag\tn_foreign_contigs\n"
                 "binA\t2\tTrue\t1\nbinB\t1\tFalse\t0\n")
    out = parse_contig_taxonomy_bins(f)
    assert out['binA']['contig_chimera'] is True and out['binA']['n_foreign'] == '1'
    assert out['binB']['contig_chimera'] is False


def test_parse_ssu_report(tmp_path):
    f = tmp_path / "ssu.tsv"
    f.write_text("genome\tagree\tssu_taxon\nbinA\tFalse\tAlphaproteobacteria\n"
                 "binB\tTrue\tGammaproteobacteria\n")
    out = parse_ssu_report(f)
    assert out['binA']['agree'] == 'False'
    assert out['binB']['ssu_taxon'] == 'Gammaproteobacteria'


def test_parse_ssu_report_aggregates_multiple_rows(tmp_path):
    f = tmp_path / "ssu.tsv"
    f.write_text("genome\tagree\tssu_taxon\n"
                 "binA\tTrue\tGammaproteobacteria\n"
                 "binA\tFalse\tAlphaproteobacteria\n"
                 "binB\tNA\tNA\n"
                 "binB\tTrue\tBacilli\n")
    out = parse_ssu_report(f)
    assert out['binA']['agree'] == 'False'
    assert out['binA']['ssu_taxon'] == 'Gammaproteobacteria;Alphaproteobacteria'
    assert out['binB']['agree'] == 'True'
    assert out['binB']['ssu_taxon'] == 'Bacilli'


def test_build_report_join_and_concerns():
    universe = {'binA', 'binB'}
    checkm2 = {'binA': {'completeness': 95.0, 'contamination': 12.0},   # contam>10
               'binB': {'completeness': 80.0, 'contamination': 2.0}}
    taxonomy = {'binA': {'classification': 'd__Bacteria;p__Pseudomonadota',
                         'rank_taxon': 'Pseudomonadota'}}
    mimag = {'binB': {'n_16S': 1, 'rrna_complete': True, 'trna_complete': True,
                      'mimag_rna_ok': True}}
    chimera = {'binA': 'High', 'binB': 'Clean'}
    contig = {'binA': {'contig_chimera': True, 'n_foreign': '1'}}
    ssu = {'binB': {'agree': 'False', 'ssu_taxon': 'X'}}

    rows = build_report(universe, checkm2=checkm2, taxonomy=taxonomy, mimag=mimag,
                        chimera=chimera, contig=contig, ssu=ssu)
    by = {r['genome']: r for r in rows}

    a = by['binA']
    # contamination 12% -> Fail tier (and a concern)
    assert a['tier'] == 'Fail' and a['chimera_risk'] == 'High'
    # concerns gather the red flags for binA
    assert 'contam>10' in a['concerns'] and 'chimera_High' in a['concerns']
    assert 'contig_chimera' in a['concerns']

    b = by['binB']
    assert b['mimag_rna_ok'] is True
    assert 'ssu_disagree' in b['concerns']        # agree == 'False'
    assert b['tier'] == 'MQ'


def test_build_report_missing_sources_are_na():
    rows = build_report({'binX'}, checkm2={'binX': {'completeness': 91.0,
                                                    'contamination': 1.0}})
    r = rows[0]
    assert r['taxonomy'] == 'NA' and r['chimera_risk'] == 'NA'
    assert r['ssu_agree'] == 'NA' and r['concerns'] == ''   # no red flags


def test_final_report_parser_defaults():
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers()
    final_report.register_parser(sub)
    args = parser.parse_args(['final-report', '--chimera-report', 'c.tsv'])
    assert args.tax_level == 'phylum'
    assert str(args.output) == 'final_report.tsv'
