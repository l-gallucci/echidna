"""Tests for echidna.utils.names (canonical cross-tool genome-name matching)."""

from echidna.utils.names import canonical_name, match_report


def test_canonical_name_strips_extension_and_lowercases():
    assert canonical_name('binA.fa') == 'bina'
    assert canonical_name('BinA.fasta') == 'bina'
    assert canonical_name('bin_A.fna.gz') == 'bin_a'
    assert canonical_name('bin.1') == 'bin.1'          # not a known ext


def test_canonical_name_strips_custom_suffix():
    assert canonical_name('binA.contigs', ['.contigs']) == 'bina'
    assert canonical_name('genome_genomic.fa', ['_genomic']) == 'genome'
    # suffix only after the extension is stripped
    assert canonical_name('g.fasta', ['_genomic']) == 'g'


def test_match_report_matched_and_unmatched():
    bins = ['binA.fa', 'binB.fa', 'binC.fa']
    source = ['BINA', 'binB.fasta']          # different case / extension
    matched, unmatched = match_report(bins, source)
    assert set(matched) == {'binA.fa', 'binB.fa'}
    assert unmatched == ['binC.fa']


def test_match_report_with_suffix():
    bins = ['binA.fa']
    source = ['binA.contigs']
    matched, unmatched = match_report(bins, source, ['.contigs'])
    assert matched == ['binA.fa'] and unmatched == []
