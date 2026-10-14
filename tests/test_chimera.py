"""Tests for echidna.utils.chimera."""

import math
import pytest

from echidna.utils.chimera import (
    compute_gc,
    windowed_gc,
    gc_statistics,
    extract_taxonomy_level,
    parse_gtdbtk_summary,
    parse_gunc_output,
    parse_gunc_contig_assignments,
    gc_outlier_contigs,
    strong_gc_outlier,
    has_taxonomy_warning,
    assess_chimera_verdict,
)


def test_compute_gc():
    assert compute_gc("GGCC") == 1.0
    assert compute_gc("ATAT") == 0.0
    assert compute_gc("ACGT") == 0.5
    assert math.isnan(compute_gc("NNNN"))


def test_windowed_gc_short_sequence():
    # shorter than window -> single whole-sequence value
    out = windowed_gc("ACGTACGT", window=5000, step=2500)
    assert len(out) == 1 and out[0][0] == 0


def test_windowed_gc_sliding():
    seq = "GC" * 5000          # 10 kbp, all GC
    out = windowed_gc(seq, window=1000, step=1000)
    assert len(out) >= 9
    assert all(abs(gc - 1.0) < 1e-9 for _, gc in out)


def test_gc_statistics():
    s = gc_statistics([0.4, 0.6])
    assert s["mean"] == pytest.approx(0.5)
    assert s["delta"] == pytest.approx(0.2)
    assert s["n"] == 2
    empty = gc_statistics([float("nan")])
    assert empty["n"] == 0


@pytest.mark.parametrize("level,expected", [
    ("genus", "Lactobacillus"),
    ("phylum", "Firmicutes"),
    ("species", "Unclassified"),   # empty s__
])
def test_extract_taxonomy_level(level, expected):
    cls = "d__Bacteria;p__Firmicutes;c__Bacilli;o__Lacto;f__Lacto;g__Lactobacillus;s__"
    assert extract_taxonomy_level(cls, level) == expected


def test_parse_gtdbtk_summary(tmp_path):
    f = tmp_path / "gtdbtk.bac120.summary.tsv"
    f.write_text("user_genome\tclassification\tmsa_percent\twarnings\n"
                 "binA\td__Bacteria;p__Firmicutes\t95.0\t\n")
    out = parse_gtdbtk_summary(f)
    assert out["binA"]["classification"].startswith("d__Bacteria")
    assert out["binA"]["msa_percent"] == 95.0


def test_parse_gunc_output(tmp_path):
    f = tmp_path / "GUNC.progenomes_2.1.maxCSS_level.tsv"
    f.write_text(
        "genome\tclade_separation_score\tcontamination_portion\t"
        "n_effective_surplus_clades\tpass.GUNC\n"
        "binA\t0.82\t0.4\t2.1\tFalse\n"
        "binB\t0.10\t0.02\t0.3\tTrue\n")
    out = parse_gunc_output(tmp_path)
    assert out["binA"]["pass_gunc"] is False
    assert out["binA"]["css"] == 0.82
    assert out["binB"]["pass_gunc"] is True


def test_parse_gunc_contig_assignments(tmp_path):
    # binA: two contigs resolve to different phyla -> chimera (n_clades=2)
    (tmp_path / "binA.contig_assignments.tsv").write_text(
        "contig\ttax_level\tassignment\tcount_of_genes_assigned\n"
        "c1\tphylum\tProteobacteria\t40\n"
        "c1\tphylum\tFirmicutes\t2\n"          # stray, below min_genes -> ignored
        "c2\tphylum\tBacteroidota\t25\n"
        "c3\tphylum\tProteobacteria\t10\n")
    # binB: all contigs one phylum -> clean (n_clades=1)
    (tmp_path / "binB.contig_assignments.tsv").write_text(
        "contig\ttax_level\tassignment\tcount_of_genes_assigned\n"
        "c1\tphylum\tFirmicutes\t50\n"
        "c2\tphylum\tFirmicutes\t30\n")

    out = parse_gunc_contig_assignments(tmp_path, rank="phylum", min_genes=3)
    assert out["binA"]["n_clades"] == 2
    assert set(out["binA"]["clades"]) == {"Proteobacteria", "Bacteroidota"}
    assert out["binB"]["n_clades"] == 1


def test_parse_gunc_contig_assignments_min_genes_filters_strays(tmp_path):
    (tmp_path / "binC.contig_assignments.tsv").write_text(
        "contig\ttax_level\tassignment\tcount_of_genes_assigned\n"
        "c1\tphylum\tProteobacteria\t40\n"
        "c2\tphylum\tFirmicutes\t1\n")        # single stray gene, filtered
    out = parse_gunc_contig_assignments(tmp_path, rank="phylum", min_genes=3)
    assert out["binC"]["n_clades"] == 1       # stray did not invent a clade


def test_gc_outlier_contigs():
    gc = {"a": 0.60, "b": 0.35, "c": 0.50}
    length = {"a": 100000, "b": 3000, "c": 50000}
    out = gc_outlier_contigs(gc, length)
    assert out["high"] == ("a", 0.60, 100000)
    assert out["low"] == ("b", 0.35, 3000)
    assert gc_outlier_contigs({"only": 0.5}, {"only": 1000}) is None


def test_has_taxonomy_warning():
    assert has_taxonomy_warning({"warnings": "something"}) is True
    assert has_taxonomy_warning({"msa_percent": 5.0}) is True
    assert has_taxonomy_warning({"msa_percent": 95.0, "warnings": ""}) is False


# ---- verdict from two signals: GUNC + per-contig taxonomy (MMseqs) ----------

def test_verdict_both_signals_chimeric():
    # GUNC fail AND taxonomy chimeric -> Chimeric
    v, reasons = assess_chimera_verdict(gunc_fail=True, tax_chimeric=True,
                                        gunc_css=0.82, contig_foreign_frac=0.2)
    assert v == "Chimeric"
    assert any("GUNC fail" in r for r in reasons)
    assert any("Per-contig taxonomy" in r for r in reasons)


def test_verdict_one_signal_is_check():
    assert assess_chimera_verdict(gunc_fail=True, tax_chimeric=False)[0] == "Check"
    assert assess_chimera_verdict(gunc_fail=False, tax_chimeric=True)[0] == "Check"


def test_verdict_neither_is_clean():
    assert assess_chimera_verdict(gunc_fail=False, tax_chimeric=False)[0] == "Clean"


def test_verdict_missing_signal_uses_the_available_one():
    # only one signal present and firing -> Check (cannot confirm both)
    assert assess_chimera_verdict(gunc_fail=True, tax_chimeric=None)[0] == "Check"
    assert assess_chimera_verdict(gunc_fail=None, tax_chimeric=True)[0] == "Check"
    # only one present and clean -> Clean
    assert assess_chimera_verdict(gunc_fail=False, tax_chimeric=None)[0] == "Clean"
    # neither available -> Unknown
    assert assess_chimera_verdict(gunc_fail=None, tax_chimeric=None)[0] == "Unknown"


def test_gc_not_in_verdict():
    # GC / contamination / GTDB do not even enter assess_chimera_verdict
    import inspect
    params = set(inspect.signature(assess_chimera_verdict).parameters)
    assert params == {'gunc_fail', 'tax_chimeric', 'gunc_css', 'contig_foreign_frac'}


# ---- strong GC outlier (reported column only, length-filtered) --------------

def test_strong_gc_outlier_ignores_short_contigs():
    gc = {'long_norm': 0.50, 'long_dev': 0.68, 'tiny': 0.95}
    length = {'long_norm': 100_000, 'long_dev': 80_000, 'tiny': 200}
    out = strong_gc_outlier(gc, length, min_len=5000, max_dev=0.10)
    assert out is not None
    assert out['contig'] == 'long_dev'          # tiny 200bp ignored despite 95% GC
    assert out['gc'] == 0.68


def test_strong_gc_outlier_none_when_within_mean():
    gc = {'a': 0.50, 'b': 0.52, 'c': 0.49}
    length = {'a': 100_000, 'b': 90_000, 'c': 80_000}
    assert strong_gc_outlier(gc, length, min_len=5000, max_dev=0.10) is None


def test_check_chimeras_verdict_wiring(tmp_path):
    import argparse
    from echidna.commands import check_chimeras
    f = tmp_path / "ct.bins.tsv"
    f.write_text("bin\tchimera_flag\tforeign_frac\tn_foreign_contigs\n"
                 "binA\tTrue\t0.1200\t2\nbinB\tFalse\t0.0000\t0\n")
    out = check_chimeras._load_contig_foreign_frac(f)
    assert out['binA'] == pytest.approx(0.12) and out['binB'] == 0.0

    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers()
    check_chimeras.register_parser(sub)
    args = parser.parse_args(['check-chimeras', '--bins-dir', 'bins',
                              '--gunc-dir', 'gunc', '--contig-taxonomy', 'ct.bins.tsv',
                              '--name-strip', '.contigs'])
    assert str(args.contig_taxonomy) == 'ct.bins.tsv'
    assert args.min_foreign_frac == 0.05 and args.name_strip == '.contigs'
    # the removed GUNC per-contig flag no longer parses
    with pytest.raises(SystemExit):
        parser.parse_args(['check-chimeras', '--bins-dir', 'b', '--gunc-tax-rank', 'order'])
