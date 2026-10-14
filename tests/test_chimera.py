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
    has_taxonomy_warning,
    assess_chimera_risk,
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


def _risk(**kw):
    base = dict(gc_delta=0.0, gc_cv=0.0, n_contigs=1,
                contamination=None, windowed_gc_delta=None)
    base.update(kw)
    return assess_chimera_risk(**base)


def test_assess_chimera_risk_clean():
    risk, reasons = _risk()
    assert risk == "Clean" and reasons == []


def test_assess_chimera_risk_gc_high():
    # inter-contig GC >10% is now +1 -> Low on its own (never Medium)
    risk, reasons = _risk(gc_delta=0.12, n_contigs=3)
    assert risk == "Low"
    assert any("GC heterogeneity" in r for r in reasons)


def test_assess_chimera_risk_gunc_fail_is_high():
    # GUNC fail is +3 -> High on its own (authoritative gene-level signal)
    risk, reasons = _risk(gunc_fail=True, gunc_css=0.82)
    assert risk == "High"
    assert any("GUNC fail" in r for r in reasons)


def test_assess_chimera_risk_gunc_css_elevated_is_low():
    risk, reasons = _risk(gunc_fail=False, gunc_css=0.5)   # +1
    assert risk == "Low"
    assert any("clade separation" in r.lower() for r in reasons)


def test_contamination_only_over_10_and_hq_tolerated():
    assert _risk(contamination=12.0)[0] == "Low"     # +1 -> Low
    assert _risk(contamination=7.0)[0] == "Clean"    # <=10% tolerated -> 0


def test_gtdbtk_warning_half_point():
    risk, _ = _risk(taxonomy_warning=True)            # +0.5 -> Low
    assert risk == "Low"


def test_assess_chimera_risk_high_combo():
    # GC (+1) + GUNC fail (+3) -> High
    risk, _ = _risk(gc_delta=0.2, n_contigs=2, gunc_fail=True, gunc_css=0.9)
    assert risk == "High"


def test_windowed_gc_alone_is_low_not_medium():
    # a circular contig with only high windowed GC (+0.5) stays Low
    risk, reasons = _risk(windowed_gc_delta=0.169, n_contigs=1)
    assert risk == "Low"
    assert any("Windowed GC" in r for r in reasons)


def test_contig_foreign_frac_tiers():
    # length-weighted per-contig foreign fraction (MMseqs via classify-contigs)
    assert _risk(contig_foreign_frac=0.12)[0] == "High"    # >=10% +3
    assert _risk(contig_foreign_frac=0.07)[0] == "Medium"  # >=5%  +2
    assert _risk(contig_foreign_frac=0.02)[0] == "Low"     # >0%   +0.5
    risk, reasons = _risk(contig_foreign_frac=0.0)
    assert risk == "Clean" and not any("foreign" in r for r in reasons)


def test_gc_outlier_note_in_reason():
    outlier = {"high": ("contigA", 0.65, 500_000),
               "low": ("contigB", 0.35, 1_200)}
    risk, reasons = _risk(gc_delta=0.12, n_contigs=2, gc_outlier=outlier)
    joined = " ".join(reasons)
    assert "contigA" in joined and "contigB" in joined
    assert "1,200bp" in joined          # short outlier length shown for review


def test_check_chimeras_loads_contig_foreign_frac(tmp_path):
    import argparse
    from echidna.commands import check_chimeras
    f = tmp_path / "ct.bins.tsv"
    f.write_text("bin\tchimera_flag\tforeign_frac\tn_foreign_contigs\n"
                 "binA\tTrue\t0.1200\t2\nbinB\tFalse\t0.0000\t0\n")
    out = check_chimeras._load_contig_foreign_frac(f)
    assert out['binA'] == pytest.approx(0.12) and out['binB'] == 0.0

    # parser gained --contig-taxonomy / --name-strip and dropped the GUNC per-contig flags
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers()
    check_chimeras.register_parser(sub)
    args = parser.parse_args(['check-chimeras', '--bins-dir', 'bins',
                              '--contig-taxonomy', 'ct.bins.tsv',
                              '--name-strip', '.contigs,_genomic'])
    assert str(args.contig_taxonomy) == 'ct.bins.tsv'
    assert args.name_strip == '.contigs,_genomic'
    with pytest.raises(SystemExit):
        parser.parse_args(['check-chimeras', '--bins-dir', 'b', '--gunc-tax-rank', 'order'])
