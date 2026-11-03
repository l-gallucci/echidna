#!/usr/bin/env python3
"""Riproduzioni dei difetti osservati, NON test di correttezza della pipeline.

Esecuzione: python reproduce_findings.py /percorso/echidna
Usa solo stdlib, fixture sintetiche e mock dichiarati; nessun tool bioinformatico.
Un PASS conferma il comportamento difettoso descritto nell'audit alla revisione
esaminata. Dopo una correzione le relative asserzioni devono essere aggiornate.
Non scrive nel repository esaminato.
"""
import argparse
import ast
import contextlib
import csv
import io
import json
import re
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path
from unittest.mock import patch

sys.dont_write_bytecode = True
repo = Path(sys.argv[1]).resolve()
sys.path.insert(0, str(repo))
from echidna.commands import process_large_contigs as plc
from echidna.commands import filter_checkm2 as org
from echidna.commands import filter_misassemblies as filt
from echidna.commands import check_read_chimeras as readchim
from echidna.commands import classify_clipping as clipping
from echidna.commands import misassembly_wf, mag_quality_wf, run_mimag_rna
from echidna.commands import check_chimeras, classify_contigs, remove_hq_contigs
from echidna.commands.setup_workflow import set_yaml_scalar
from echidna.utils import bam_stats, zero_coverage, contig_taxonomy, final_report
from echidna.utils import misassembly
from echidna.utils.read_taxonomy import compute_taxonomy_stats

results = []


def probe(name, function):
    try:
        with tempfile.TemporaryDirectory(prefix="echidna_audit_") as temp:
            with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                evidence = function(Path(temp))
        results.append({"id": name, "confirmed": True, "evidence": evidence})
    except Exception as exc:
        results.append({"id": name, "confirmed": False, "error": repr(exc)})


def parser_for(module, argv):
    parser = argparse.ArgumentParser()
    module.register_parser(parser.add_subparsers(dest="command"))
    return parser.parse_args(argv)


def hash_mismatch(t):
    assembly = t / "assembly.fa"
    assembly.write_text(">contig_7\n" + "ACGT" * 100 + "\n")
    misasm = t / "misasm"
    misasm.mkdir()
    (misasm / "test-clipping.txt").write_text(
        "contig_name\tsample_name\tposition\tclipping_coverage\ncontig_7\ts\t200\t100\n")
    qc = t / "qc"
    qc.mkdir()
    hashed = plc.safe_filename("contig_7")
    (qc / "quality_report.tsv").write_text(
        f"Name\tCompleteness\tContamination\n{hashed}\t95\t1\n")
    args = parser_for(plc, ["process-large-contigs", "-a", str(assembly), "-m", str(misasm),
                           "-o", str(t / "out"), "--max-length", "100", "--skip-checkm2",
                           "--checkm2-results", str(qc)])
    plc.run(args)
    assert (t / "out/HQ_extracted" / f"{hashed}.fa").exists()
    assert not (t / "out/assembly_for_filtering.fa").read_text()
    return {"contig": "contig_7", "checkm2_name": hashed,
            "has_clipping": True, "observed": "EXTRACT_HQ; absent from filtering assembly"}


def no_large(t):
    a = t / "a.fa"
    a.write_text(">c\nACGT\n")
    m = t / "misasm"
    m.mkdir()
    args = parser_for(plc, ["process-large-contigs", "-a", str(a), "-m", str(m), "-o", str(t / "o")])
    plc.run(args)
    assert not (t / "o/assembly_for_filtering.fa").exists()
    return "Exit normale senza assembly_for_filtering.fa necessario al passo successivo"


def missing_quality(t):
    a = t / "a.fa"
    a.write_text(">c1\n" + "A" * 200 + "\n>c2\n" + "C" * 200 + "\n")
    m = t / "m"
    m.mkdir()
    q = t / "q"
    q.mkdir()
    (q / "quality_report.tsv").write_text(
        f"Name\tCompleteness\tContamination\n{plc.safe_filename('c1')}\t10\t0\n")
    args = parser_for(plc, ["process-large-contigs", "-a", str(a), "-m", str(m), "-o", str(t / "o"),
                           "--max-length", "100", "--skip-checkm2", "--checkm2-results", str(q)])
    plc.run(args)
    assert ">c2" not in (t / "o/assembly_for_filtering.fa").read_text()
    assert not list((t / "o/HQ_extracted").glob("*.fa"))
    return "c2 non presente nel report CheckM2 scompare sia dall'assembly aggiornato sia da HQ_extracted"


def wrapper_windows(t):
    args = parser_for(misassembly_wf, ["misassembly-wf", "-a", "a.fa", "-b", "m.bam",
                                     "-o", "out", "--kraken2-output", "r.tsv"])
    plan = misassembly_wf.build_plan(args)
    rc = next(s for s in plan if s.key == "readchim")
    assert "--assembly" not in rc.cmd
    p = readchim.ContigReadProfile("c", length=0)
    for i in range(20):
        p.add_read(i * 1000, 2 if i < 10 else 3)
    rows, windows = readchim.analyze_profiles({"c": p}, {}, 10000, 5000, 10, 5, True)
    assert rows and not windows
    assert "--taxonomy-windows" in next(s for s in plan if s.key == "classify").cmd
    return "readchim omette --assembly, non produce finestre, classify richiede quel file"


def window_schema(t):
    f = t / "windows.tsv"
    readchim.write_tsv([{"contig": "c", "window_start": 100, "window_end": 200,
                         "is_taxon_shift": True}], readchim.WINDOW_FIELDS, f)
    assert clipping.load_taxonomy_shifts(f) == {}
    return "is_taxon_shift=True scritto dal produttore; il lettore cerca is_shift e perde il segnale"


def long_read_stats(t):
    sam = "r\t0\tc\t1000\t60\t9000M100S\t*\t0\t0\t" + "ACGT" * 2275 + "\t*\tSA:Z:other,1,+,100M,60,0;\n"
    class FakeProc:
        def __init__(self):
            self.stdout = io.StringIO(sam)
            self.stderr = io.StringIO("")
            self.returncode = 0
        def wait(self):
            return 0
    with patch.object(bam_stats.subprocess, "Popen", return_value=FakeProc()):
        stats = bam_stats.collect_position_stats(Path("mock.bam"), {"c": [(10000, 1)]},
                                                  500, 20, {"c": 1}, {"c": 20000}, 0, 0)
    assert stats[("c", 10000)]["n_clipped_reads"] == 0
    return "Read 9000M100S in POS=1000 termina e clippa in 10000, ma non viene conteggiata a 10000"


def unpaired(t):
    acc = bam_stats._PosAccumulator(10000, "c")
    for i in range(10):
        acc.add_read(0 if i % 2 else 16, 0, "ACGT" * 25, "100M", 10000, 0, 0)
    d = acc.to_dict("c", 20000, 10)
    assert d["discordant_fraction"] == 1
    label, _, _ = clipping.classify_position(d)
    assert label == "chimera_candidate"
    return {"unpaired_reads": 10, "discordant_fraction": 1, "label": label}


def blast_failure(t):
    exc = subprocess.CalledProcessError(2, ["blastn"], stderr="synthetic DB error")
    with patch.object(zero_coverage.subprocess, "run", side_effect=exc):
        hits = zero_coverage.blast_regions_vs_reads(Path("q.fa"), Path("db"))
    verdict = zero_coverage.assign_verdict(hits.get("region", 0), float("nan"), False)
    assert hits == {} and verdict == "unsupported"
    return {"mock_blast_exit": 2, "verdict": verdict}


def missing_rna(t):
    records = [{"name": "x", "tier": "HQ"}]
    org.apply_mimag_rna_gate(records, {})
    assert records[0]["tier"] == "HQ"
    return records


def stale_tiers(t):
    bins = t / "bins"
    bins.mkdir()
    (bins / "x.fa").write_text(">c\nACGT\n")
    out = t / "out"
    for tier in ("HQ", "MQ"):
        org.link_or_copy_genomes([{"name": "x", "tier": tier}], bins, out,
                                 ["HQ", "MQ"], ["fa"], {}, {}, mode="copy")
    assert (out / "Selected/HQ/x.fa").exists() and (out / "Selected/MQ/x.fa").exists()
    return "x.fa rimane in HQ dopo una nuova selezione MQ nella stessa directory"


def unresolved(t):
    bins = t / "bins"
    bins.mkdir()
    (bins / "x.fa").write_text(">c\nACGT\n")
    summary = contig_taxonomy.summarize_bin_lineage({"x": ["c"]}, {}, "order", {"c": 4})
    _, rows = classify_contigs.build_rows(summary, {})
    table = t / "taxonomy.bins.tsv"
    classify_contigs.write_tsv(rows, classify_contigs.BIN_FIELDS, table)
    output = t / "chimera.tsv"
    check_chimeras.run_chimera_check(bins, output, ["fa"], None, None, None, table,
                                    .05, 5000, .1, "phylum", [], False)
    row = list(csv.DictReader(output.open(), delimiter="\t"))[0]
    assert row["chimera_verdict"] == "Clean" and rows[0]["n_contigs_resolved"] == 0
    return {"n_contigs_resolved": 0, "foreign_frac": 0, "verdict": "Clean"}


def missing_reads(t):
    stats = compute_taxonomy_stats(Counter({-1: 100}), {}, True)
    assert stats["n_classified"] == 100 and stats["classified_fraction"] == 1
    return stats


def remove_by_bin_filename(t):
    hq = t / "hq"
    hq.mkdir()
    (hq / "bin.42.fa").write_text(">contig_7\nAAAA\n>contig_9\nCCCC\n")
    a = t / "a.fa"
    a.write_text(">contig_7\nAAAA\n>contig_9\nCCCC\n>contig_42\nGGGG\n")
    ids = remove_hq_contigs.load_hq_ids_from_dir(hq)
    out = t / "out.fa"
    remove_hq_contigs.process_assembly(a, ids, out)
    content = out.read_text()
    assert ">contig_42" not in content and ">contig_7" in content and ">contig_9" in content
    return "HQ bin.42.fa contiene contig_7 e contig_9; viene rimosso contig_42 estraneo al MAG"


def final_tier(t):
    cm = {"x": {"completeness": 95., "contamination": 1.}}
    records = [{"name": "x", "tier": "HQ"}]
    org.apply_mimag_rna_gate(records, {"x": {"mimag_rna_ok": False}})
    rows = final_report.build_report({"x"}, checkm2=cm, mimag={"x": {"mimag_rna_ok": False}})
    assert records[0]["tier"] == "MQ" and rows[0]["tier"] == "HQ"
    return {"organize_tier": "MQ", "final_report_tier": "HQ"}


def yaml_regex(t):
    try:
        set_yaml_scalar('sample_regex: "x"\n', "sample_regex", r"^(S\d+)_", quote=True)
    except re.error as exc:
        return str(exc)
    raise AssertionError("Expected re.error")


def extension_plan(t):
    args = parser_for(mag_quality_wf, ["mag-quality-wf", "--bins-dir", "bins", "-o", "out", "--genome-ext", "fna"])
    mi = next(s for s in mag_quality_wf.build_plan(args) if s.key == "mimag")
    mi_args = parser_for(run_mimag_rna, mi.cmd[1:])
    assert mi_args.extension == "fa"
    return "--genome-ext fna non inoltrato a MIMAG: il comando successivo cerca *.fa"


def verdict_not_consumed(t):
    args = parser_for(mag_quality_wf, ["mag-quality-wf", "--bins-dir", "bins", "-o", "out"])
    final = next(s for s in mag_quality_wf.build_plan(args) if s.key == "organize2")
    assert not any("chimera" in token for token in final.cmd)
    return "organize2 non riceve il report chimera; selezione finale basata su CheckM2 e RNA"


def verdict_schema(t):
    f = t / "chim.tsv"
    f.write_text("bin\tchimera_verdict\nx\tChimeric\n")
    assert filt.load_chimera_flagged(f) == set()
    return "filter-misassemblies accetta il report di check-chimeras ma legge solo chimera_risk"


def zero_clips(t):
    a = t / "a.fa"
    a.write_text(">c\n" + "ACGT" * 1000 + "\n")
    m = t / "m"
    m.mkdir()
    (m / "x-clipping.txt").write_text("contig_name\tsample_name\tposition\tclipping_coverage\n")
    try:
        filt.process_assembly(a, None, m, t / "out.fa", "split", 2500, False)
    except SystemExit as exc:
        assert exc.code == 1
        return "File clipping valido senza eventi: filter-misassemblies termina con exit 1"
    raise AssertionError("Expected exit 1")


def anvi_clipping_schema(t):
    # Schema ufficiale anvi'o 9: contig length pos relative_pos cov clipping clipping_ratio.
    f = t / "MisAsm-clipping.txt"
    f.write_text("contig\tlength\tpos\trelative_pos\tcov\tclipping\tclipping_ratio\n"
                 "contig_001\t1665603\t500999\t0.30079136504917436\t120\t120\t1.0\n")
    data = misassembly.load_clipping_data(t)
    filtered = misassembly.load_clipping_data(t, min_coverage=50)
    assert data['contig_001'] == [(500999, 1)] and not filtered
    return "120 reads clipped nello schema anvi'o 9 diventano coverage=1; soglia 50 rimuove l'evento"


def anvi_zero_schema(t):
    f = t / "MisAsm-zero_cov.txt"
    f.write_text("contig\tlength\trange\trange_size\n"
                 "contig_001\t1665603\t100500-101000\t500\n")
    assert zero_coverage.parse_zero_cov_files(t) == []
    return "Intervallo zero-cov ufficiale 100500-101000 scartato perché il parser attende due interi"


for name, fn in [
    ("F01_hash_mismatch", hash_mismatch), ("F02_no_large_output", no_large),
    ("F03_missing_checkm2_contig", missing_quality), ("F04_wrapper_missing_windows", wrapper_windows),
    ("F05_window_column_mismatch", window_schema), ("F06_long_read_overlap", long_read_stats),
    ("F07_unpaired_discordance", unpaired), ("F08_blast_failure_as_no_hits", blast_failure),
    ("F09_missing_rna_keeps_HQ", missing_rna), ("F10_stale_tier_files", stale_tiers),
    ("F11_unresolved_as_clean", unresolved), ("F12_missing_read_as_classified", missing_reads),
    ("F13_remove_wrong_contig", remove_by_bin_filename), ("F14_final_tier_mismatch", final_tier),
    ("F15_yaml_regex", yaml_regex), ("F16_extension_not_forwarded", extension_plan),
    ("D01_verdict_report_only", verdict_not_consumed), ("F17_verdict_schema", verdict_schema),
    ("F18_zero_clipping_events", zero_clips),
    ("F19_anvio_clipping_schema", anvi_clipping_schema), ("F20_anvio_zero_cov_schema", anvi_zero_schema),
]:
    probe(name, fn)

print(json.dumps({"repo": str(repo), "results": results,
                  "confirmed": sum(r["confirmed"] for r in results), "total": len(results)}, indent=2))
sys.exit(0 if all(r["confirmed"] for r in results) else 1)
