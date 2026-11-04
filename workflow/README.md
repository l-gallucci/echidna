# Echidna Snakemake workflow

Runs the MAG quality-control tail of Echidna with **one conda env per rule**, so
CheckM2, GTDB-Tk, GUNC and the rRNA/tRNA tools (which pin conflicting
dependencies) never share an env. Snakemake builds each env automatically on
first run.

## Install Snakemake

```bash
conda install -c conda-forge -c bioconda snakemake-minimal
# or: pip install snakemake
```

## Configure

The quick way — let Echidna do the edits:

```bash
echidna setup-workflow --genomes-dir bins/ --genome-ext fa --threads 40 --install
echidna fetch-databases --gtdbtk-data /data/gtdbtk_r220 \
                       --gunc-db /data/gunc/gunc_db_progenomes2.1.dmnd   # or --download
```

`setup-workflow` points each `envs/*.yaml` at this checkout (`pip -> -e <repo>`)
and writes the paths/threads into `config.yaml`; `--install` builds every
per-rule conda env. `fetch-databases` registers existing databases into the
config, or downloads missing ones with each tool's own downloader.

Or edit by hand: set `genomes_dir`, `genome_ext`, the database paths
(`gtdbtk_data`, `gunc_db`), `threads`, `sample_regex` / `tax_level` in
`config.yaml`, and in each `envs/*.yaml` replace `- echidnaxomics` with
`-e /absolute/path/to/this/repo` until Echidna is on PyPI.

## Run

Two modes, set by `use_named_envs` in config.yaml:

```bash
# mode A — Snakemake builds one conda env per rule from envs/*.yaml
snakemake --use-conda --cores 20 -s workflow/Snakefile --configfile workflow/config.yaml

# mode B — use YOUR existing named envs (use_named_envs: true, set by
#          `setup-workflow --gunc-env ... --checkm2-env ...`); NO --use-conda,
#          each rule is wrapped in `conda run -n <env>`
snakemake --cores 20 -s workflow/Snakefile --configfile workflow/config.yaml

# preview the plan without running
snakemake -n -s workflow/Snakefile --configfile workflow/config.yaml
```

## DAG

```
genomes_dir ─┬─> run-checkm2 ───┐
             ├─> run-gtdbtk ────┤
             ├─> run-mimag-rna ─┼─> organize-mags (HQ/MQ/LQ + rename + MIMAG rRNA/tRNA gate)
             └─> run-gunc ──────┘
                      └──────────> check-chimeras (GUNC + per-contig taxonomy)
```

`run-mimag-rna` (barrnap + tRNAscan-SE) is optional: set `mimag_rna: false` in
`config.yaml` to skip it, and HQ is then based on completeness/contamination
only. Binning and anvi'o misassembly detection are assumed done upstream; this
workflow starts from a directory of MAG FASTAs.

## Outputs (under `results/`)

| Path | Produced by |
|---|---|
| `checkm2/quality_report.tsv` | run-checkm2 |
| `gtdbtk/` | run-gtdbtk |
| `gunc/` | run-gunc |
| `mimag_rna/mimag_rna_summary.tsv` | run-mimag-rna (if `mimag_rna: true`) |
| `organized/Selected/{HQ,MQ,LQ}/` | organize-mags (renamed genomes) |
| `organized/merged_quality.tsv` | organize-mags |
| `chimera_report.tsv` | check-chimeras |


Current scope: the supplied Snakefile can reuse a sequence-linked per-contig
taxonomy table through `classify-contigs`; it does not generate raw MMseqs
classifications, run SSU checks or final-report. Without the optional table its
chimera verdict uses GUNC alone. Organized MAGs are not automatically excluded
by that verdict. MIMAG uses the configured fixed kingdom; mixed domains need
separate handling.

## Optional overlap resolution and taxonomy reuse

```bash
echidna setup-workflow --genomes-dir drep/dereplicated_genomes/ \
  --resolve-overlaps exclude-mags --quality-report checkm2_all/quality_report.tsv \
  --contig-taxonomy all_contigs.classified.tsv --taxonomy-rank order \
  --drep-dir drep/ --candidates-dir all_original_mags/
```

Selection uses existing completeness/contamination for ranking, tests each
candidate against every retained MAG for shared full contigs, and requires
sufficient taxonomic evidence and coherence. It writes under
`<outdir>/overlap_resolution/`. All downstream CheckM2, GTDB-Tk, GUNC and optional
RNA analyses run on the resolved `.fa` collection. Bin-level taxonomy fractions
are recalculated from the per-contig evidence for those genomes.

Set `overlap_resolution` to `remove-contigs` for shared-copy removal, or blank
to disable selection. `overlap_cluster_table` can replace `overlap_drep_dir` for
other tools. `min_taxonomy_resolved_frac` defaults to 0.5 of all MAG bases and
`min_foreign_frac` to 0.05 of resolved bases. The latter threshold is also passed
to check-chimeras. A checkpoint tracks changed representatives, candidates and
tables and discovers the resulting FASTA names before scheduling downstream QC.
If no MAG passes, inspect the selection reports; downstream jobs stop.

`contig_taxonomy` alone can also enable the second chimera signal without any
overlap selection. Use the per-contig `classify-contigs` table with sequence
fingerprints, **not `.bins.tsv`**. Older tables can be regenerated from existing
classifier output without rerunning the classifier. See
[taxonomy reuse](../docs/shared-contigs.md#reuse-contig-taxonomy-when-selecting-replacements).


## Installed distribution and reruns

The distribution is `echidnaxomics`; the command remains `echidna`. Wheels ship
the templates. Run `echidna setup-workflow --export-workflow ./my-workflow` to
export a writable copy to a new directory. After publication, exported envs
install `echidnaxomics`; before publication use `--repo` with the checkout to
point their pip entries at local source. No new command name is introduced.

Default config is resolved beside the Snakefile, so it can be launched from a
different working directory. Relative input/output paths still refer to that
working directory. FASTA inputs are tracked individually, including content
changes that do not change their parent directory's mtime. Snakemake owns the
declared tool output directories and `organized/Selected`, and removes them
before rerunning their producers. Do not keep unrelated files there. Standalone
`organize-mags` still refuses a populated selection rather than deleting it.

Integration tests execute Snakemake with fixture upstream reports and real
`organize-mags`/`check-chimeras` commands. They cover RNA on/off, paths with
spaces, HQ-to-MQ tier changes and reruns after editing an input FASTA; they do
not execute the expensive scientific tools or their databases.
