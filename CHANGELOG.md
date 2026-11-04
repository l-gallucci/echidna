# Changelog

All notable changes to Echidna will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.0.0] - 2026-11-04

First GitHub Release; distribution `echidnaxomics`, command/import `echidna`.
See [release notes](docs/releases/v1.0.0.md). Earlier version entries are retained
as project history.

### Fixed
- Preserve one-taxon clean signals when the foreign-fraction threshold is zero.
- Parse native anvi’o clipping counts and zero-coverage ranges by column name,
  preserving legacy schemas; distinguish no events from invalid/missing files.
- Match hashed CheckM2 names to original contigs; retain contigs with missing
  quality rows and write an identity map. Produce outputs with no large contigs.
- Pass assembly/extension through workflow steps; share the taxonomy-window
  column contract and write header-only outputs when needed.
- Compute clipping-window depth from aligned CIGAR bases, include overlapping
  long reads and avoid treating unpaired HiFi reads as discordant pairs.
- Propagate external process failures instead of interpreting them as no hits.
- Treat unresolved contig taxonomy and unmatched reads as missing evidence;
  require at least ten classified reads for read-level risk.
- Demote HQ without RNA evidence when the RNA gate is requested; refuse reused
  populated selection outputs rather than leaving obsolete tier files.
- Remove HQ contigs by exact FASTA identities instead of bin filenames.
- Preserve regex/backslash/quote values when writing YAML configuration.

### Added
- Reuse sequence-linked `classify-contigs` assignments to gate retained MAGs,
  trimmed subsets and dereplication-cluster alternatives by taxonomic coverage
  and foreign-base fraction, with final checks on the actual exported FASTAs.
- Add `classify-contigs --taxonomy-table` to recompute bin membership summaries
  without rerunning a classifier. Reject ambiguous raw contig-ID collisions;
  the MMseqs wrapper deduplicates identical shared IDs before classification.
- Add optional pre-QC overlap/taxonomy selection to `mag-quality-wf` and the
  Snakemake DAG; all downstream QC follows resolved genomes, including fresh
  CheckM2 and GUNC results. Setup flags record equivalent configuration.
- Optional shared-contig resolution via --exclude-mags / --remove-shared-contigs,
  explicit quality ranking, preserved originals and removed-record exports.
  dRep secondary clusters or generic cluster tables can supply whole-genome
  replacements from original candidates; candidates and final FASTAs are
  checked globally. Trimmed quality is invalidated and marked for reassessment.
- Standalone check-shared-contigs: exact whole-contig sharing across MAG FASTAs,
  including renamed sequences and reverse complements; separate conflicting-ID
  reports, pairwise overlap and per-MAG summaries. No external dependencies or
  automatic changes to bin membership.
- Distribution metadata renamed to echidnaxomics; the executable and import
  remain echidna, with no echidnaxomics command.
- Isolated source/wheel builds, bundled workflow templates exportable with
  setup-workflow --export-workflow, CFF citation metadata, CI checks and a local
  Bioconda recipe renderer. No release or submission is performed automatically.
- Snakemake execution tests with fixture upstream reports and real downstream
  commands, covering RNA on/off, paths with spaces and reruns across tier changes.
- Optional read-taxonomy roll-up via --taxonomy-rank/--nodes-dmp.
- final-report organized_tier via --organized-quality; explicit circular header
  annotations, contig counts and an independent near-complete candidate flag.
- Maintained wiki source and instructions for publishing the separate GitHub wiki.

### Compatibility notes
- SSU comparisons preserve GTDB clade suffixes (_A/_B); this can expose previously
  hidden disagreements. Existing sequence names and tier thresholds are unchanged.
- Assembly filters require per-contig read-risk tables, not per-MAG verdicts.
- --hq-list requires exact FASTA identifiers. Repeated organization needs a new
  output directory. Added report columns must be read by header.
- setup-workflow supports installed templates via --export-workflow; --repo
  still selects an editable checkout. Chimera screening remains report-only.
- Snakemake tracks individual FASTAs and owns the declared tool output directories
  and organized/Selected; reruns recreate them to avoid stale tiers.

## [0.0.4] - 2026-10-14

### Changed (BREAKING output)
- `check-chimeras` no longer scores points; it now emits a 3-state VERDICT from
  TWO signals only: GUNC (`pass.GUNC=False`) and per-contig taxonomy (MMseqs via
  `--contig-taxonomy`, foreign LENGTH fraction >= `--min-foreign-frac`, default
  0.05). both chimeric -> `Chimeric`; one -> `Check`; neither -> `Clean`; with a
  signal missing the available one is used (a single firing -> `Check`), else
  `Unknown`. Output column `chimera_risk` is replaced by `chimera_verdict`.
- CheckM2 contamination, GTDB-Tk taxonomy/warning and GC are REPORTED as columns
  only and no longer affect the verdict. GC is now a `gc_strong_outlier` column:
  a long contig whose GC deviates from the bin's length-weighted mean by more
  than `--gc-max-dev` (default 0.10), ignoring contigs shorter than
  `--gc-min-length` (default 5000 bp) - so a short high-GC contig no longer
  flags anything. `--gc-window`/`--gc-step` removed.
- `final-report`: the chimera column is now `chimera_verdict`; concerns use
  `chimeric` / `chimera_check`.

### Added
- `utils/chimera.strong_gc_outlier` and `assess_chimera_verdict`.

## [0.0.3] - 2026-10-14

### Changed
- `check-chimeras` scoring retuned so Medium means a real taxonomic/contamination
  concern, not accumulated weak signals. GTDB-Tk warning is now +1 (was +0.5),
  GC inter-contig +0.5/+0.25 (was +1/+0.5) and windowed GC +0.25 (was +0.5).
  Resulting tiers: contamination alone = Low; contamination + GTDB warning =
  Medium; GC + minor foreign + GTDB warning = Medium; the same with foreign
  >= 5% = High; GUNC fail or foreign >= 10% = High; GC alone = Low. A minor
  foreign contig still scores (+0.5 -> at least Low) so it is never hidden, and
  is always named in `foreign_contigs` regardless of tier.

## [0.0.2] - 2026-10-14

### Changed
- `check-chimeras` scoring reworked to be abundance-aware and less restrictive:
  - CheckM2 contamination only scores above 10 % (+1); <=10 % is tolerated (HQ).
  - GUNC is trusted as authoritative (already proportion-aware): `pass.GUNC=False`
    +3 and nothing else; `pass.GUNC=True & CSS>0.45` +1.
  - The old binary GUNC per-contig hard flag is removed. The per-contig signal
    now comes from classify-contigs/MMseqs via `--contig-taxonomy`, weighted by
    the foreign LENGTH fraction: >=10 % +3, >=5 % +2, >0 % +0.5 (a few genes /
    one short contig stays minor, likely HGT).
  - GC inter-contig +1/+0.5, windowed GC +0.5, GTDB-Tk warning +0.5.
  - Thresholds: High >=3, Medium >=2, Low >=0.5. New output column
    `contig_foreign_frac` (replaces `gunc_contig_clades`); `--gunc-tax-rank`/
    `--gunc-min-genes` removed.
- `classify-contigs` now weights the majority taxon and the foreign fraction by
  contig LENGTH; new `foreign_frac` column and `--min-foreign-frac` (default
  0.05) gating `chimera_flag` (foreign contigs are always listed).

### Added
- `utils/names.canonical_name`: shared genome-name normalisation (strip FASTA
  extension + optional `--name-strip` suffixes) for cross-tool joins.
  `check-chimeras` now joins all sources canonically and WARNS about bins with
  no match in a source (silent N/A made visible).
- `final-report`: `gunc_mmseqs_concordance` column (True/False/NA) = whether GUNC
  (pass.GUNC) and MMseqs (classify-contigs) agree on chimeric-or-not.

### Fixed
- `setup.py` read the version with a non-multiline regex, so the distribution
  metadata was always `0.0.0`; it now reads `echidna/__init__.py` correctly.

## [0.0.1] - 2026-10-14

First public release of Echidna (developed earlier in-house as MaMISA; the
package, import name and CLI command are now `echidna`).

### Added
- Misassembly cleaning (Branch A, PacBio HiFi): `run-anvi-misassemblies`
  (wrapper for `anvi-script-find-misassemblies`, `--just-do-it` on by default),
  `process-large-contigs`, `check-read-chimeras`, `classify-clipping`,
  `filter-misassemblies`.
- MAG quality + taxonomy (Branch B): `run-checkm2`, `run-gtdbtk`, `run-gunc`
  (with `--contig-taxonomy-output`), `run-mimag-rna` (with `--outseq`),
  `organize-mags` tiering/renaming.
- Chimera / integrity screening: `check-chimeras` (GC + CheckM2 contamination +
  GTDB-Tk + GUNC, incl. GUNC per-contig multi-clade), `classify-contigs`
  (per-contig taxonomy via Kraken2/Kaiju + NCBI nodes.dmp, or MMseqs2 via GTDB
  lineage; names the minority/foreign contigs like MAGpurify2),
  `run-mmseqs-contigs`, `check-ssu-taxonomy` (16S vs GTDB-Tk, blastn).
- `final-report`: one per-MAG table joining CheckM2, GTDB-Tk, MIMAG rRNA/tRNA,
  check-chimeras, classify-contigs and check-ssu-taxonomy, with a `concerns`
  column.
- Chained workflows `misassembly-wf` / `mag-quality-wf` (with `--dry-run`,
  `--from`, `--skip`), a Snakemake `all_mag_quality` target, and
  `docs/workflows.md` with data-flow diagrams.

### Notes
- MMseqs2 lineage in `classify-contigs` is located by content (the last
  ';'-separated field), robust to the variable column count of
  `mmseqs createtsv --tax-lineage`.
- `setup.py` reads the version from `echidna/__init__.py` (single source of
  truth).
- `python_requires >= 3.8` for older tool envs (e.g. GTDB-Tk).

## Earlier (as MaMISA)

Initial internal development under the name MaMISA
- `filter-misassemblies` command for intelligent misassembly handling
- `remove-hq-contigs` command for HQ genome removal
- `filter-checkm2` command for quality-based genome filtering
- `run-gtdbtk` command for taxonomy classification
- Comprehensive documentation and examples
- Streaming FASTA processing for memory efficiency
- Dry-run mode for all commands
- Detailed statistics reporting

### Features
- HQ genome preservation with misassembly detection
- Contig splitting at misassembly positions
- CheckM2 integration with MiMAG-compliant tiers
- GTDB-Tk workflow automation
- Multi-format ID extraction and matching
- Progress reporting and logging

### Documentation
- Complete README with workflow examples
- Installation guide
- Command reference
- Citation information

## Roadmap

- HTML report generation
- Multi-threading support for large assemblies
- Integration with additional QC tools
- Docker container
- Test suite
- Continuous integration
