# Changelog

All notable changes to Echidna will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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
