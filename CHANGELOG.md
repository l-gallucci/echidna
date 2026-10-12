# Changelog

All notable changes to Echidna will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [2.0.0] - 2026-10-12

### Changed (BREAKING)
- Renamed the project MaMISA -> Echidna. The package, import name, CLI command
  and the suggested conda env are now `echidna` (were `mamisa`), and the wf flag
  is `--echidna-env`. No backward-compatible `mamisa` alias is kept; reinstall
  with `pip install -e .` in each env.

### Added
- `run-anvi-misassemblies`: wrapper for `anvi-script-find-misassemblies`
  (`--just-do-it` on by default).
- `misassembly-wf` and `mag-quality-wf`: chained per-branch workflow commands
  (like `gtdbtk classify_wf`), with `--dry-run`, `--from`, `--skip`; plus a
  Snakemake `all_mag_quality` target and `docs/workflows.md` with data-flow
  diagrams.
- `classify-contigs`: per-contig taxonomy check independent of GUNC
  (Kraken2/Kaiju via NCBI nodes.dmp, or MMseqs2 via GTDB lineage). Flags bins
  whose contigs span >= 2 taxa at a chosen rank and names the minority (foreign)
  contigs, like MAGpurify2's taxonomy module.
- `run-mmseqs-contigs`: runs MMseqs2 taxonomy on bin contigs (adds
  `--tax-lineage 1`) to feed `classify-contigs` with a GTDB-consistent lineage.
- `check-ssu-taxonomy`: classifies each MAG's 16S/SSU (blastn) and flags
  disagreement with the GTDB-Tk genome taxonomy; `run-mimag-rna --outseq`
  extracts the 16S.
- `final-report`: one per-MAG table joining CheckM2, GTDB-Tk, MIMAG rRNA/tRNA,
  check-chimeras, classify-contigs and check-ssu-taxonomy, with a `concerns`
  column.

### Changed
- `check-chimeras`: GUNC per-contig multi-clade is now a hard flag (risk forced
  to High at any `--gunc-tax-rank`); GC weights rebalanced so GC never reaches
  Medium alone. `run-gunc` gains `--contig-taxonomy-output`.
- `python_requires` lowered to `>= 3.8` for older tool envs (e.g. GTDB-Tk).

## [1.0.0] - 2024-XX-XX

### Added
- Initial release (as MaMISA)
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
