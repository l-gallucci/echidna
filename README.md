<p align="center">
  <img src="https://raw.githubusercontent.com/l-gallucci/echidna/main/assets/LOGO2.jpg" alt="Echidna Logo" width="400">
</p>

# Echidna - metagenomic assembly, MAG integrity and chimera screening

A dependency-free Python toolkit for metagenomic assembly quality control: from a
raw assembly to quality-filtered, taxonomically classified genomes.

Echidna itself is pure Python (no dependencies) and installs in seconds. It
orchestrates standard tools (anvi'o, CheckM2, GTDB-Tk, GUNC, barrnap,
tRNAscan-SE, Kraken2, samtools, BLAST+), each in its own conda env.

## Quickstart

```bash
git clone https://github.com/l-gallucci/echidna.git
cd echidna
pip install -e .
echidna --help
```

That covers the dependency-free commands. Heavy commands additionally call
external tools, each installed in its own conda env - see
[docs/installation.md](docs/installation.md).

## Commands at a glance

Assembly QC and misassembly handling:

| Command | What it does |
|---|---|
| `run-anvi-misassemblies` | Detect soft-clipping misassembly positions (anvi'o wrapper, PacBio HiFi) |
| `process-large-contigs` | Extract, QC and filter large contigs before binning |
| `classify-clipping` | Label each clipping position (repeat / deletion / chimeric_join / ...) from BAM evidence |
| `filter-misassemblies` | Split or remove contigs at misassembly positions |
| `remove-hq-contigs` | Remove known HQ-genome contigs from an assembly |
| `check-zero-coverage` | Validate no-coverage regions via BLAST + k-mers |

Chimera detection:

| Command | What it does |
|---|---|
| `check-read-chimeras` | Chimeras from read-level taxonomy (Kraken2 or Kaiju + BAM) |
| `check-chimeras` | Per-bin chimera verdict from GUNC + per-contig taxonomy (CheckM2 / GTDB-Tk / GC reported) |
| `run-gunc` | Gene-level chimerism/contamination (GUNC wrapper) |
| `run-mmseqs-contigs` | MMseqs2 taxonomy on bin contigs -> per-contig TSV for `classify-contigs` |
| `classify-contigs` | Flag bins whose contigs span >=2 taxa (per-contig taxonomy foreign fraction) |
| `check-ssu-taxonomy` | Check each MAG's 16S agrees with its GTDB-Tk taxonomy (blastn) |

MAG quality and taxonomy:

| Command | What it does |
|---|---|
| `run-checkm2` | Completeness/contamination (CheckM2 wrapper) |
| `run-gtdbtk` | Taxonomy classification (GTDB-Tk wrapper) |
| `run-mimag-rna` | MIMAG rRNA (5S/16S/23S) + tRNA criterion (barrnap + tRNAscan-SE) |
| `organize-mags` | Split genomes into HQ/MQ/LQ, rename `<sample>__<taxon>__<orig>`, optional full-MIMAG HQ gate |
| `final-report` | Join per-MAG signals (CheckM2, GTDB-Tk, MIMAG, chimera, 16S) into one table |

One-command workflows (run a whole branch, each step in its own env):

| Command | What it does |
|---|---|
| `misassembly-wf` | Whole misassembly-cleaning branch (anvi -> process -> readchim -> classify -> filter) |
| `mag-quality-wf` | Whole MAG quality + taxonomy branch (checkm2 -> organize -> gtdbtk -> gunc -> mimag -> chimeras) |

Setup and environment helpers:

| Command | What it does |
|---|---|
| `setup-workflow` | Configure the Snakemake workflow (env yamls + config), optionally build envs |
| `fetch-databases` | Register existing databases, or download missing ones |
| `check-envs` | Report tool versions per conda env + Echidna compatibility |

Every subcommand has `-h/--help`.

## Documentation

- [docs/installation.md](docs/installation.md) - install, the one-env-per-tool
  model, external tools, databases, tests
- [docs/workflows.md](docs/workflows.md) - the two branches (misassembly cleaning,
  MAG quality), per-step input/output and file types, with data-flow diagrams
- [docs/pipeline.md](docs/pipeline.md) - **Quick start** (run everything yourself,
  copy-paste), run options (standalone, `*-wf` wrappers or Snakemake), complete
  workflow, step-by-step guide
- [docs/hpc.md](docs/hpc.md) - running on a SLURM cluster
- [docs/command-reference.md](docs/command-reference.md) - every flag of every
  command
- [docs/references.md](docs/references.md) - citations (Trigodet et al. 2025 for
  misassemblies, MIMAG, and the driven tools)
- [workflow/README.md](workflow/README.md) - the Snakemake workflow

## Citation

If you use Echidna, please also cite the underlying methods you used - see
[docs/references.md](docs/references.md). In particular, misassembly detection is
built on anvi'o `anvi-script-find-misassemblies` and the signals described in
Trigodet et al. 2025.

## License

MIT License - see LICENSE.

## Contributing

Fork, branch, open a pull request. Issues and feature requests:
https://github.com/l-gallucci/echidna/issues
