# Echidna

![Echidna logo](img/logo.jpg){ .project-logo }

**Metagenomic assembly, MAG integrity and chimera screening.**

Echidna is a Python command-line toolkit for assessing metagenomic assemblies
and metagenome-assembled genomes. It coordinates established bioinformatics
tools and combines their outputs into quality and integrity reports.

The installable distribution is **`echidnaxomics`**. The command and Python
import are **`echidna`**.

## Start here

1. [Install Echidna and the external tools](installation.md).
2. [Choose an assembly or MAG workflow](workflows.md).
3. Follow the [pipeline tutorial](pipeline.md), or use the
   [command reference](command-reference.md).
4. Read the [guide to interpreting results](guide.md) (Italian).

## Two workflows

| Workflow | Input | Purpose |
|---|---|---|
| Assembly integrity | Assembly FASTA and HiFi reads mapped to the same assembly | Detect misassembly signals and prepare an assembly for binning |
| MAG quality and taxonomy | One genome FASTA per bin | Assess completeness, contamination, taxonomy, RNA and chimera signals |

Binning is performed externally. After splitting an assembly, remap reads before
computing depth for a new binning run.

## Install from source

```bash
git clone https://github.com/l-gallucci/echidna.git
cd echidna
python -m pip install -e .
echidna --help
```

Until a release is published, install from this checkout. The registry package
named `echidna` belongs to an unrelated project.

The Python core uses the standard library. Commands that run CheckM2, GTDB-Tk,
GUNC, anvi'o or other tools also need their separate environments and databases;
see [installation](installation.md). The supplied [Snakemake workflow](snakemake.md)
starts from existing bins and covers part of the MAG workflow.

## Shared contigs between MAGs

Check an existing MAG collection for shared complete contig sequences and
conflicting FASTA identifiers with `echidna check-shared-contigs`.
See [examples and report definitions](shared-contigs.md). This standalone
check produces reports by default. Optional quality-based resolution exports
a new MAG set, with same-cluster replacements when original candidates are
available; the input FASTAs remain unchanged.

## Reading the results

Quality tiers are Echidna's operational thresholds and are distinct from strict
MIMAG categories. Chimera screening reports evidence for review; it does not
automatically exclude MAGs from quality folders. A `Clean` result means no
anomaly in the signals measured, rather than proof that a genome is error-free.
See the [interpretation guide](guide.md) for missing evidence, circularity and
reference-database concordance.

## Project

- [Source code and issues](https://github.com/l-gallucci/echidna)
- [References](references.md)
- [Citation and release information](releasing.md)
- [Changelog](changelog.md)

Echidna is distributed under the [MIT license](https://github.com/l-gallucci/echidna/blob/main/LICENSE).
