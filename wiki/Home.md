# Echidna — interpreting results

A guide to metagenomic assembly integrity, MAG quality assessment and chimera screening.
This page describes Echidna 1.0.0, prepared as the first GitHub Release.
See [GitHub Releases](https://github.com/l-gallucci/echidna/releases) for published
versions and the [1.0.0 release notes](../docs/releases/v1.0.0.md) for features
and validation.

## Two branches

| Branch | Input | Purpose |
|---|---|---|
| Assembly | Assembly and HiFi reads mapped back to that same assembly | Find misassemblies, process long contigs and prepare for binning |
| MAG | One FASTA per bin | Assess completeness, contamination, taxonomy, RNA and chimera signals |

The core uses only the Python standard library. Commands that call anvi'o,
CheckM2, GTDB-Tk, GUNC, MMseqs2, barrnap, tRNAscan-SE, samtools or BLAST need
those tools and their databases installed separately. The project's conda
environments separate their dependencies.

The current automatic workflow setup installs Echidna in each environment that
runs an Echidna wrapper. From a checkout:

```bash
python -m pip install -e .
echidna --version
echidna --help
```

The declared Python requirement for the core is Python >=3.8. Check external
tool compatibility separately. Local tests and CI configuration do not certify
every operating system or external tool environment.

## Identifiers and formats

Use unique MAG and contig names. BAM and FASTA files must refer to the same
reference and use matching identifiers. Read IDs must match the classifier
output. Complete the checks before renaming MAGs, and preserve the original
identifiers used in the reports.

The parsers recognize both native anvi'o 9 schemas and legacy Echidna columns:

```text
clipping: contig length pos relative_pos cov clipping clipping_ratio
zero-cov: contig length range range_size
```

In the native format, `clipping` is the number of clipped reads;
`relative_pos` is not coverage. `range` is a zero-based, half-open interval.
Header-only files are valid results with no events; missing or malformed files
are errors. See the [official anvi'o format](https://anvio.org/help/9/programs/anvi-script-find-misassemblies/).

`process-large-contigs` writes `contig_name_map.tsv` to link hashed CheckM2
names to original contigs. A contig without a CheckM2 row is retained as
unassessed. `remove-hq-contigs --hq-dir` reads contigs from MAG FASTAs; numbers
in bin filenames are not contig identifiers. Use exact FASTA IDs with `--hq-list`.

## Assembly branch

Dependency order:

```text
original assembly + BAM of the same reads
  → anvi'o: clipping and zero coverage
  → process-large-contigs
  → check-read-chimeras, when read classifications are available
  → classify-clipping: report-only annotation
  → filter-misassemblies
  → new mapping, depth calculation and external binning
```

Preview the wrapper without running external tools:

```bash
echidna misassembly-wf --assembly assembly.fa --bam hifi.bam \
  --output-dir assembly_qc --dry-run
```

The wrapper passes the FASTA to the read-level check, whose producer writes
a window table even when no windows are reported. `classify-clipping` reads
the resulting `is_taxon_shift` field. Its report does not automatically decide
where to split: the filter uses the original clipping sites, optionally
selected with `--min-clip-coverage`.

Local depth is a per-base window average calculated from aligned CIGAR blocks,
rather than a count of read starts. Long reads starting outside the window
are included if they span it. Pair statistics are N/A for unpaired reads such
as HiFi. Primary-read SA tags provide the join signal without counting
supplementary alignments twice.

After splitting, for example `contig_A → contig_A_split1 + contig_A_split2`,
the old BAM still refers to `contig_A`. Remap reads to the modified FASTA before
calculating depth for new binning. The original BAM remains appropriate for
checks against the original reference.

`HQ_extracted/` and `assembly_for_filtering.fa` are complementary outputs:
extracted contigs are absent from the residual assembly. Preserve and assess
both. An assembly without long contigs still produces the assembly needed by
the next step.

## MAG branch

```bash
echidna mag-quality-wf --bins-dir bins --genome-ext fa \
  --output-dir mag_qc --tiers HQ,MQ \
  --gtdbtk-data /data/gtdbtk_release \
  --gunc-db /data/gunc/reference.dmnd \
  --mmseqs-db /data/mmseqs_gtdb/gtdb \
  --split-by-domain --dry-run
```

Remove `--dry-run` after preparing the environments, inputs and databases.
The default step order is:
`checkm2 → organize1 → gtdbtk → gunc → [mmseqs → contigtax] → mimag → chimeras → organize2`.

LQ is not carried forward with the default `--tiers HQ,MQ`; select the desired
tiers explicitly. Extensions are passed to the MIMAG and chimera steps too.
Names remain unchanged until final organization.

Chimera checks are **report-only**: flagged MAGs can still appear in quality
folders. Choose explicitly which to exclude and which to review.
`check-ssu-taxonomy` and `final-report` remain separate commands. The Snakemake
DAG can reuse a per-contig taxonomy table through `classify-contigs`, but does
not perform the initial MMseqs classification.

## Shared contigs between MAGs

The standalone `check-shared-contigs` command checks whether two or more MAGs
contain the same complete contig sequence, including renamed contigs and
reverse complements:

```bash
echidna check-shared-contigs --bins-dir bins/ --output-dir overlap/
```

Add `--recursive` for nested HQ/MQ/LQ folders. Multiple directories can follow
`--bins-dir`. Reports distinguish shared sequences from identical IDs assigned
to different sequences, and include MAG pairs and shared fractions. By default,
FASTA files, membership and chimera verdicts are unchanged. Identical sequences
can represent repeated elements; sharing alone is not proof of contamination.

With `--quality-report`, choose `--exclude-mags` or `--remove-shared-contigs`
to export resolved copies into `retained/`, with excluded MAGs and removed
contigs saved separately. Default priority is decreasing completeness, then
increasing contamination. `--rank-by quality-score` instead uses completeness
minus five times contamination.

With `--drep-dir` (or `--cluster-table`) and `--candidates-dir`, the command
looks for alternatives from the same cluster. Each candidate is checked against
the entire retained set, and the final collection is checked again. Trimmed
MAGs have unassessed final quality and need new CheckM2 results; replacements
need their own downstream checks. Original FASTAs remain unchanged.
See [comparison and resolution rules](../docs/shared-contigs.md).

### Cluster recovery with taxonomy checks

`check-shared-contigs --contig-taxonomy <classify-contigs.tsv>` also checks
the internal taxonomic coherence of retained MAGs and recovered candidates.
Supply the **per-contig** table containing sequence fingerprints, lengths
and rank, not the `.bins.tsv` summary. Candidates must pass both the sharing
and taxonomy checks. If one fails, the next candidate is tried. The command
does not construct hybrid MAGs by joining contigs from different genomes.

By default, at least 50% of all MAG bases must resolve at the selected rank;
minority taxa must account for less than 5% of resolved bases. Both thresholds
are configurable. Insufficient evidence does not mean coherence. Sequence
matching allows reuse after renaming or reverse complementation; conflicting
assignments become unresolved.

`mag-quality-wf --resolve-overlaps exclude-mags` (or `remove-contigs`) performs
this selection before CheckM2, GTDB-Tk, GUNC and RNA screening. It requires
`--quality-report` for ranking and `--contig-taxonomy` for verification.
`--drep-dir`/`--cluster-table` with `--candidates-dir` enables replacements.
Final bin taxonomy is reaggregated from per-contig assignments without
rerunning the classifier. The Snakemake DAG offers the same optional selection
through the `overlap_*` configuration keys.

To add fingerprints to older outputs, rerun only `classify-contigs` using the
existing raw classifier output and its corresponding original FASTAs.
See the [procedure and workflow](../docs/workflows.md#overlap-resolution-with-existing-contig-taxonomy).

## Consistent databases and read comparisons

When comparing GUNC, MMseqs2, GTDB-Tk and SSU signals, prefer references derived
from the **same taxonomy and release**, where available. A common taxonomy
reduces differences in nomenclature and coverage; it does not make the two
methods independent.

The same `.dmnd` file cannot be passed to MMseqs2: each tool requires its own
formats and preparation. Use an appropriate GUNC release and a corresponding
MMseqs database, recording release, provenance and checksums. Running both
steps supplies the two signals; running the same classifier twice does not
provide corroboration. Alternative databases can be used in separate runs
for sensitivity analysis.

For reads, comparison at a common rank is available:

```bash
echidna check-read-chimeras --bam hifi.bam --assembly assembly.fa \
  --kraken2-output reads.kraken2 --kraken2-report reads.report \
  --taxonomy-rank order --nodes-dmp /data/kraken/taxonomy/nodes.dmp \
  --output-dir reads_qc
```

`nodes.dmp` must come from the classifier's database. Without `--taxonomy-rank`,
the original taxids are compared. Rank roll-up reduces false conflicts between
a species taxid and one assigned to its genus/order. At least ten classified
reads are required for read-level risk assessment. Reads missing from the
classifier do not become an artificial taxon. Windows assess classified reads.

Useful additional evidence includes consistent SA support, mapping identity
and quality, support on both sides of a junction, and inspection in a viewer.
Mixed taxonomy is a screening signal, not definitive proof of a chimera.

## Quality, circularity and GTDB suffixes

The internal thresholds intentionally remain:

| Tier | Completeness | Contamination |
|---|---|---|
| HQ | >=90% | <=5% |
| MQ | >=70% | <=10% |
| LQ | >=50% | <=10% |

The default HQ extraction threshold in `process-large-contigs` is instead
50%/10%. These operational labels are distinct from the MIMAG standard.
With `--mimag-rna-dir`, an HQ MAG without complete RNA evidence, including a
missing result, is demoted to MQ. The RNA value remains blank when unassessed.

The report distinguishes:

- `circular_annotated`: explicit FASTA annotation for a single-contig MAG;
  N/A when absent. A `circular=true` description is preserved during long-contig
  extraction. `circular=false` overrides the filter's length heuristic.
- `near_complete_candidate`: completeness >=99% and contamination <=5%,
  independently of circularity. This highlights quality; it does not declare
  a complete genome.
- `n_contigs`: FASTA contig count. These structural fields require
  `final-report --genomes-dir`.

In assembly filtering, the length heuristic remains a fallback for
`--split-hq-circular` when an annotation is absent. It is not proof of
circularity; assembler definitions can differ.

SSU comparisons **preserve `_A`, `_B` and other GTDB suffixes**. This correction
does not rename FASTA records. A comparison that previously agreed only because
a suffix was removed can now disagree; interpret it against the reference
release. Do not assume automatic SILVA/GTDB name equivalence: prefer consistent
references or a curated taxonomy map.

## Verdicts and final report

| Evidence | MAG verdict |
|---|---|
| GUNC fail and foreign fraction >=0.05 | Chimeric |
| Only one signal flags an anomaly | Check |
| No anomaly in the measured signals | Clean |
| No measured signal | Unknown |

The foreign fraction is weighted by the length of **resolved** contigs. Zero
resolved contigs means missing evidence, not a measured fraction of zero.
Legacy tables without `n_contigs_resolved` remain accepted for compatibility,
but cannot support that availability check. Prefer current `classify-contigs`
outputs.

GC, CheckM2 contamination and GTDB warnings are informative columns. One
measured, non-flagging signal can yield Clean; read this as "no anomaly in that
signal." The consensus rule is unchanged. The table above uses the default
foreign-fraction threshold; `--min-foreign-frac` can change it.

```bash
echidna final-report --genomes-dir bins \
  --checkm2-report mag_qc/01_checkm2/quality_report.tsv \
  --gtdbtk-dir mag_qc/03_taxonomy --mimag-rna-dir mag_qc/05_mimag \
  --chimera-report mag_qc/06_chimera/chimera_report.tsv \
  --organized-quality mag_qc/07_final/merged_quality.tsv \
  --output mag_qc/final_report.tsv
```

`tier` retains its historical meaning: default CheckM2 thresholds.
`organized_tier` reports the effective selection, including custom thresholds
and RNA gating, when `--organized-quality` is supplied. Identify report genomes
by their original names rather than renamed final-output filenames. If overlap
resolution was used, point `--genomes-dir` to `mag_qc/00_overlap/retained/`,
so structural fields describe the resolved collection.

Added columns do not rename existing ones; TSV consumers should read headers.
An empty `concerns` field does not certify MAG quality. BLAST retains the
Trigodet/anvi'o workflow criterion: a successful search without hits can yield
`unsupported`, while execution errors stop the command.

## Resuming and validation

Use new output directories for runs with changed inputs or thresholds.
To prevent stale tiers, `organize-mags` refuses a populated selection directory;
it does not delete existing files. `--from` assumes valid outputs from earlier
steps. To repeat final organization, use a new output or first move the previous
selection aside.

The distribution is named **`echidnaxomics`**; the command and import remain
**`echidna`**. Install from source until this distribution is available through
your chosen registry; the registry package named `echidna` is unrelated.
Wheels include Snakemake templates: `echidna setup-workflow --export-workflow ./my-workflow`
copies them into a new writable directory. `--repo /path/to/checkout` remains
available. With an exported workflow, pass an explicit `--config` path to
`fetch-databases`. Environment templates must be checked against the tool
releases you use.

The DAG tracks individual FASTAs and tool output directories. On reruns,
Snakemake recreates declared outputs, including `organized/Selected`, to avoid
stale tiers. Keep only workflow outputs in those directories. DAG tests use
simulated upstream reports and real downstream Echidna commands; they are not
validation against real databases. Author metadata, DOI archiving and the
Bioconda recipe are described in the [release guide](../docs/releasing.md).

```bash
PYTHONDONTWRITEBYTECODE=1 python -m pytest -p no:cacheprovider
```

Tests include clean controls, native-format handoffs and an optional smoke test
using a synthetic BAM with real samtools. They do not replace a full pipeline
run with installed databases. The historical correction record and retained
methodological choices are in the [repository report](https://github.com/l-gallucci/echidna/blob/main/docs/audit-2026-11-03/FIXES.md).
