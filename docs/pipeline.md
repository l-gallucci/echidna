# Pipeline

Binning and anvi'o misassembly detection are assumed done upstream for the
MAG-quality tail; the misassembly branch starts from an assembly plus a BAM.

For the two branches at a glance - which step needs which input file, and the
data-flow diagrams - see [workflows.md](workflows.md). The misassembly branch is
PacBio HiFi only.

## Quick start - run everything yourself

The fast path if you want full control and no Snakemake. Each heavy tool runs in
its own conda env (`conda run -n <env> ...`); set the envs up once
([installation.md](installation.md)) and point Echidna at databases you already
have. Fill in the paths and paste.

```bash
# --- 0. one-time: databases you already downloaded (validated, written to config)
echidna fetch-databases \
    --gtdbtk-data /data/gtdbtk_r220 \
    --gunc-db     /data/gunc/gunc_db_gtdb214.dmnd \
    --checkm2-db  /data/checkm2/uniref100.KO.1.dmnd        # CheckM2: a .dmnd FILE
# don't have them? add --download (see installation.md "Databases")

ASM=assembly.fa          # your assembly
BAM=mapping.bam          # reads mapped to $ASM (sorted + indexed); HiFi for branch A
BINS=bins/               # your bins (after binning)
GTDBTK_DB=/data/gtdbtk_r220
GUNC_DB=/data/gunc/gunc_db_gtdb214.dmnd
MMSEQS_DB=/data/mmseqs_gtdb/gtdb
T=40

# --- Branch A: misassembly cleaning (PacBio HiFi only; skip if not HiFi) ---
conda run -n anvio-9 echidna run-anvi-misassemblies -b $BAM -o misasm/MisAsm -T $T
echidna process-large-contigs -a $ASM -m misasm/ -o 01_large/ --threads $T
# optional read-level chimera signal (needs Kraken2/Kaiju on the reads):
echidna check-read-chimeras --bam $BAM --assembly $ASM --kraken2-output reads.kraken2 \
    --kraken2-report report.txt -o 02_chimera/
echidna classify-clipping --bam $BAM -m misasm/ \
    --taxonomy-windows 02_chimera/chimera_read_windows.tsv -o 03_classified.tsv
echidna filter-misassemblies -a 01_large/assembly_for_filtering.fa -m misasm/ \
    -g 01_large/HQ_extracted/ -o 04_clean.fa --mode split --safe
# -> remap the source reads to 04_clean.fa, recompute depth, then bin into $BINS

# --- Branch B: MAG quality + taxonomy (starts from $BINS) ---
conda run -n checkm2 echidna run-checkm2 --genome-dir $BINS -o results/checkm2 --threads $T
conda run -n gtdbtk  env GTDBTK_DATA_PATH=$GTDBTK_DB \
                     echidna run-gtdbtk --genome-dir $BINS -o results/gtdbtk --cpus $T
conda run -n gunc    echidna run-gunc --genome-dir $BINS -o results/gunc \
                     --file-suffix .fa --threads $T --db-file $GUNC_DB
conda run -n mimag   echidna run-mimag-rna --genome-dir $BINS -o results/mimag_rna --threads $T
# per-contig taxonomy = the 2nd chimera signal (optional but recommended):
conda run -n mmseqs2 echidna run-mmseqs-contigs --bins-dir $BINS --db $MMSEQS_DB \
                     -o results/contig_taxonomy.tsv --threads $T
echidna classify-contigs --bins-dir $BINS --mmseqs-output results/contig_taxonomy.tsv \
                     -o results/contigs.classified.tsv                # writes .bins.tsv
echidna check-chimeras --bins-dir $BINS --gunc-dir results/gunc \
                     --contig-taxonomy results/contigs.classified.tsv.bins.tsv \
                     --checkm2-report results/checkm2/quality_report.tsv \
                     --gtdbtk-dir results/gtdbtk -o results/chimera_report.tsv
echidna organize-mags --checkm2-root results/checkm2 --genomes-dir $BINS \
                     --gtdbtk-dir results/gtdbtk -o results/filtered --rename --copy \
                     --mimag-rna-dir results/mimag_rna
# optional: one joined per-MAG table
echidna final-report --checkm2-report results/checkm2/quality_report.tsv \
                     --gtdbtk-dir results/gtdbtk --mimag-rna-dir results/mimag_rna \
                     --chimera-report results/chimera_report.tsv \
                     --contig-taxonomy results/contigs.classified.tsv.bins.tsv \
                     -o results/final_report.tsv
```

Prefer one command per branch? See Option C (`misassembly-wf` / `mag-quality-wf`)
below. Want the full DAG on a cluster? Option B (Snakemake). Per-step detail with
explanations is in the step-by-step guide further down.

## Run options

Both ways use the same conda envs (see [installation.md](installation.md)).

### Option A - run commands yourself (simplest)

Call each command in its env:

```bash
conda run -n checkm2 echidna run-checkm2 --genome-dir bins/ -o results/checkm2 --threads 40
conda run -n gtdbtk  env GTDBTK_DATA_PATH=/data/gtdbtk_r220 \
                     echidna run-gtdbtk  --genome-dir bins/ -o results/gtdbtk --cpus 40
conda run -n gunc    echidna run-gunc    --genome-dir bins/ -o results/gunc \
                     --file-suffix .fa --threads 40 --db-file /data/gunc/gunc_db_progenomes2.1.dmnd
conda run -n mimag   echidna run-mimag-rna --genome-dir bins/ -o results/mimag_rna --threads 40  # optional
conda run -n echidna  echidna organize-mags --checkm2-root results/checkm2 --genomes-dir bins/ \
                     --gtdbtk-dir results/gtdbtk -o results/filtered --rename --copy \
                     --mimag-rna-dir results/mimag_rna          # drop this flag to skip the rRNA/tRNA HQ gate
conda run -n mmseqs2  echidna run-mmseqs-contigs --bins-dir bins/ --db /data/mmseqs_gtdb/gtdb \
                     -o results/contig_taxonomy.tsv --threads 40                  # optional 2nd chimera signal
conda run -n echidna  echidna classify-contigs --bins-dir bins/ \
                     --mmseqs-output results/contig_taxonomy.tsv -o results/contigs.classified.tsv  # writes .bins.tsv
conda run -n echidna  echidna check-chimeras --bins-dir bins/ --gunc-dir results/gunc \
                     --contig-taxonomy results/contigs.classified.tsv.bins.tsv \
                     --checkm2-report results/checkm2/quality_report.tsv --gtdbtk-dir results/gtdbtk \
                     -o results/chimera_report.tsv
```

The verdict uses two signals: GUNC (`pass.GUNC`) and the per-contig taxonomy
foreign fraction from `classify-contigs` (`--contig-taxonomy`). Both chimeric ->
`Chimeric`, one -> `Check`, neither -> `Clean`. CheckM2 contamination, GTDB-Tk
and the GC outlier are reported as columns only. Give at least one of
`--gunc-dir` / `--contig-taxonomy` (with only GUNC, a failing bin is `Check`).

MIMAG HQ also requires 5S/16S/23S rRNA and tRNAs for at least 18 amino acids,
which CheckM2 does not check. `run-mimag-rna` adds that evidence and
`organize-mags --mimag-rna-dir` demotes HQ genomes that fail it to MQ. Skip it and
HQ is completeness/contamination only.

For a mixed bacteria + archaea bin set, barrnap and tRNAscan-SE need the right
kingdom per genome. Run it once with `--split-by-domain --gtdbtk-dir <taxonomy>`:
Echidna reads the GTDB-Tk per-domain summaries and applies `bac`/`arc` automatically
(no manual splitting of the bins):

```bash
conda run -n mimag echidna run-mimag-rna \
    --selected-dir results/filtered/Selected --tiers HQ \
    -o results/mimag_rna --threads 40 \
    --split-by-domain --gtdbtk-dir results/gtdbtk
```

### Option B - Snakemake orchestrates it

One command runs the whole DAG. Configure once, then run:

```bash
conda install -c conda-forge -c bioconda snakemake-minimal   # in base or its own env

# tell the workflow which envs + databases to use (sanity-checks them)
echidna setup-workflow --genomes-dir bins/ --genome-ext fa --threads 40 \
    --echidna-env echidna --checkm2-env checkm2 --gtdbtk-env gtdbtk --gunc-env gunc

# then run (named-env mode: Snakemake wraps each rule in `conda run -n`)
snakemake --cores 40 -s workflow/Snakefile --configfile workflow/config.yaml
```

Don't have the envs and want Snakemake to build them for you instead? Use
`echidna setup-workflow --genomes-dir bins/ --install` and run with
`snakemake --use-conda ...`. Details in
[`../workflow/README.md`](../workflow/README.md). For clusters see
[hpc.md](hpc.md).

### Option C - one-command branch wrappers (no Snakemake)

Each branch also has a single subcommand that runs all its steps in order, each
in the right env, resumable with `--from` / `--skip`. No Snakemake, no DAG file -
just one process that shells out per step.

```bash
# misassembly-cleaning branch (PacBio HiFi): anvi -> process -> readchim
# (only with Kraken2/Kaiju) -> classify -> filter
echidna misassembly-wf \
    --assembly assembly.fa --bam mapping.bam -o 01_misasm/ --threads 40 \
    --anvio-env anvio-9 --checkm2-env checkm2 \
    --kraken2-output kraken2_reads.txt --kraken2-report kraken2_report.txt   # optional readchim step

# MAG quality + taxonomy branch: checkm2 -> organize -> gtdbtk -> gunc ->
# mimag -> chimeras -> organize. Genome names are frozen after the first
# organize so GUNC / check-chimeras / run-mimag-rna all match.
echidna mag-quality-wf \
    --bins-dir bins/ -o 02_magqc/ --threads 40 --tiers HQ,MQ \
    --checkm2-env checkm2 --gtdbtk-env gtdbtk --gunc-env gunc --mimag-env mimag \
    --gtdbtk-data /data/gtdbtk_r220 --gunc-db /data/gunc/gunc_db_progenomes2.1.dmnd \
    --mmseqs-db /data/mmseqs_gtdb/gtdb          # optional: adds per-contig chimera signal
```

Both print their plan with `--dry-run`, resume a failed run with `--from <step>`,
and drop steps with `--skip <step>[,<step>]`. See
[command-reference.md](command-reference.md) for every flag (step names included).

## Complete workflow (misassembly branch to taxonomy)

```
assembly.fa
    |
    v
[1] anvi-script-find-misassemblies   (anvi'o)
    |                           Detect soft-clipping positions in the original BAM.
    |                           Produces *-clipping.txt files.
    v
[2] process-large-contigs       Separate very large contigs, run CheckM2 on each,
    |                           extract clean candidates, return updated assembly.
    v
[3] check-read-chimeras         (optional, recommended)
    |                           BAM + Kraken2 per-read taxonomy.
    |                           Flags contigs whose reads come from >=2 organisms.
    |
    +-------------------------->[3b] check-chimeras   (optional)
    |                                Per-bin verdict from GUNC + per-contig taxonomy.
    |
    v
[4] classify-clipping           (optional, recommended)
    |                           Label each clipping position from BAM evidence:
    |                           end_artefact / repeat_collapse / deletion_artefact /
    |                           chimera_candidate / chimeric_join / sv_candidate.
    |
    v
[5] filter-misassemblies        Split or remove misassembled contigs.
    |                           Preserves HQ genomes. Integrates chimera report.
    |
    v
[6] binning                     (external - MetaBAT2, MaxBin2, etc.)
    |
    v
[7] run-checkm2                 (CheckM2) completeness/contamination
    |
    v
[8] organize-mags               Split bins into HQ / MQ / LQ quality tiers.
    |
    v
[9] run-gtdbtk                  Taxonomic classification of selected genomes.
```

## Step-by-step guide

### Step 1 - Detect misassemblies with anvi'o (external)

```bash
anvi-script-find-misassemblies \
    -b mapping.bam \
    -o misassemblies/MisAsm \
    -T 40 --just-do-it
```

Produces `MisAsm-clipping.txt` in `misassemblies/`. See
[references.md](references.md) (Trigodet et al. 2025) for the method.

### Step 2 - Process large contigs

Very long contigs (> 300 kbp by default) often represent single complete genomes
and should be handled separately before binning. This command separates large and
regular contigs, runs CheckM2 on each large contig, classifies each as extract HQ,
keep for splitting, or low quality, and returns an updated assembly.

```bash
echidna process-large-contigs \
    --assembly assembly.fa \
    --misassemblies misassemblies/ \
    --output-dir 01_large_contigs/ \
    --max-length 300000 \
    --min-completeness 50 \
    --max-contamination 10 \
    --threads 40
```

Key output: `01_large_contigs/assembly_for_filtering.fa` (use in Step 5) and
`HQ_extracted/` (clean HQ genomes).

### Step 3 - Detect chimeric contigs (read-level taxonomy)

Streams the BAM once, cross-references every read against Kraken2 per-read
taxonomy, builds a per-contig taxonomic profile, and flags contigs whose reads
originate from more than one organism.

```bash
echidna check-read-chimeras \
    --bam mapping.bam \
    --assembly assembly.fa \
    --kraken2-output kraken2_reads.txt \
    --kraken2-report kraken2_report.txt \
    -o 02_chimera/ \
    --min-mapq 20 --window 10000 --window-step 5000
```

Outputs `chimera_read_report.tsv` (per-contig) and `chimera_read_windows.tsv`
(sliding-window detail for classify-clipping).

Classifier choice: `--kraken2-output` (k-mer, fastest) or `--kaiju-output`
(protein-level, more sensitive on divergent/novel organisms). Give exactly one.
For names use `--kraken2-report` with Kraken2, or `--kaiju-names names.dmp` with
Kaiju. BLAST is not used here on purpose: per-read BLAST against nt is orders of
magnitude slower and short reads give low-specificity hits.

### Step 3b - Detect chimeric MAGs (GUNC + per-contig taxonomy, optional)

Run after binning for a gene-consistency + composition check on complete bins.
The verdict uses exactly two signals: GUNC gene-level (`pass.GUNC == False` ->
chimeric) and the per-contig taxonomy foreign LENGTH fraction from
`classify-contigs` (`>= --min-foreign-frac` -> chimeric). Both chimeric ->
`Chimeric`, one -> `Check`, neither -> `Clean`. CheckM2 contamination, GTDB-Tk
taxonomy and a strong GC outlier are reported as columns only - they do not
change the verdict.

```bash
# the per-contig taxonomy signal: run-mmseqs-contigs -> classify-contigs
echidna run-mmseqs-contigs --bins-dir 06_bins/ --db /data/mmseqs_gtdb/gtdb \
    -o 02_chimera/contig_taxonomy.tsv --threads 40
echidna classify-contigs --bins-dir 06_bins/ \
    --mmseqs-output 02_chimera/contig_taxonomy.tsv \
    -o 02_chimera/contigs.classified.tsv                      # writes <out>.bins.tsv

echidna check-chimeras \
    --bins-dir 06_bins/ \
    --output 02_chimera/chimera_report.tsv \
    --gunc-dir 07_gunc/ \
    --contig-taxonomy 02_chimera/contigs.classified.tsv.bins.tsv \
    --gtdbtk-dir 09_taxonomy/ \
    --checkm2-report 07_checkm2/quality_report.tsv
```

The verdict needs at least one of `--gunc-dir` / `--contig-taxonomy`; with only
GUNC a failing bin is `Check`, not `Chimeric`.

Pre-dRep filter: exclude rows with `chimera_verdict == Chimeric` (GUNC and
per-contig taxonomy both chimeric); send `Check` to manual review (the
informational columns - gunc_pass, contig_foreign_frac, gc_strong_outlier - say
why).

### Step 4 - Classify clipping positions

Characterise each anvi'o-reported clipping position using BAM evidence so you can
prioritise which splits are biologically meaningful.

```bash
echidna classify-clipping \
    --bam mapping.bam \
    --misassemblies misassemblies/ \
    --taxonomy-windows 02_chimera/chimera_read_windows.tsv \
    --output 03_classified/clipping_classified.tsv
```

Classification labels:

| Label | Biological meaning |
|---|---|
| `end_artefact` | Near contig terminus; assembly edge noise |
| `repeat_collapse` | Coverage spike + low-entropy clipped bases; collapsed repeat |
| `deletion_artefact` | Local depth drop; internal deletion or coverage collapse |
| `chimera_candidate` | Discordant pairs + large inserts, possibly a taxonomy shift |
| `chimeric_join` | Clipped reads' supplementary alignment (`SA:Z:`) maps to one other contig: the two contigs are joined here (Trigodet et al. 2025 read-level signal) |
| `sv_candidate` | SV signal (discordant) without depth anomaly or taxon shift |
| `low_confidence` | Ambiguous or insufficient evidence |

### Step 5 - Filter misassemblies

```bash
# Conservative: use all anvi'o-reported positions (--safe, the default)
echidna filter-misassemblies \
    --assembly 01_large_contigs/assembly_for_filtering.fa \
    --misassemblies misassemblies/ \
    --hq-genomes 01_large_contigs/HQ_extracted/ \
    --output 04_clean/assembly_clean.fa \
    --mode split --preserve-hq-with-issues \
    --chimera-report 02_chimera/chimera_read_report.tsv \
    --split-hq-circular --safe
```

`--safe` uses every anvi'o-reported position; `--min-clip-coverage N` is the
selective alternative (only split where at least N reads are clipped).

### Step 6 - Binning (external)

```bash
minimap2 -ax map-hifi 04_clean/assembly_clean.fa reads.hifi.fastq.gz | samtools sort -o clean_mapping.bam
samtools index clean_mapping.bam
jgi_summarize_bam_contig_depths --outputDepth depth.txt clean_mapping.bam
metabat2 -i 04_clean/assembly_clean.fa -a depth.txt -o 05_bins/bin -t 40
```

### Step 7 - Quality assessment with CheckM2

```bash
conda run -n checkm2 echidna run-checkm2 --genome-dir 05_bins/ -o 06_checkm2/ --threads 40
```

### Step 8 - Organise genomes by quality

```bash
echidna organize-mags \
    --checkm2-root 06_checkm2/ --genomes-dir 05_bins/ \
    --output 07_filtered/ --tiers HQ,MQ --symlink
```

MiMAG quality thresholds (defaults):

| Tier | Completeness | Contamination |
|---|---|---|
| HQ | >= 90% | <= 5% |
| MQ | >= 70% | <= 10% |
| LQ | >= 50% | <= 10% |

The defaults split the MIMAG medium band into MQ/LQ. For a wider medium band use `--mq-comp-min 50 --lq-comp-min 0`.
These remain internal inclusive thresholds, not strict MIMAG boundary rules.
The optional RNA gate uses `--mimag-rna-dir` (see Option A).

### Step 9 - Taxonomic classification with GTDB-Tk

```bash
conda run -n gtdbtk echidna run-gtdbtk \
    --selected-dir 07_filtered/Selected/ \
    --output 08_taxonomy/ --cpus 40 --tiers HQ,MQ
```

## Understanding misassembly detection

Echidna uses soft-clipping information from read mapping to detect misassemblies. A
clipping position is a genomic location where reads are predominantly soft-clipped:
strong evidence that two unrelated sequences were joined during assembly.

`classify-clipping` adds mechanistic insight to each position using BAM signals:

| Signal | What it detects |
|---|---|
| depth_ratio | Coverage spike (repeat collapse) or drop (deletion) |
| discordant_fraction | Reads whose mates map to a different contig or wrong orientation |
| large_insert_fraction | Pairs with abnormally large insert sizes (> mean + 3 sd) |
| clipped_base_entropy | Repetitive vs diverse sequence at the break point |
| sa_partner_contig | Clipped reads whose supplementary alignment maps to one other contig: contig join (`chimeric_join`) |

When `check-read-chimeras` output is provided, taxonomy shifts near a clipping
position provide an additional, strong signal for chimera classification.

BAM alignment categories used throughout:
- Depth counting: primary + secondary (supplementary excluded, `-F 2048`)
- Pair statistics: primary only (FLAG `0x100 == 0`)
- Supplementary: excluded entirely


Current operational note: run anvi’o BEFORE `process-large-contigs`, since the latter consumes its clipping output. `classify-clipping` is report-only; quality organization does not apply a chimera exclusion policy. See [the maintained wiki](../wiki/Home.md).
