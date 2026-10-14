# Command reference

Every subcommand also has `-h/--help`. Flags below are the stable ones.

## process-large-contigs

```
Required:
  -a, --assembly PATH         Input assembly FASTA
  -m, --misassemblies PATH    Directory with *-clipping.txt files
  -o, --output-dir PATH       Output directory

Thresholds:
  --max-length INT            Large contig threshold (default: 300000)
  --min-completeness FLOAT    Min completeness for HQ (default: 50)
  --max-contamination FLOAT   Max contamination for HQ (default: 10)

CheckM2:
  --threads INT               Threads (default: 1)
  --skip-checkm2              Skip CheckM2 (requires --checkm2-results)
  --checkm2-results PATH      Existing CheckM2 output directory
```

## check-read-chimeras

```
Required:
  --bam PATH                  Sorted, indexed BAM file
  -o, --output-dir PATH       Output directory
  one classifier (mutually exclusive):
    --kraken2-output PATH     Kraken2 per-read classification output, or
    --kaiju-output PATH       Kaiju per-read output (protein-level, more
                              sensitive on divergent taxa)

Optional:
  --kraken2-report PATH       Kraken2 report (adds taxon names; with Kraken2)
  --kaiju-names PATH          NCBI names.dmp (adds taxon names; with Kaiju)
  --assembly PATH             Assembly FASTA (for contig lengths if not in BAM)
  --min-mapq INT              Min mapping quality (default: 20)
  --min-reads INT             Min reads per contig to report (default: 10)
  --window INT                Sliding window size in bp (default: 10000)
  --window-step INT           Window step in bp (default: 5000)
  --window-threshold INT      Min contig length for windowed analysis (default: 50000)
  --exclude-unclassified      Exclude taxid=0 reads from diversity calculations
  --dry-run
```

## check-chimeras

Emits a per-bin 3-state **verdict** from TWO signals only:

- **GUNC** gene-level: `pass.GUNC == False` means chimeric.
- **Per-contig taxonomy** (MMseqs via classify-contigs, `--contig-taxonomy`): the
  foreign LENGTH fraction >= `--min-foreign-frac` (default 0.05) means chimeric.

| GUNC | Per-contig taxonomy | Verdict |
|---|---|---|
| chimeric | chimeric | `Chimeric` |
| one of them chimeric | | `Check` (review) |
| neither | | `Clean` |

With a signal missing the verdict uses the one available (a single firing ->
`Check`); with neither available -> `Unknown`.

CheckM2 contamination, GTDB-Tk taxonomy/warning and GC are **reported as columns
only** and do NOT change the verdict. The GC column (`gc_strong_outlier`) names a
long contig whose GC deviates from the bin's length-weighted mean by more than
`--gc-max-dev`, ignoring contigs shorter than `--gc-min-length` (so a short
high-GC contig flags nothing).

All sources are joined by a canonical genome name (FASTA extension plus any
`--name-strip` suffix removed); bins with no match in a source are WARNED about
(so a silent N/A becomes visible).

```
Required:
  --bins-dir PATH             Directory containing bin FASTA files

Verdict signals:
  --gunc-dir PATH             GUNC output directory (pass.GUNC drives the verdict)
  --contig-taxonomy PATH      classify-contigs <output>.bins.tsv (per-contig taxonomy)
  --min-foreign-frac FLOAT    Foreign LENGTH fraction at/above which the per-contig
                              signal counts as chimeric (default: 0.05)

Informational only (reported columns, not in the verdict):
  --gtdbtk-dir PATH           GTDB-Tk output directory
  --checkm2-report PATH       CheckM2 quality_report.tsv
  --gc-min-length INT         Ignore contigs shorter than this for the GC column (5000)
  --gc-max-dev FLOAT          Report a GC outlier past this deviation from bin mean (0.10)
  --name-strip LIST           Suffixes to strip from genome names when joining sources
  --taxonomy-level STR        GTDB rank shown in gtdbtk_taxonomy (default: phylum)
  --extensions LIST           Genome extensions (default: fa,fasta,fna)
  -o, --output PATH / --dry-run
```

Output columns: `bin`, `n_contigs`, `gc_mean_pct`, `gc_strong_outlier`,
`checkm2_contamination`, `gtdbtk_taxonomy`, `gtdbtk_warning`, `gunc_pass`,
`gunc_css`, `contig_foreign_frac`, `chimera_verdict`, `reasons`.

The verdict needs at least one of `--gunc-dir` / `--contig-taxonomy`. Pre-dRep
filter: exclude `chimera_verdict == Chimeric`; send `Check` to manual review (the
informational columns say why).

## classify-contigs

Per-contig taxonomy check, independent of GUNC. Assigns each contig a taxon
(Kraken2 or Kaiju on the contigs, or MMseqs2 via `run-mmseqs-contigs`), resolves
it to `--rank`, and per bin counts the distinct rank-level taxa. Two or more
distinct taxa flags the bin as a taxonomic chimera. Like MAGpurify2's taxonomy
module, the bin's **majority** rank-taxon is treated as its identity and the
**minority (foreign) contigs are named** as the likely contaminants. Contigs
classified only above the rank are left unresolved (they cannot prove a split).
Report-only; orthogonal to the GUNC/GC signals in check-chimeras.

Two resolution modes, picked automatically:
- **lineage mode** (MMseqs2 run with `--tax-lineage`, which `run-mmseqs-contigs`
  does): the rank is read from the per-contig lineage string, so **no nodes.dmp
  is needed** and ranks stay in the DB's own system (GTDB when the DB is GTDB) -
  the recommended path, consistent with GTDB-Tk.
- **taxid mode** (Kraken2/Kaiju, or an MMseqs2 TSV with no lineage): each contig's
  taxid is walked up the NCBI tree with `--nodes-dmp` (required here).

Produce the input first, e.g. `kraken2 --db DB --output contigs.kraken2 all_contigs.fa`
or `echidna run-mmseqs-contigs ...` (the sequence-id column is the contig name).

```
Required:
  --bins-dir PATH             Directory containing bin FASTA files
  one classifier (mutually exclusive):
    --mmseqs-output PATH      MMseqs2 taxonomy TSV (run-mmseqs-contigs); GTDB-capable
    --kraken2-output PATH     Kraken2 per-contig output (nucleotide, DB-biased)
    --kaiju-output PATH       Kaiju per-contig output (protein, more sensitive)

Optional:
  --nodes-dmp PATH            NCBI nodes.dmp - required for Kraken2/Kaiju and for an
                              MMseqs2 TSV without a lineage column; not needed in
                              lineage mode
  -o, --output PATH           Per-contig TSV (default: contig_taxonomy.tsv);
                              a per-bin summary is written as <output>.bins.tsv
  --names-dmp PATH            NCBI names.dmp for readable names (taxid mode only)
  --rank STR                  Rank to compare contigs at (default: order; 'domain'
                              maps to NCBI superkingdom in taxid mode)
  --min-foreign-frac FLOAT    Min foreign LENGTH fraction (0-1) for chimera_flag=True
                              (default: 0.05). Below it the foreign contigs are still
                              listed and foreign_frac reported, just not flagged - a
                              few genes / one short contig is likely HGT.
  --extensions LIST           Genome extensions (default: fa,fasta,fna)
  --dry-run
```

Per-contig columns: `contig`, `bin`, `contig_taxid`, `rank_taxon`, `is_foreign`.
Per-bin columns: `bin`, `n_contigs`, `n_contigs_resolved`, `n_taxa`,
`majority_taxon`, `chimera_flag`, `foreign_frac`, `n_foreign_contigs`,
`foreign_contigs`, `taxa`. The majority taxon and `foreign_frac` are weighted by
contig LENGTH (bp); `foreign_frac` is what `check-chimeras --contig-taxonomy`
consumes.

Rank choice: `order` is a good default; `class` separates class-level
contaminants (e.g. an Alphaproteobacteria contig in a Gammaproteobacteria bin);
phylum is too coarse, genus/species gets noisy on short contigs. In taxid mode
Kraken2/Kaiju use NCBI taxonomy, so ranks will not line up one-to-one with
GTDB-Tk - prefer MMseqs2 + a GTDB DB (lineage mode) for GTDB consistency.

## run-mmseqs-contigs

Runs MMseqs2 taxonomy on the contigs of a set of bins and writes the per-contig
TSV that `classify-contigs --mmseqs-output` consumes. MMseqs2 is the per-contig
classifier MAGpurify2 uses: protein-level 2bLCA, sensitive on divergent/novel
taxa, and its DB can be GTDB-derived so ranks match GTDB-Tk. Heavy dependency, so
it lives in its own `mmseqs2` env: `conda run -n mmseqs2 echidna run-mmseqs-contigs ...`.

Needs a prebuilt MMseqs2 taxonomy DB (`--db`), e.g.
`conda run -n mmseqs2 mmseqs databases GTDB /data/mmseqs_gtdb tmp`. `--db` is the
MMseqs2 DB **prefix** (the handle, e.g. `/data/mmseqs_gtdb/gtdb`), not the folder
or a single `.dbtype`/`.index` file.

```
Required:
  --bins-dir PATH             Directory with bin FASTA files (contigs concatenated)
  --db PATH                   MMseqs2 taxonomy sequence DB (mmseqs databases/createtaxdb)
  -o, --output PATH           Output per-contig taxonomy TSV

Optional:
  -T, --threads INT           Threads (default: 1)
  --tmp-dir PATH              Working dir for MMseqs2 DBs/tmp (default: <output>.mmseqs_work)
  --extensions LIST           Genome extensions (default: fa,fasta,fna)
  --mmseqs-exe NAME           MMseqs2 executable (default: mmseqs)
  --keep-tmp                  Keep the working dir (default: removed)
  --mmseqs-args STR           Extra args passed to `mmseqs taxonomy` (quoted)
  --dry-run                   Print the MMseqs2 commands without running
```

Runs `mmseqs createdb` -> `mmseqs taxonomy --tax-lineage 1` -> `mmseqs createtsv`.
The `--tax-lineage 1` is added automatically so the output carries a lineage
column, letting `classify-contigs` resolve any rank without a nodes.dmp. Note:
contig headers must be unique across bins (they are concatenated into one query
FASTA).

## check-ssu-taxonomy

Classifies each MAG's 16S rRNA genes (blastn vs an SSU reference DB)
and checks they agree with the genome-level GTDB-Tk taxonomy. A 16S resolving to a
different taxon than the genome flags misbinning or a contaminating rRNA -
orthogonal to the GUNC/GC/MMseqs per-contig signals. Report-only.

Only recovered ~full-length 16S records are judged (rRNA is often
fragmented/unbinned), and it is confirmatory, not exhaustive. The 16S records come from
`run-mimag-rna --outseq` (barrnap writes `<genome>.rrna.fna` into
`<output>/rrna_fasta/`). blastn is in the echidna env; barrnap runs in the mimag
env - so extraction and classification stay in their own envs.

Build the blastn DB first (GTDB `ssu_all.fna` recommended, 1:1 with GTDB-Tk):
`makeblastdb -in ssu_all.fna -dbtype nucl -out /data/ssu_gtdb/ssu`.

```
Required:
  --ssu-fasta-dir PATH        Directory of <genome>.rrna.fna (run-mimag-rna --outseq)
  --ssu-db PREFIX             blastn SSU DB prefix (makeblastdb -dbtype nucl)
  --gtdbtk-dir PATH           GTDB-Tk output (genome-level taxonomy to compare)

Optional:
  --ssu-taxonomy {gtdb,silva} How to read the DB hit lineage (default: gtdb).
                              gtdb = d__/p__/c__ prefixes, compares 1:1 with
                              GTDB-Tk; silva = positional names, agreement is
                              best-effort (SILVA names differ from GTDB)
  --rank STR                  Rank to compare 16S vs genome (default: order)
  --min-length INT            Min 16S length to classify (default: 1000 bp)
  -o, --output PATH           Output TSV (default: ssu_taxonomy.tsv)
  --threads INT / --task STR  blastn threads / task (default: 1 / megablast)
  --blastn-exe NAME           blastn executable (default: blastn)
  --keep-tmp                  Keep the combined query FASTA and blast TSV
```

Output columns: `genome`, `ssu_id`, `ssu_type`, `has_16s`, `ssu_len`,
`ssu_pident`, `ssu_taxon`, `genome_taxon`, `agree` (True/False/NA), `note`
(`no_16s`, `partial_16s`, `no_blast_hit`, `low_identity`). Each recovered 16S
gets its own row; 18S/23S/5S records are ignored. `agree=False` is the flag to
review.

DB choice: GTDB SSU keeps names aligned with GTDB-Tk so `agree` is meaningful;
SILVA is broader (better 16S coverage) but its names differ (e.g. Proteobacteria
vs GTDB Pseudomonadota), so at phylum the comparison is noisy - use class/order,
and read `ssu_taxon`/`genome_taxon` directly rather than trusting `agree`.

## classify-clipping

```
Required:
  --bam PATH                  Sorted, indexed BAM file
  -m, --misassemblies PATH    Directory with *-clipping.txt files
  -o, --output PATH           Output classification TSV

Optional:
  --min-mapq INT              Min mapping quality (default: 20)
  --window INT                Read window around each position in bp (default: 500)
  --min-clip-coverage N       Only classify positions with >=N clipped reads
  --taxonomy-windows TSV      chimera_read_windows.tsv from check-read-chimeras
  --taxonomy-flank BP         Extend shift windows by this many bp (default: 5000)
  --assembly FASTA            Assembly FASTA - required for --self-blast
  --self-blast                Self-BLAST repeat detection (needs --assembly + BLAST+)
  --dry-run
```

The `chimeric_join` label: when the soft-clipped parts of reads at a position map
(via their `SA:Z:` supplementary-alignment tag) consistently to one other contig,
the two contigs are joined there. This read-level join signal (Trigodet et al.
2025) is computed automatically from the BAM; see output columns
`sa_partner_contig` / `sa_partner_fraction`.

## filter-misassemblies

```
Required:
  -a, --assembly PATH         Input assembly FASTA
  -m, --misassemblies PATH    Directory with *-clipping.txt files

Optional:
  -g, --hq-genomes PATH       Directory with HQ genome files
  -o, --output PATH           Output filtered assembly
  -l, --min-length INT        Minimum fragment length (default: 2500)
  --mode {remove,split}       Misassembly handling (default: split)
  --preserve-hq-with-issues   Split HQ contigs with clipping instead of removing

Clipping selection (mutually exclusive):
  --safe                      Use ALL anvi'o-reported positions (default)
  --min-clip-coverage N       Only split where >=N reads are clipped

Chimera awareness:
  --chimera-report PATH       TSV from check-chimeras or check-read-chimeras
  --chimera-risk-threshold    Min risk level to act on (default: Medium)

HQ circular handling:
  --split-hq-circular         Force-split HQ circular contigs with clipping
  --hq-circular-min-length    Length threshold for circular detection (default: 200000)

Other:
  --dry-run
  --stats PATH                Save statistics to TSV
```

## remove-hq-contigs

```
Required:
  -a, --assembly PATH         Input assembly FASTA
  --hq-dir PATH               Directory with HQ genome files
    OR
  --hq-list PATH              Text file with HQ contig IDs (one per line)

Optional:
  -o, --output PATH           Output filtered assembly
  -l, --min-length INT        Minimum contig length (default: 0)
  --dry-run
  --stats PATH                Save statistics to TSV
```

## organize-mags  (formerly filter-checkm2)

Splits genomes into `Selected/HQ`, `Selected/MQ`, `Selected/LQ` by MiMAG tier, and
can rename each output `<sample>__<taxon>__<original>` using the sample parsed from
the filename and the GTDB-Tk taxonomy. `filter-checkm2` still works as a deprecated
alias.

```
Required:
  --checkm2-root PATH         Root directory with CheckM2 results
  --genomes-dir PATH          Directory with genome files
  -o, --output PATH           Output directory

Quality thresholds:
  --hq-comp-min FLOAT         HQ min completeness (default: 90)
  --hq-cont-max FLOAT         HQ max contamination (default: 5)
  --mq-comp-min FLOAT         MQ min completeness (default: 70)
  --mq-cont-max FLOAT         MQ max contamination (default: 10)
  --lq-comp-min FLOAT         LQ min completeness (default: 50)
  --lq-cont-max FLOAT         LQ max contamination (default: 10)

Renaming (sample + taxonomy):
  --rename                    Rename outputs as <sample>__<taxon>__<original>
  --gtdbtk-dir PATH           GTDB-Tk summaries, for the <taxon> part of the name
  --sample-regex RE           Capture group for the sample id
                              (default: ^([^._]+) = token before first . or _)
  --tax-level LEVEL           Rank used in name: domain..species (default: genus;
                              falls back up ranks, NoTax if no GTDB hit)

MIMAG rRNA/tRNA HQ gate:
  --mimag-rna-dir PATH        run-mimag-rna output (mimag_rna_summary.tsv);
                              demotes HQ genomes lacking 5S/16S/23S rRNA or
                              >=18 tRNA amino acids to MQ

Other:
  --tiers LIST                Tiers to select (default: HQ,MQ,LQ)
  --extensions LIST           Genome extensions (default: fa,fasta,fna,...)
  --symlink | --copy          Link or copy files (default: symlink)
  --name-map / --strip-prefix / --strip-suffix / --name-prefix / --name-suffix
  --dry-run
```

MiMAG note: the defaults split the MIMAG medium band into MQ/LQ. For strict MIMAG
tiers (HQ comp>90/cont<5, MQ comp>=50/cont<10, LQ comp<50) run with
`--mq-comp-min 50 --mq-cont-max 10 --lq-comp-min 0 --lq-cont-max 10`. Full MIMAG HQ
also needs rRNA/tRNA (`--mimag-rna-dir`).

## run-checkm2

```
Required (one of):
  --selected-dir PATH         Directory with HQ/, MQ/, LQ/ subdirectories
  --genome-dir PATH           Single directory with genome files

Required:
  -o, --output PATH           Output directory

Optional:
  --extension STR             Genome file extension (default: fa)
  --threads INT               Threads for CheckM2 (default: 1)
  --database PATH             CheckM2 diamond database (.dmnd FILE)
  --force                     Overwrite existing CheckM2 output
  --tiers LIST                Tiers to process (default: HQ,MQ,LQ)
  --checkm2-args STR          Extra arguments passed to checkm2
```

## run-gtdbtk

```
Required (one of):
  --selected-dir PATH         Directory with HQ/, MQ/, LQ/ subdirectories
  --genome-dir PATH           Single directory with genome files

Required:
  -o, --output PATH           Output directory

Optional:
  --extension STR             Genome file extension (default: fa)
  --cpus INT                  CPUs for GTDB-Tk (default: 1)
  --place-species             Pass --place_species (GTDB-Tk >=2.7): place in the
                              pplacer tree even when skani classifies
  --skip-existing             Skip classify_wf if a gtdbtk.*.summary.tsv already
                              exists in the output (reuse previous results)
  --mash-db PATH              DEPRECATED/ignored - GTDB-Tk removed Mash in v2.5.0
                              (now uses skani)
  --tiers LIST                Tiers to process (default: HQ,MQ,LQ)
  --gtdbtk-args STR           Extra arguments passed to gtdbtk
```

GTDB-Tk version note (latest 2.7.2): since v2.5.0 GTDB-Tk screens with skani and no
longer uses Mash, so `--mash_db` is gone. The old `--skip_ani_screen` was replaced
by `--place_species` in v2.7.0. Echidna parses both the old (`fastani_*`) and new
(`closest_genome_*`) summary column names, so reports from any recent GTDB-Tk
version work.

## run-mimag-rna

Checks the MIMAG rRNA (5S/16S/23S) + tRNA (>=18 amino acids) high-quality criterion
that CheckM2 cannot evaluate, via barrnap + tRNAscan-SE. Writes
`mimag_rna_summary.tsv`; feed it to `organize-mags --mimag-rna-dir` to gate HQ.

```
Required (one of):
  --selected-dir PATH         Directory with HQ/, MQ/, LQ/ subdirectories
  --genome-dir PATH           Single directory with genome files

Required:
  -o, --output PATH           Output directory

Optional:
  --extension STR             Genome file extension (default: fa)
  --kingdom {bac,arc,euk}     Kingdom for barrnap/tRNAscan-SE (default: bac;
                              with --split-by-domain, the fallback only)
  --split-by-domain           Pick bac/arc per genome from GTDB-Tk (needs
                              --gtdbtk-dir); handles mixed bacteria+archaea sets
  --gtdbtk-dir PATH           GTDB-Tk output dir for --split-by-domain
  --threads INT               Threads for barrnap (default: 1)
  --tiers LIST                Tiers to process (default: HQ,MQ,LQ)
  --keep-intermediate         Keep per-genome barrnap/tRNAscan output files
  --outseq                    Also write detected rRNA genes to
                              <output>/rrna_fasta/<genome>.rrna.fna; the 16S
                              feeds check-ssu-taxonomy
```

Mixed bacteria + archaea: run once with `--split-by-domain --gtdbtk-dir <taxonomy>`.
GTDB-Tk writes one summary per domain (`gtdbtk.bac120.summary.tsv`,
`gtdbtk.ar53.summary.tsv`), so each genome gets barrnap/tRNAscan with the correct
kingdom automatically; genomes absent from GTDB-Tk fall back to `--kingdom`.

## run-gunc

Runs GUNC gene-level chimerism/contamination detection, complementing
`check-chimeras`. Feed the output back with `check-chimeras --gunc-dir`. The
command runs an env-health check and refuses to launch a GUNC broken by an
incompatible pandas/numpy/python (GUNC is old - see installation.md).

```
Required (one of):
  --selected-dir PATH         Directory with HQ/, MQ/, LQ/ subdirectories
  --genome-dir PATH           Single directory with genome files

Required:
  -o, --output PATH           Output directory

Optional:
  --file-suffix STR           Genome file suffix (default: .fa)
  --threads INT               Threads (default: 1)
  --db-file PATH              GUNC diamond database (.dmnd); else uses $GUNC_DB
  --tiers LIST                Tiers to process (default: HQ,MQ,LQ)
  --contig-taxonomy-output    Also write per-contig taxonomy (<genome>.contig_assignments.tsv);
                              check-chimeras consumes it for the multi-clade signal
  --gunc-args STR             Extra arguments passed to gunc run
```

## setup-workflow

Configures the Snakemake workflow: points `workflow/envs/*.yaml` at this checkout,
writes paths into `workflow/config.yaml`, validates database kinds, and (when env
names are given) sanity-checks them. `--install` builds every per-rule conda env.

```
Optional:
  --repo PATH                 Echidna checkout (default: auto-detected)
  --workflow-dir PATH         Workflow dir (default: <repo>/workflow)
  --genomes-dir PATH          Input genomes directory (config)
  --genome-ext STR            Genome extension (config)
  --outdir PATH / --threads INT
  --gtdbtk-data DIR           GTDB-Tk data DIRECTORY
  --gunc-db FILE              GUNC .dmnd FILE
  --checkm2-db FILE           CheckM2 .dmnd FILE
  --sample-regex RE / --tax-level LEVEL
  --echidna-env / --checkm2-env / --gtdbtk-env / --gunc-env / --mimag-env NAME
                              Use these existing envs (sets use_named_envs)
  --install                   Build per-rule conda envs (snakemake)
  --no-check                  Skip the env sanity check
  --dry-run
```

## fetch-databases

Registers databases you already have (validated, written to config) or downloads
missing ones with each tool's own downloader. Path KIND is enforced: GTDB-Tk =
directory, GUNC and CheckM2 = `.dmnd` file.

```
Optional:
  --config PATH               config.yaml to update (default: workflow/config.yaml)
  --db-dir PATH               Where to download (default: ./databases)
  --tools LIST                Which DBs (default: gtdbtk,gunc,checkm2)
  --download                  Download any DB not supplied as an existing path
  --gtdbtk-data DIR           Existing GTDB-Tk data DIRECTORY
  --gunc-db FILE              Existing GUNC .dmnd FILE
  --gunc-db-source {progenomes_2.1,progenomes_3,gtdb_95,gtdb_214}
                              Which GUNC DB to download (GUNC >=1.1.1; default:
                              progenomes_2.1; use gtdb_214 for environmental MAGs)
  --checkm2-db FILE           Existing CheckM2 .dmnd FILE
  --gtdbtk-env / --gunc-env / --checkm2-env NAME   Env to run each downloader in
  --dry-run
```

## check-envs

Sanity-checks existing conda envs: reports each tool's version, whether Echidna is
importable there, and (for GUNC) pandas/numpy/python compatibility. Prints the
exact install/fix command for anything that fails.

```
Optional (at least one):
  --echidna-env NAME           Env with the light CLI tools + Echidna
  --checkm2-env NAME          Env with CheckM2
  --gtdbtk-env NAME           Env with GTDB-Tk
  --gunc-env NAME             Env with GUNC
  --mimag-env NAME            Env with barrnap + tRNAscan-SE
  --repo PATH                 Echidna checkout used in fix hints
```

## run-anvi-misassemblies

Wrapper around anvi'o `anvi-script-find-misassemblies` (soft-clipping detection,
PacBio HiFi only). Run inside the anvi'o env. `--just-do-it` is on by default
(anvi'o refuses to proceed without it non-interactively).

```
Required:
  -b, --bam PATH              Sorted, indexed HiFi BAM mapped to the assembly
  -o, --output-prefix PATH    Writes <prefix>-clipping.txt (parent dir = misassemblies dir)

Optional:
  -T, --threads INT           Threads (default: 1)
  --clipping-ratio FLOAT      anvi clipping ratio (default: 0.8)
  --min-dist-to-end INT       anvi min distance to contig end (default: 100)
  --no-just-do-it             Do not pass --just-do-it to anvi'o
  --anvi-args STR             Extra args passed to anvi-script-find-misassemblies
```

## misassembly-wf

Runs the whole misassembly-cleaning branch (PacBio HiFi only), each step in the
right env. Steps: `anvi` -> `process` -> `readchim` (only if Kraken2/Kaiju
taxonomy is given) -> `classify` -> `filter`.

```
Required:
  -a, --assembly PATH         Assembly FASTA (reads were mapped to this)
  -b, --bam PATH              Sorted, indexed HiFi BAM
  -o, --output-dir PATH       Base output directory

Optional:
  --threads INT               Threads (default: 1)
  --anvio-env NAME            Env with anvi'o + Echidna (default: anvio-9)
  --checkm2-env NAME          Env with CheckM2 (default: checkm2)
  --anvi-prefix STR           Prefix for clipping files (default: MisAsm)
  --max-length / --min-completeness / --max-contamination   process-large-contigs tuning
  --kraken2-output PATH  --kraken2-report PATH    enable check-read-chimeras (Kraken2)
  --kaiju-output PATH    --kaiju-names PATH        enable check-read-chimeras (Kaiju)
  --from STEP                 Resume at a step
  --skip STEP[,STEP]          Skip steps
  --dry-run                   Print commands without running
```

## mag-quality-wf

Runs the whole MAG quality + taxonomy branch from a bin set, each step in the
right env. Steps: `checkm2` -> `organize1` -> `gtdbtk` -> `gunc` -> `mimag` ->
`chimeras` -> `organize2`. Genome names are frozen after `organize1` so GUNC,
check-chimeras and run-mimag-rna match; GUNC is run with
`--contig-taxonomy-output`.

```
Required:
  --bins-dir PATH             Directory of bin FASTA files (from binning)
  -o, --output-dir PATH       Base output directory

Optional:
  --threads INT               Threads (default: 1)
  --tiers LIST                Tiers to carry downstream (default: HQ,MQ)
  --genome-ext STR            Bin extension without dot (default: fa)
  --echidna-env NAME           Env for pure-Echidna steps (default: current env)
  --checkm2-env / --gtdbtk-env / --gunc-env / --mimag-env   tool envs
  --gtdbtk-data PATH          GTDB-Tk data dir (sets GTDBTK_DATA_PATH)
  --gunc-db PATH              GUNC .dmnd (else $GUNC_DB)
  --checkm2-db PATH           CheckM2 .dmnd (else CheckM2 default)
  --split-by-domain           run-mimag-rna picks bac/arc per genome from GTDB-Tk
  --skip-mimag                Skip the rRNA/tRNA MIMAG HQ gate
  --from STEP                 Resume at a step
  --skip STEP[,STEP]          Skip steps
  --dry-run                   Print commands without running
```

## final-report

Joins the per-MAG signals the other commands produce into one table - one row
per MAG. Reads only existing TSVs; runs nothing, needs no external tools. Every
input is optional and joined by genome name (FASTA extension stripped); missing
fields are `NA`.

```
Optional (pass whatever you have):
  --checkm2-report PATH       A single CheckM2 quality_report.tsv
  --checkm2-root PATH         Directory to search for quality_report.tsv (recursive)
  --gtdbtk-dir PATH           GTDB-Tk output directory
  --mimag-rna-dir PATH        run-mimag-rna output (mimag_rna_summary.tsv)
  --chimera-report PATH       check-chimeras chimera_report.tsv
  --contig-taxonomy PATH      classify-contigs <output>.bins.tsv
  --ssu-report PATH           check-ssu-taxonomy ssu_taxonomy.tsv
  --genomes-dir PATH          Define the genome set from a directory (else the
                              union of names across the inputs)
  --extensions LIST           Genome extensions for --genomes-dir (default: fa,fasta,fna)
  --tax-level STR             GTDB rank shown in rank_taxon (default: phylum)
  -o, --output PATH           Output TSV (default: final_report.tsv)
```

Columns: `genome`, `completeness`, `contamination`, `tier` (CheckM2), `taxonomy`,
`rank_taxon` (GTDB-Tk), `n_16S`, `rrna_complete`, `trna_complete`, `mimag_rna_ok`
(MIMAG), `chimera_verdict` (check-chimeras), `contig_chimera`, `n_foreign_contigs`
(classify-contigs), `ssu_agree`, `ssu_taxon` (check-ssu-taxonomy),
`gunc_mmseqs_concordance` (True/False/NA: do GUNC pass.GUNC and MMseqs
classify-contigs agree on chimeric-or-not), `concerns`.

`concerns` collects the red flags - `contam>10`, `chimera_High`, `chimera_Medium`,
`contig_chimera`, `ssu_disagree`, `mimag_rna_incomplete` - so you can sort/filter
on it. An empty `concerns` means no red flag fired; it is a summary, not a formal
pass (apply your own thresholds). The run prints a tally of each concern.

```
echidna final-report \
  --checkm2-report checkm2/quality_report.tsv \
  --gtdbtk-dir taxonomy/ --mimag-rna-dir mimag/ \
  --chimera-report chimera_report.tsv \
  --contig-taxonomy contig_taxonomy.tsv.bins.tsv \
  --ssu-report ssu_taxonomy.tsv \
  -o final_report.tsv
```
