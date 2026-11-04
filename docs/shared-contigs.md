# Shared contigs between MAGs

Use `check-shared-contigs` to check whether genome FASTAs contain the same
full-length contig sequence. It also reports reused contig IDs whose sequences
differ. The command uses the Python standard library and needs no database.

```bash
echidna check-shared-contigs --bins-dir bins/ --output-dir overlap/
```

By default this is a standalone report. Optional resolution flags export a new
MAG collection; original FASTAs are always read-only. The command does not
automatically change any chimera verdict. An overlap alone does
not establish contamination: repeated elements can be identical across genomes.
Review the evidence before changing bin membership.

## Inputs

Provide one FASTA file per MAG. By default `.fa`, `.fna`, `.fasta` and their
`.gz`/`.gzip` forms are recognized, case-insensitively. Other extensions can be
selected with `--extensions fa,fna,fasta`.

For a directory containing quality tiers:

```bash
echidna check-shared-contigs --bins-dir results/Selected/ \
  --recursive --output-dir overlap_selected/
```

For MAGs from multiple binners or collections:

```bash
echidna check-shared-contigs --bins-dir metabat_bins/ semibin_bins/ \
  --output-dir overlap_binners/
```

Single-directory bin labels are their relative file paths, including the
extension, such as `HQ/bin.1.fa`. Multiple input directories receive prefixes
`set1/`, `set2/`, etc., in the supplied order, to keep matching filenames
separate. `bin_summary.tsv` records each label's source path. Different symlink
paths are compared as different input MAGs; including the same path twice via
repeated or overlapping roots is rejected.

Use a new or empty output directory. Inputs are fully parsed before reports
are written. Empty FASTAs, empty sequences, duplicate IDs inside one MAG,
sequence data before a header, or non-IUPAC-DNA characters are errors rather
than a clean result.

## What counts as shared?

The complete sequences must match exactly after removing sequence whitespace,
uppercasing and accounting for reverse complement. SHA256 hashes index the
canonical sequences; the program does not retain all sequences in memory.
FASTA descriptions are ignored, and IDs use the exact first header token.
There is no extraction of numbers or truncation of assembler IDs.

| Input relationship | Report |
|---|---|
| Same ID, same complete sequence | Shared sequence |
| Different IDs, same complete sequence | Shared sequence |
| Reverse complements | Shared sequence |
| Same ID, different sequence | ID conflict; not a shared sequence on that basis |
| Same sequence repeated only within one MAG | Not shared between MAGs |
| Partial overlap, small sequence differences, or a changed circular start | Not matched |

Matching identifies identical full-contig sequences. It does not prove that
they originated from the same assembly contig. Reused IDs across different
samples are especially prone to conflicts. Ambiguous IUPAC bases are compared
literally; `N` is not a wildcard. No alignment, ANI calculation, fragment matching
or genome dereplication is performed.

## Reports

| File | Contents |
|---|---|
| `shared_contigs.tsv` | One row per occurrence of a sequence present in at least two MAGs: hash, number of MAGs, bin label, original contig ID and length |
| `id_conflicts.tsv` | One row per occurrence of an ID associated with different canonical sequences across MAGs |
| `bin_pairs.tsv` | Only pairs sharing a complete sequence: number of distinct shared sequences, shared bases, overlap fractions and sequence-set Jaccard |
| `bin_summary.tsv` | Every MAG, including those without overlap: source path, contig count, total length, shared contigs/bases/fraction and ID-conflict count |
| `summary.json` | Global counts and matching/denominator definitions |

When no overlap is found, `shared_contigs.tsv` and `bin_pairs.tsv` contain only
their headers; all MAGs still appear in `bin_summary.tsv`. Conflict tables are
also header-only when there are no conflicting IDs. A successful report returns
exit status 0 whether or not sharing was detected. Invalid inputs return a
nonzero status. With resolution enabled, the root reports describe the input
collection and `after/` describes the exported retained collection.

Fractions range from 0 to 1:

- In `bin_summary.tsv`, `shared_fraction_bp` is the fraction of all contig-record
  bases in that MAG belonging to sequences also found in another MAG. Each record
  is counted once even when the sequence is present in several other MAGs.
- In `bin_pairs.tsv`, `shared_unique_bp` counts each shared sequence once.
  `shared_fraction_a` and `shared_fraction_b` divide it by the total length of
  distinct sequences in the corresponding MAG (`unique_bp_a`, `unique_bp_b`).
- `sequence_jaccard` is the number of distinct shared sequences divided by the
  size of the union of the two sequence sets; it is not base-weighted or ANI.

Distinct IDs with an identical sequence within a MAG contribute separately to
its contig-record totals, but only once to pairwise sequence-set statistics.
An ID-conflict group can coexist with a shared-sequence group: for example, if
bins A and B have identical `ctg1` and bin C has a different `ctg1` sequence.

## Difference from bin refinement

[Binette](https://binette.readthedocs.io/en/stable/) performs bin refinement:
it constructs candidate bins using intersections, differences and unions of
contig sets and evaluates candidates with CheckM2. It is a separate tool and
is not required by `check-shared-contigs`.

Echidna does not construct hybrid MAGs by joining different input genomes.
It preserves a whole MAG, removes shared sequence copies, or selects a whole
alternative from the same dereplication cluster. A high completeness estimate
alone is not evidence of taxonomic coherence; the optional taxonomy gate below
tests reused contig classifications before accepting a genome.


## Resolve overlaps by quality

Two mutually exclusive flags enable changes in a new output collection:

```bash
# Preserve priority MAGs; exclude entire conflicting MAGs.
echidna check-shared-contigs --bins-dir bins/ -o resolved_mags/ \
  --quality-report checkm2/quality_report.tsv --exclude-mags

# Preserve shared contigs in priority MAGs; remove those copies from other MAGs.
echidna check-shared-contigs --bins-dir bins/ -o resolved_contigs/ \
  --quality-report checkm2/quality_report.tsv --remove-shared-contigs
```

The default priority is higher completeness, then lower contamination, then the
exact bin label for a deterministic tie. `--rank-by quality-score` instead
prioritizes `completeness - 5 * contamination`, with completeness, contamination
and label as tie-breakers. This is an Echidna priority rule, not dRep's full
representative score. No additional quality thresholds are applied.

Every MAG is checked against all higher-priority MAGs already retained. A chain
of overlaps is not treated as one inseparable cluster: A and C may both survive
when only A–B and B–C overlap. Selection is deterministic and greedy; it does
not maximize the number of MAGs or a global quality objective.

In contig-removal mode all records matching an already assigned sequence,
including renamed or reverse-complement copies, are removed from the
lower-priority MAG. A MAG with no remaining contigs is excluded rather than
exported as an empty FASTA. ID conflicts alone never trigger sequence removal.
Shared repeated elements may be biologically valid: these flags enforce an
operational non-overlapping collection, not proof of the contigs' true source.

### Quality table

`--quality-report` is required for either resolution flag. Accepts native
CheckM2 `Name`, `Completeness`, `Contamination` TSVs, or CSV/TSV tables with
`genome` (also `name` or `bin`), `completeness` and `contamination`. Headers are
case-insensitive. Completeness must be finite and between 0 and 100;
contamination must be finite and non-negative. Missing or duplicate matched
quality entries are errors. FASTA length is not used as a completeness estimate.

Genome identifiers can be exact bin labels, unique filenames or unique stems.
Ambiguous filenames across input directories require exact labels. Original
identifiers should be preserved when using dRep metadata; taxonomic prefixes
or other renaming are not guessed or removed.

### Replace excluded or trimmed representatives from a cluster

When original pre-dereplication FASTAs are still available:

```bash
echidna check-shared-contigs --bins-dir drep/dereplicated_genomes/ \
  --quality-report checkm2_all_candidates/quality_report.tsv \
  --exclude-mags --drep-dir drep/ \
  --candidates-dir all_original_mags/ -o resolved_with_replacements/
```

`--drep-dir` reads `data_tables/Cdb.csv` and its `secondary_cluster` column.
It does not use primary clusters, rerun dRep or recompute ANI. The input
`--bins-dir` is the already selected representative collection: at most one
input MAG per declared cluster is allowed. `Wdb.csv` is not needed because the
selected FASTAs identify the starting representatives.

For other dereplication tools, use `--cluster-table clusters.tsv`, containing
`genome` and `cluster` columns (CSV or TSV; `secondary_cluster` and `cluster_id`
are also accepted). Use the intended final cluster level. The table and the
candidate FASTA directories must be supplied together; cluster tables do not
contain genome sequences. Add `--recursive` if the supplied directories are
nested.

Available alternatives from clusters represented in the starting collection
must have their own quality entries. Replacements receive `candidate:` prefixes
in report identifiers; exact labels can resolve ambiguous metadata names.
Copies of the selected MAG in the candidate pool are recognized by filename
and canonical sequence multiplicities and are not counted as alternatives.
Unavailable members listed in the cluster table are reported separately; the
best available supplied candidate is not necessarily the best missing member.

Replacements are evaluated after the initial collection has been resolved:

1. For each affected cluster, consider the original and available alternatives
   in the selected quality order.
2. Check each candidate's full sequences against **every other retained MAG**,
   including already accepted replacements and provisionally trimmed MAGs.
3. Accept the first compatible whole MAG from that same cluster. In contig mode,
   it replaces the provisional trimmed genome, preserving the candidate whole.
4. Retry unresolved clusters after other replacements have freed contigs. The
   original representative can be restored if it becomes compatible.
5. In contig mode, recover original contigs that have become unclaimed after
   replacements, in priority order. This can also recover a nonempty remnant of
   a MAG that initially lost all its contigs.
6. Read and check the actual exported FASTAs again before publishing the output
   directory. The final `after/shared_contigs.tsv` must be header-only.

An accepted replacement is not repeatedly optimized afterward; this remains a
greedy selection. If no compatible alternative exists, exclusion mode leaves
the affected cluster without a representative, while contig mode keeps its
trimmed MAG when nonempty. These outcomes appear in `unresolved_clusters.tsv`.

### Resolved outputs and follow-up quality checks

| Output | Meaning |
|---|---|
| `retained/` | Final collection, plain `.fa` files with complete original headers |
| `excluded_mags/` | Full originals excluded or superseded by cluster replacements |
| `removed_contigs/` | Records removed from MAGs that remain as trimmed MAGs |
| `resolution.tsv` | Source/output mapping, actions, input quality, final quality status and blockers |
| `replacement_attempts.tsv` | Every tested original/alternative, retry round, result and blocking MAGs |
| `unresolved_clusters.tsv` | Affected clusters without a compatible whole replacement |
| `unavailable_cluster_members.tsv` | Relevant cluster members lacking a supplied FASTA |
| `after/` | Recomputed overlap reports for the final retained FASTAs |

Unique genome basenames are preserved with a `.fa` suffix; filename collisions
receive deterministic prefixes recorded in `resolution.tsv`. Sequence case and
line wrapping are normalized. Input FASTAs are never overwritten. Output is
prepared in a temporary directory and only exposed after the final overlap
check passes; data or metadata validation errors do not publish a partial set.

**Removing contigs invalidates that MAG's prior quality estimates.** A trimmed
MAG has `quality_status=requires_reassessment` and blank `final_completeness`
and `final_contamination`; input values remain explicitly labeled as input
values. For retained whole MAGs, final quality fields are the supplied estimates,
not a new prediction. Reassess the final collection, for example:

```bash
conda run -n checkm2 echidna run-checkm2 \
  --genome-dir resolved_contigs/retained/ \
  --output resolved_contigs/checkm2_reassessed/ --extension fa --threads 20
```

Trimmed genomes and replacements also need the applicable downstream
RNA/taxonomy/chimera checks. Do not transfer the removed representative's
reports onto a different genome. `requires_downstream_qc` records these cases.
A passed final overlap check is not a certification of biological quality.

Native dRep table definitions: [dRep advanced use](https://drep.readthedocs.io/en/latest/advanced_use.html).

## Reuse contig taxonomy when selecting replacements

Add `--contig-taxonomy` to gate **all retained MAGs**, including initial
representatives, trimmed genomes and whole cluster alternatives. It consumes
the **per-contig** output of `classify-contigs`, not its `.bins.tsv` summary.
Taxonomic membership and minority fractions are recomputed for each actual
candidate; the old representative's bin summary is never transferred.

```bash
echidna check-shared-contigs --bins-dir drep/dereplicated_genomes/ \
  --exclude-mags --quality-report checkm2_all/quality_report.tsv \
  --drep-dir drep/ --candidates-dir all_original_mags/ \
  --contig-taxonomy all_contigs.classified.tsv --taxonomy-rank order \
  --min-taxonomy-resolved-frac 0.5 --min-foreign-frac 0.05 \
  -o resolved_checked/
```

The classifier output must cover the original candidate pool, not just the
previously selected representatives. A candidate is tested against every other
retained MAG for full-contig sharing, then checked for internal taxonomic
coherence. If either check fails, the next candidate is tried. No taxonomy-based
union or transfer of contigs between genomes is performed.

Acceptance requires:

- At least `--min-taxonomy-resolved-frac` of **all MAG bases** must resolve at
  the table's rank (default `0.5`; configurable up to `1`). Missing,
  unclassified, above-rank and conflicting assignments do not count as resolved.
- Minority taxa must represent **less than** `--min-foreign-frac` of the
  **resolved bases** (default `0.05`, matching `check-chimeras`). With a threshold
  of zero, any nonzero minority is rejected; a single taxon remains acceptable.

These are explicit operational thresholds, not a claim of biological purity.
Coherence at order does not establish coherence at genus/species. Use the
desired rank when generating the table; `--taxonomy-rank` can assert it but
cannot infer a finer rank from an already aggregated classification. Use the
same reference taxonomy/database for the compared contig classifications.
This gate does not compare the candidate's majority label to the old MAG's label;
cluster membership supplies the replacement group, and coherence is internal.

Tables now include `sequence_sha256`, `contig_length` and `taxonomy_rank`.
Matching uses the same canonical full sequence as the overlap check, so renamed
and reverse-complemented contigs can reuse evidence. A reused ID with a different
sequence cannot inherit another contig's taxonomy. Conflicting labels for one
sequence are unresolved and reported, never resolved by last-row precedence.

For **older ID-only classify-contigs outputs**, rerun only the inexpensive
aggregation using the existing raw classifier output and the original matching
FASTA collection:

```bash
echidna classify-contigs --bins-dir all_original_mags/ \
  --mmseqs-output existing_contig_taxonomy.tsv --rank order \
  -o all_contigs.classified.tsv
```

Kraken2/Kaiju inputs with their taxonomy trees remain supported. The raw
classifier output must refer to those exact contig sequences. Before raw
ID-based aggregation, ambiguous IDs shared by different sequences are rejected;
classify separate assemblies separately or use globally unique IDs. Already
sequence-linked tables can be reused across such IDs safely. The MMseqs2
wrapper classifies identical shared IDs once and rejects conflicting sequences.
Unknown candidate contigs require additional classification; absence of evidence
does not produce a passing taxonomy result.

`taxonomy_checks.tsv` records each evaluated subset and the final exported
FASTA checks, including rank, resolved/all-base coverage, foreign fraction,
majority, missing bases and conflicting bases. `replacement_attempts.tsv`
distinguishes `shared_contigs`, `incoherent_taxonomy`, `insufficient_taxonomy`
and acceptance. `resolution.tsv` adds taxonomy status and measured fractions.
If no MAG passes, retained/ is empty and summary.json explicitly records zero
evaluated retained MAGs; the workflow stops before downstream QC.

Both shared-contig and taxonomy checks are repeated on the exported FASTAs.
Contig recovery after another cluster's replacement is also gated: newly freed
contigs are not restored if they would make a previously accepted subset fail.
Final GUNC, CheckM2, GTDB-Tk and RNA checks still apply. Taxonomy coherence is one
signal, not a replacement for the existing two-signal chimera verdict.

To recompute summaries on any new bin collection without rerunning a classifier:

```bash
echidna classify-contigs --bins-dir resolved_checked/retained/ \
  --taxonomy-table all_contigs.classified.tsv --rank order \
  -o resolved_checked/final_contig_taxonomy.tsv
```

The [MAG workflows](workflows.md#overlap-resolution-with-existing-contig-taxonomy)
can perform selection, reuse taxonomy and run final QC in one ordered workflow.
