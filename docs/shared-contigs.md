# Shared contigs between MAGs

Use `check-shared-contigs` to check whether genome FASTAs contain the same
full-length contig sequence. It also reports reused contig IDs whose sequences
differ. The command uses the Python standard library and needs no database.

```bash
echidna check-shared-contigs --bins-dir bins/ --output-dir overlap/
```

This is a standalone report. It does not change, remove, merge or reassign any
input sequence and does not alter the chimera verdict. An overlap alone does
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
nonzero status.

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

## Relation to Binette

[Binette](https://binette.readthedocs.io/en/stable/) performs bin refinement:
it constructs candidate bins using intersections, differences and unions of
contig sets and evaluates candidates with CheckM2. It is a separate tool and
is not required by `check-shared-contigs`.

Use this Echidna command to inspect shared sequences in an existing MAG set.
Use Binette separately when you want to refine bin sets from multiple binners;
its input contig identities must correspond to the assembly being refined.
