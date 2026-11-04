#!/usr/bin/env bash
# Two explicit phases: assembly QC, then MAG QC after external mapping/binning.
# Heavy tools run in the environments configured by the maintained *-wf commands.
set -euo pipefail

MODE="${1:-}"
if [[ $# -gt 0 ]]; then shift; fi
WORKDIR="${WORKDIR:-echidna_workflow}"
THREADS="${THREADS:-20}"

case "$MODE" in
  assembly)
    : "${ASSEMBLY:?Set ASSEMBLY to the original assembly FASTA}"
    : "${MAPPING_BAM:?Set MAPPING_BAM to sorted/indexed PacBio HiFi mapping on that assembly}"
    echidna misassembly-wf \
      --assembly "$ASSEMBLY" --bam "$MAPPING_BAM" \
      -o "$WORKDIR/misassembly" --threads "$THREADS" \
      --anvio-env "${ANVIO_ENV:-anvio-9}" --checkm2-env "${CHECKM2_ENV:-checkm2}" "$@"
    cat <<EOF
Assembly QC finished. Review the reports in $WORKDIR/misassembly.
Next, map reads to the cleaned assembly again, calculate coverage and run your
binner externally. Old BAM coordinates/depth cannot be reused for split contigs.
Keep HQ_extracted genomes and residual-assembly bins available; inspect combined
collections for shared contigs before interpretation. Run the mags phase when
the bin FASTAs are ready. No binning is performed by this example.
EOF
    ;;
  mags)
    : "${BINS_DIR:?Set BINS_DIR to the MAG FASTA directory}"
    : "${GTDBTK_DATA_PATH:?Set GTDBTK_DATA_PATH to the GTDB-Tk database directory}"
    : "${GUNC_DB:?Set GUNC_DB to the GUNC diamond database file}"
    args=(--bins-dir "$BINS_DIR" --genome-ext "${GENOME_EXT:-fa}"
          -o "$WORKDIR/mag_quality" --threads "$THREADS"
          --gtdbtk-data "$GTDBTK_DATA_PATH" --gunc-db "$GUNC_DB"
          --checkm2-env "${CHECKM2_ENV:-checkm2}" --gtdbtk-env "${GTDBTK_ENV:-gtdbtk}"
          --gunc-env "${GUNC_ENV:-gunc}" --mimag-env "${MIMAG_ENV:-mimag}")
    if [[ -n "${CONTIG_TAXONOMY:-}" ]]; then
      args+=(--contig-taxonomy "$CONTIG_TAXONOMY" --taxonomy-rank "${TAXONOMY_RANK:-order}")
    elif [[ -n "${MMSEQS_DB:-}" ]]; then
      args+=(--mmseqs-db "$MMSEQS_DB" --mmseqs-env "${MMSEQS_ENV:-mmseqs2}")
    fi
    if [[ -n "${RESOLVE_OVERLAPS:-}" ]]; then
      : "${QUALITY_REPORT:?Selection requires prior quality for representatives and candidates}"
      : "${CONTIG_TAXONOMY:?Selection requires sequence-linked per-contig taxonomy covering candidates}"
      args+=(--resolve-overlaps "$RESOLVE_OVERLAPS" --quality-report "$QUALITY_REPORT")
      if [[ -n "${DREP_DIR:-}" ]]; then args+=(--drep-dir "$DREP_DIR"); fi
      if [[ -n "${CLUSTER_TABLE:-}" ]]; then args+=(--cluster-table "$CLUSTER_TABLE"); fi
      if [[ -n "${CANDIDATES_DIR:-}" ]]; then args+=(--candidates-dir "$CANDIDATES_DIR"); fi
    fi
    echidna mag-quality-wf "${args[@]}" "$@"
    ;;
  *)
    cat >&2 <<EOF
Usage: $0 assembly|mags [additional flags for the matching Echidna workflow]

ASSEMBLY=assembly.fa MAPPING_BAM=hifi.bam $0 assembly
BINS_DIR=bins GTDBTK_DATA_PATH=/data/gtdbtk GUNC_DB=/data/gunc.dmnd $0 mags

Optional MAG selection: RESOLVE_OVERLAPS=exclude-mags (or remove-contigs),
QUALITY_REPORT=all_quality.tsv, CONTIG_TAXONOMY=all_contigs.classified.tsv,
DREP_DIR=drep (or CLUSTER_TABLE=clusters.tsv), CANDIDATES_DIR=all_original_mags.
The per-contig taxonomy is reused; CheckM2 and downstream QC run on final MAGs.
EOF
    exit 2
    ;;
esac
