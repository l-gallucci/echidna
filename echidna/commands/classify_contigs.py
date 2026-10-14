#!/usr/bin/env python3
"""
Echidna - classify-contigs command

Per-contig taxonomy check for bins, independent of GUNC. A bin whose contigs
resolve to two or more taxa at a chosen rank is a taxonomic chimera: contigs of
different taxonomy do not belong in one MAG at any rank.

This is the same idea as MAGpurify2's taxonomy module (per-contig taxonomy ->
bin consensus -> flag the minority/foreign contigs), kept report-only and
dependency-light: it consumes a contig-level Kraken2 or Kaiju classification and
the NCBI taxonomy tree (nodes.dmp), both of which ship with those databases.

How to produce the contig classification
----------------------------------------
Run the classifier on the bin contigs (not reads), one record per contig:

  # concatenate the bins (or classify a single assembly), then:
  kraken2 --db DB --threads 40 --output contigs.kraken2 all_contigs.fa
  # or
  kaiju  -t nodes.dmp -f db.fmi -i all_contigs.fa -o contigs.kaiju

The "sequence id" column Kraken2/Kaiju emit is the contig name, which is matched
back to the bin FASTAs here.
"""

import csv
import sys
import argparse
from pathlib import Path
from typing import Dict, List, Optional

from ..utils.read_taxonomy import parse_kraken2_output, parse_kaiju_output, parse_names_dmp
from ..utils.contig_taxonomy import (
    parse_nodes_dmp,
    parse_mmseqs_taxonomy,
    parse_mmseqs_lineage,
    assign_contigs_to_bins,
    collect_contig_lengths,
    summarize_bin_taxa,
    summarize_bin_lineage,
    VALID_RANKS,
)
from ..utils.validation import validate_dir_exists, validate_file_exists
from ..utils.logging import log_info, log_warning, log_error, print_header, print_section


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

CONTIG_FIELDS = ['contig', 'bin', 'contig_taxid', 'rank_taxon', 'is_foreign']
BIN_FIELDS = ['bin', 'n_contigs', 'n_contigs_resolved', 'n_taxa',
              'majority_taxon', 'chimera_flag', 'foreign_frac',
              'n_foreign_contigs', 'foreign_contigs', 'taxa']


def build_rows(summary: Dict[str, Dict], contig_taxid: Dict[str, int],
               min_foreign_frac: float = 0.0):
    """Flatten summarize output into per-contig + per-bin rows.

    `chimera_flag` is True only when the bin spans >= 2 taxa AND the foreign
    length fraction reaches `min_foreign_frac` (so a tiny foreign tail is listed
    but not flagged). `foreign_frac` is reported either way for downstream use
    (e.g. check-chimeras weighting)."""
    contig_rows: List[Dict] = []
    bin_rows: List[Dict] = []
    for bin_name in sorted(summary):
        rec = summary[bin_name]
        majority = rec['majority_taxon']
        for contig in sorted(rec['per_contig']):
            rank_taxid, rank_name = rec['per_contig'][contig]
            is_foreign = (rank_taxid is not None and rank_taxid != majority)
            contig_rows.append({
                'contig': contig,
                'bin': bin_name,
                'contig_taxid': contig_taxid.get(contig, 0),
                'rank_taxon': rank_name,
                'is_foreign': is_foreign,
            })
        frac = rec.get('foreign_frac', 0.0)
        bin_rows.append({
            'bin': bin_name,
            'n_contigs': rec['n_contigs'],
            'n_contigs_resolved': rec['n_contigs_resolved'],
            'n_taxa': rec['n_taxa'],
            'majority_taxon': rec['majority_name'],
            'chimera_flag': bool(rec['chimera'] and frac >= min_foreign_frac),
            'foreign_frac': f"{frac:.4f}",
            'n_foreign_contigs': len(rec['foreign_contigs']),
            'foreign_contigs': ';'.join(c for c, _ in rec['foreign_contigs']),
            'taxa': ';'.join(rec['taxa']),
        })
    return contig_rows, bin_rows


def write_tsv(rows: List[Dict], fields: List[str], path: Path):
    with open(path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter='\t',
                                extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)
    log_info(f"Written: {path} ({len(rows):,} rows)")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def register_parser(subparsers):
    parser = subparsers.add_parser(
        'classify-contigs',
        help='Per-contig taxonomy check: flag bins whose contigs span >=2 taxa '
             '(independent of GUNC; Kraken2/Kaiju + NCBI nodes.dmp)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
What it does
------------
Assigns each contig its classifier taxon, walks it up the NCBI tree to the
chosen --rank, and per bin counts the distinct rank-level taxa. >= 2 distinct
taxa -> the bin is flagged as a taxonomic chimera. The bin's MAJORITY taxon is
its identity; the minority (foreign) contigs are listed as the likely
contaminants. Contigs classified only ABOVE the rank are left unresolved (they
cannot prove a split).

Pick --rank by how fine a split you want to catch (order is a good default:
phylum is too coarse, genus/species gets noisy on short contigs).

This is report-only (no contigs removed), orthogonal to the GUNC/GC signals in
check-chimeras.

Examples
--------
  # Kraken2 classification of the bin contigs
  echidna classify-contigs \\
    --bins-dir bins/ \\
    --kraken2-output contigs.kraken2 \\
    --nodes-dmp /data/kraken2_db/taxonomy/nodes.dmp \\
    --names-dmp /data/kraken2_db/taxonomy/names.dmp \\
    --rank order \\
    --output contig_taxonomy.tsv

  # Kaiju (protein-level; more sensitive on divergent taxa)
  echidna classify-contigs \\
    --bins-dir bins/ \\
    --kaiju-output contigs.kaiju \\
    --nodes-dmp /data/kaiju_db/nodes.dmp \\
    --names-dmp /data/kaiju_db/names.dmp \\
    --rank order

  # MMseqs2 with a GTDB DB (via run-mmseqs-contigs, which adds --tax-lineage):
  # the lineage column is used, so NO nodes.dmp is needed and ranks are GTDB's.
  echidna classify-contigs \\
    --bins-dir bins/ \\
    --mmseqs-output contig_tax.mmseqs.tsv \\
    --rank class
        """
    )

    parser.add_argument('--bins-dir', type=Path, required=True,
                        help='Directory containing bin/MAG FASTA files')
    parser.add_argument('-o', '--output', type=Path,
                        default=Path('contig_taxonomy.tsv'),
                        help='Per-contig output TSV (default: contig_taxonomy.tsv). '
                             'A per-bin summary is written alongside as '
                             '<output>.bins.tsv')

    taxo = parser.add_mutually_exclusive_group(required=True)
    taxo.add_argument('--kraken2-output', type=Path,
                      help='Kraken2 per-contig output (kraken2 --output on the contigs)')
    taxo.add_argument('--kaiju-output', type=Path,
                      help='Kaiju per-contig output (kaiju -o on the contigs); '
                           'protein-level, more sensitive on divergent taxa')
    taxo.add_argument('--mmseqs-output', type=Path,
                      help='MMseqs2 taxonomy TSV (mmseqs createtsv on a taxonomy '
                           'result, or run-mmseqs-contigs output); protein-level, '
                           'GTDB-capable - best for consistency with GTDB-Tk')

    parser.add_argument('--nodes-dmp', type=Path,
                        help='NCBI taxonomy nodes.dmp (ships with the Kraken2/Kaiju DB). '
                             'Required for Kraken2/Kaiju and for an MMseqs2 TSV without a '
                             'lineage column; NOT needed for MMseqs2 run with '
                             '--tax-lineage (e.g. a GTDB DB), where the lineage is used.')
    parser.add_argument('--names-dmp', type=Path,
                        help='NCBI taxonomy names.dmp for human-readable taxon names '
                             '(taxid/nodes.dmp path only)')
    parser.add_argument('--rank', default='order', choices=list(VALID_RANKS),
                        help='Rank at which to compare per-contig taxa (default: order)')
    parser.add_argument('--min-foreign-frac', type=float, default=0.05,
                        help='Min foreign LENGTH fraction (0-1) for chimera_flag=True '
                             '(default: 0.05). Below it the foreign contigs are still '
                             'listed and foreign_frac reported, but the bin is not '
                             'flagged - a few genes / one short contig is likely HGT.')
    parser.add_argument('--extensions', default='fa,fasta,fna',
                        help='Comma-separated bin file extensions (default: fa,fasta,fna)')
    parser.add_argument('--dry-run', action='store_true',
                        help='Analyse but do not write output files')

    parser.set_defaults(func=run)
    return parser


def run(args):
    """Execute classify-contigs."""
    validate_dir_exists(args.bins_dir, "Bins directory")
    if args.names_dmp and not args.names_dmp.exists():
        log_error(f"names.dmp not found: {args.names_dmp}")
        sys.exit(1)

    if args.mmseqs_output is not None:
        validate_file_exists(args.mmseqs_output, "MMseqs2 taxonomy TSV")
    elif args.kaiju_output is not None:
        validate_file_exists(args.kaiju_output, "Kaiju output file")
    else:
        validate_file_exists(args.kraken2_output, "Kraken2 output file")

    extensions = [e.strip() for e in args.extensions.split(',')]

    print_header("Echidna - Per-Contig Taxonomy Check")

    print_section("Loading classification")
    contig_lineage: Dict[str, str] = {}
    if args.mmseqs_output is not None:
        contig_taxid = parse_mmseqs_taxonomy(args.mmseqs_output)
        contig_lineage = parse_mmseqs_lineage(args.mmseqs_output)
    elif args.kaiju_output is not None:
        contig_taxid = parse_kaiju_output(args.kaiju_output)
    else:
        contig_taxid = parse_kraken2_output(args.kraken2_output)
    log_info(f"Loaded taxon assignments for {len(contig_taxid):,} contigs")

    use_lineage = bool(contig_lineage)
    if not use_lineage:
        # taxid path needs nodes.dmp to walk to the rank
        if not args.nodes_dmp:
            log_error("--nodes-dmp is required for Kraken2/Kaiju, or for an MMseqs2 "
                      "TSV without a lineage column. Re-run run-mmseqs-contigs "
                      "(it adds --tax-lineage), or pass --nodes-dmp.")
            sys.exit(1)
        validate_file_exists(args.nodes_dmp, "nodes.dmp")

    print_section("Mapping contigs to bins")
    bin_contigs = assign_contigs_to_bins(args.bins_dir, extensions)
    if not bin_contigs:
        log_error(f"No bin files found in {args.bins_dir} with extensions: {extensions}")
        sys.exit(1)
    contig_len = collect_contig_lengths(args.bins_dir, extensions)

    print_section(f"Resolving per-contig taxonomy at rank '{args.rank}'")
    if use_lineage:
        log_info("Using MMseqs2 lineage column (no nodes.dmp needed)")
        summary = summarize_bin_lineage(bin_contigs, contig_lineage, args.rank,
                                        contig_len=contig_len)
    else:
        nodes = parse_nodes_dmp(args.nodes_dmp)
        names = parse_names_dmp(args.names_dmp) if args.names_dmp else None
        summary = summarize_bin_taxa(bin_contigs, contig_taxid, args.rank, nodes,
                                     names, contig_len=contig_len)

    contig_rows, bin_rows = build_rows(summary, contig_taxid,
                                       min_foreign_frac=args.min_foreign_frac)

    if not args.dry_run:
        print_section("Writing output")
        write_tsv(contig_rows, CONTIG_FIELDS, args.output)
        bins_path = Path(str(args.output) + '.bins.tsv')
        write_tsv(bin_rows, BIN_FIELDS, bins_path)
    else:
        log_warning("DRY-RUN: output files not written")

    # Summary
    print_section("CONTIG TAXONOMY SUMMARY")
    n_chimeric = sum(1 for r in summary.values() if r['chimera'])
    total = len(summary)
    print(f"  Bins analysed:            {total:>6,}")
    print(f"  Taxonomic chimeras (>=2): {n_chimeric:>6,}")
    if total:
        print(f"  Flagged rate:             {100 * n_chimeric / total:.1f}%")

    flagged = [(b, r) for b, r in summary.items() if r['chimera']]
    if flagged:
        print_section("FLAGGED BINS (majority -> foreign contigs)")
        for bin_name, rec in sorted(flagged)[:20]:
            foreign = ', '.join(c for c, _ in rec['foreign_contigs'][:5])
            more = '' if len(rec['foreign_contigs']) <= 5 else \
                f" (+{len(rec['foreign_contigs']) - 5} more)"
            print(f"  {bin_name}: majority {rec['majority_name']} | "
                  f"{rec['n_taxa']} taxa | foreign: {foreign}{more}")

    if args.dry_run:
        log_warning("\nDRY-RUN mode: no output file written")

    log_info("\n✓ Done!")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    register_parser(parser.add_subparsers(dest='command'))
    run(parser.parse_args())
