#!/usr/bin/env python3
"""
Echidna - check-ssu-taxonomy command

Classify each MAG's 16S/SSU rRNA gene (blastn vs an SSU reference DB) and check
it agrees with the genome-level GTDB-Tk taxonomy. A 16S resolving to a different
taxon than the genome flags misbinning or a contaminating rRNA - an orthogonal,
independent signal to GUNC/GC/MMseqs per-contig taxonomy.

Input 16S comes from `run-mimag-rna --outseq` (barrnap writes
<genome>.rrna.fna per genome); point --ssu-fasta-dir at that rrna_fasta/ dir.

DB: a blastn database made from an SSU reference. A GTDB SSU DB (ssu_all.fna)
compares 1:1 with GTDB-Tk (use --ssu-taxonomy gtdb). SILVA is broader but uses
different names (--ssu-taxonomy silva), so agreement is best-effort at high
ranks. Build with: makeblastdb -in ssu.fna -dbtype nucl -out <db>.

Report-only. Only MAGs with a recovered ~full-length 16S are judged.
"""

import csv
import sys
import argparse
from pathlib import Path
from typing import Dict, List, Optional

from ..utils.ssu_taxonomy import (
    extract_ssu_records, longest_ssu, check_blastn, run_blastn,
    parse_blast_best, ssu_taxon_at_rank, taxa_agree,
)
from ..utils.chimera import parse_gtdbtk_summary, extract_taxonomy_level
from ..utils.fasta import write_fasta
from ..utils.contig_taxonomy import VALID_RANKS
from ..utils.validation import validate_dir_exists
from ..utils.logging import log_info, log_warning, log_error, print_header, print_section


OUT_FIELDS = ['genome', 'has_16s', 'ssu_len', 'ssu_pident', 'ssu_taxon',
              'genome_taxon', 'agree', 'note']


def _find_gtdbtk_summaries(gtdbtk_dir: Path) -> List[Path]:
    cands = sorted(gtdbtk_dir.rglob('gtdbtk.*.summary.tsv'))
    return cands


def load_genome_taxonomy(gtdbtk_dir: Path, rank: str) -> Dict[str, str]:
    """Return {genome: taxon_at_rank} from GTDB-Tk summaries."""
    taxo: Dict[str, str] = {}
    for summary in _find_gtdbtk_summaries(gtdbtk_dir):
        for name, rec in parse_gtdbtk_summary(summary).items():
            taxo[name] = extract_taxonomy_level(rec.get('classification', ''), rank)
    return taxo


def collect_longest_16s(ssu_dir: Path) -> Dict[str, Dict]:
    """Map genome -> {'seq', 'len'} for the longest 16S per <genome>.rrna.fna."""
    out: Dict[str, Dict] = {}
    for fasta in sorted(ssu_dir.rglob('*.rrna.fna')):
        genome = fasta.name[:-len('.rrna.fna')]
        best = longest_ssu(extract_ssu_records(fasta))
        if best is not None:
            out[genome] = {'seq': best[1], 'len': len(best[1])}
        else:
            out[genome] = {'seq': None, 'len': 0}
    return out


def _match_genome_taxon(genome: str, genome_taxa: Dict[str, str]) -> Optional[str]:
    """GTDB-Tk keys may carry/omit the file extension; match tolerantly."""
    if genome in genome_taxa:
        return genome_taxa[genome]
    for key, val in genome_taxa.items():
        if key == genome or key.rsplit('.', 1)[0] == genome:
            return val
    return None


def register_parser(subparsers):
    parser = subparsers.add_parser(
        'check-ssu-taxonomy',
        help="Classify each MAG's 16S and check it agrees with GTDB-Tk (blastn)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Produce the 16S first with run-mimag-rna --outseq:
  echidna run-mimag-rna --genome-dir bins/ -o mimag/ --outseq

Build a blastn SSU DB (GTDB ssu_all.fna recommended, 1:1 with GTDB-Tk):
  makeblastdb -in ssu_all.fna -dbtype nucl -out /data/ssu_gtdb/ssu

Then:
  echidna check-ssu-taxonomy \\
    --ssu-fasta-dir mimag/rrna_fasta \\
    --ssu-db /data/ssu_gtdb/ssu --ssu-taxonomy gtdb \\
    --gtdbtk-dir taxonomy/ --rank order \\
    -o ssu_taxonomy.tsv

With SILVA (broader, but names differ from GTDB -> agreement best-effort):
  ... --ssu-db /data/silva/silva --ssu-taxonomy silva --rank class
        """
    )
    parser.add_argument('--ssu-fasta-dir', type=Path, required=True,
                        help='Directory of <genome>.rrna.fna (run-mimag-rna --outseq)')
    parser.add_argument('--ssu-db', required=True,
                        help='blastn SSU database prefix (makeblastdb -dbtype nucl)')
    parser.add_argument('--gtdbtk-dir', type=Path, required=True,
                        help='GTDB-Tk output directory (genome-level taxonomy to compare)')
    parser.add_argument('--ssu-taxonomy', choices=['gtdb', 'silva'], default='gtdb',
                        help='How to read the DB hit lineage (default: gtdb)')
    parser.add_argument('--rank', default='order', choices=list(VALID_RANKS),
                        help='Rank to compare 16S vs genome taxonomy (default: order)')
    parser.add_argument('-o', '--output', type=Path, default=Path('ssu_taxonomy.tsv'),
                        help='Output TSV (default: ssu_taxonomy.tsv)')
    parser.add_argument('--min-length', type=int, default=1000,
                        help='Minimum 16S length to classify (default: 1000 bp; '
                             'shorter = partial, reported as low-confidence)')
    parser.add_argument('--threads', type=int, default=1,
                        help='blastn threads (default: 1)')
    parser.add_argument('--task', default='megablast',
                        help='blastn task (default: megablast)')
    parser.add_argument('--blastn-exe', default='blastn',
                        help='blastn executable (default: blastn)')
    parser.add_argument('--keep-tmp', action='store_true',
                        help='Keep the combined query FASTA and blast TSV')
    parser.set_defaults(func=run)
    return parser


def run(args):
    print_header("Echidna - 16S/SSU Taxonomy Check")
    validate_dir_exists(args.ssu_fasta_dir, "SSU FASTA directory")
    validate_dir_exists(args.gtdbtk_dir, "GTDB-Tk directory")
    if not check_blastn(args.blastn_exe):
        sys.exit(1)

    print_section("Collecting 16S and genome taxonomy")
    ssu = collect_longest_16s(args.ssu_fasta_dir)
    if not ssu:
        log_error(f"No *.rrna.fna in {args.ssu_fasta_dir} "
                  "(run run-mimag-rna with --outseq first)")
        sys.exit(1)
    genome_taxa = load_genome_taxonomy(args.gtdbtk_dir, args.rank)
    log_info(f"{len(ssu):,} genomes with rRNA FASTA; "
             f"{len(genome_taxa):,} with GTDB-Tk taxonomy")

    # Build one query FASTA (header = genome) for genomes with a usable 16S.
    tmp_query = Path(str(args.output) + '.query.fna')
    tmp_blast = Path(str(args.output) + '.blast.tsv')
    query_records = []
    for genome, rec in ssu.items():
        if rec['seq'] and rec['len'] >= args.min_length:
            query_records.append((genome, rec['seq']))
    if not query_records:
        log_error(f"No 16S >= {args.min_length} bp found; nothing to classify")
        sys.exit(1)
    write_fasta(query_records, tmp_query)

    print_section(f"blastn {len(query_records):,} 16S vs {args.ssu_db}")
    rc = run_blastn(tmp_query, args.ssu_db, tmp_blast, threads=args.threads,
                    task=args.task, exe=args.blastn_exe)
    if rc != 0:
        sys.exit(rc)
    best = parse_blast_best(tmp_blast)

    print_section("Comparing 16S vs genome taxonomy")
    rows: List[Dict] = []
    n_disagree = 0
    for genome in sorted(ssu):
        rec = ssu[genome]
        g_taxon = _match_genome_taxon(genome, genome_taxa)
        if not rec['seq'] or rec['len'] < args.min_length:
            note = 'no_16s' if not rec['seq'] else f"partial_16s(<{args.min_length}bp)"
            rows.append({'genome': genome, 'has_16s': bool(rec['seq']),
                         'ssu_len': rec['len'], 'ssu_pident': 'NA',
                         'ssu_taxon': 'NA', 'genome_taxon': g_taxon or 'NA',
                         'agree': 'NA', 'note': note})
            continue
        hit = best.get(genome)
        if hit is None:
            rows.append({'genome': genome, 'has_16s': True, 'ssu_len': rec['len'],
                         'ssu_pident': 'NA', 'ssu_taxon': 'NA',
                         'genome_taxon': g_taxon or 'NA', 'agree': 'NA',
                         'note': 'no_blast_hit'})
            continue
        ssu_taxon = ssu_taxon_at_rank(hit['stitle'], args.rank, args.ssu_taxonomy)
        agree = taxa_agree(ssu_taxon, g_taxon)
        if agree is False:
            n_disagree += 1
        note = ''
        if hit['pident'] is not None and hit['pident'] < 90:
            note = 'low_identity'
        rows.append({
            'genome': genome, 'has_16s': True, 'ssu_len': rec['len'],
            'ssu_pident': (f"{hit['pident']:.1f}" if hit['pident'] is not None else 'NA'),
            'ssu_taxon': ssu_taxon or 'NA',
            'genome_taxon': g_taxon or 'NA',
            'agree': ('NA' if agree is None else agree),
            'note': note,
        })

    with open(args.output, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=OUT_FIELDS, delimiter='\t')
        w.writeheader()
        w.writerows(rows)
    log_info(f"\nWritten: {args.output} ({len(rows):,} genomes)")

    if not args.keep_tmp:
        tmp_query.unlink(missing_ok=True)
        tmp_blast.unlink(missing_ok=True)

    print_section("SSU TAXONOMY SUMMARY")
    n_classified = sum(1 for r in rows if r['ssu_taxon'] != 'NA')
    n_agree = sum(1 for r in rows if r['agree'] is True)
    print(f"  Genomes:                 {len(rows):>6,}")
    print(f"  With classified 16S:     {n_classified:>6,}")
    print(f"  16S agrees with genome:  {n_agree:>6,}")
    print(f"  16S DISAGREES (flagged): {n_disagree:>6,}")
    if n_disagree:
        print_section("DISAGREEMENTS (16S taxon vs genome taxon)")
        for r in rows:
            if r['agree'] is False:
                print(f"  {r['genome']}: 16S={r['ssu_taxon']} vs genome={r['genome_taxon']} "
                      f"(pident {r['ssu_pident']})")

    log_info("\nDone.")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    register_parser(parser.add_subparsers(dest='command'))
    run(parser.parse_args())
