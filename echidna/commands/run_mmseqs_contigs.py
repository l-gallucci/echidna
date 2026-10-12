#!/usr/bin/env python3
"""
Echidna - run-mmseqs-contigs command

Thin wrapper that runs MMseqs2 taxonomy on the contigs of a set of bins and
writes the per-contig TSV that `classify-contigs --mmseqs-output` consumes.

MMseqs2 is the per-contig classifier MAGpurify2 uses: protein-level 2bLCA, fast
and sensitive on divergent/novel environmental taxa, and its reference DB can be
GTDB-derived (via `mmseqs createtaxdb` / `mmseqs databases`), so the ranks line
up with GTDB-Tk downstream. It is a heavy dependency, so it lives in its own
`mmseqs2` conda env (like GUNC/anvi'o live in theirs):

  conda run -n mmseqs2 echidna run-mmseqs-contigs ...

Pipeline (per MMseqs2 docs):
  mmseqs createdb  <contigs.fa> queryDB
  mmseqs taxonomy  queryDB <seqTaxDB> taxResult tmp --threads N
  mmseqs createtsv queryDB taxResult <output.tsv>

<output.tsv> columns: contig  taxid  rank  name  [...].
"""

import sys
import shlex
import shutil
import subprocess
import argparse
from pathlib import Path
from typing import List, Tuple

from ..utils.fasta import read_fasta_streaming, write_fasta
from ..utils.validation import check_dependencies, validate_dir_exists, validate_file_exists
from ..utils.logging import log_info, log_error, log_warning, print_header, print_section

MMSEQS_EXE = 'mmseqs'


def check_mmseqs(exe: str = MMSEQS_EXE) -> bool:
    deps = check_dependencies([exe])
    if not deps[exe]:
        log_error(f"{exe} not found in PATH. Run inside the mmseqs2 env, e.g. "
                  f"`conda run -n mmseqs2 echidna run-mmseqs-contigs ...`")
        return False
    log_info(f"Found {exe}: {deps[exe]}")
    return True


def concat_bins(bins_dir: Path, extensions: List[str], out_fasta: Path) -> int:
    """Write every contig of every bin into one FASTA. Returns contig count.

    Contig headers are kept verbatim so classify-contigs can map them back to
    their bins; a name shared by two bins would collide here (rare for MAGs).
    """
    def records():
        bin_files: List[Path] = []
        for ext in extensions:
            bin_files.extend(bins_dir.rglob(f"*.{ext}"))
        for bin_file in sorted(set(bin_files)):
            for name, seq in read_fasta_streaming(bin_file):
                yield name, seq

    recs = list(records())
    write_fasta(recs, out_fasta)
    return len(recs)


def build_commands(exe: str, query_fasta: Path, seq_tax_db: Path,
                   output_tsv: Path, tmp_dir: Path, threads: int,
                   extra_args: List[str]) -> List[Tuple[str, List[str]]]:
    """Pure: the three MMseqs2 steps as (label, argv). Intermediate DBs live
    next to tmp_dir so nothing is left in the working directory."""
    query_db = tmp_dir / "queryDB"
    tax_result = tmp_dir / "taxResult"
    mmseqs_tmp = tmp_dir / "tmp"
    # Emit the full lineage column so classify-contigs can resolve any rank
    # without an NCBI nodes.dmp (works with a GTDB DB). Skip if the user set it.
    lineage = [] if '--tax-lineage' in extra_args else ['--tax-lineage', '1']
    return [
        ("createdb", [exe, 'createdb', str(query_fasta), str(query_db)]),
        ("taxonomy", [exe, 'taxonomy', str(query_db), str(seq_tax_db),
                      str(tax_result), str(mmseqs_tmp),
                      '--threads', str(threads)] + lineage + extra_args),
        ("createtsv", [exe, 'createtsv', str(query_db), str(tax_result),
                       str(output_tsv)]),
    ]


def register_parser(subparsers):
    parser = subparsers.add_parser(
        'run-mmseqs-contigs',
        help='Run MMseqs2 taxonomy on bin contigs -> per-contig TSV for classify-contigs',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Needs a prebuilt MMseqs2 taxonomy database (--db), e.g. a GTDB one:
  conda run -n mmseqs2 mmseqs databases GTDB /data/mmseqs_gtdb tmp

Then:
  conda run -n mmseqs2 echidna run-mmseqs-contigs \\
    --bins-dir bins/ --db /data/mmseqs_gtdb -o contig_taxonomy.mmseqs.tsv --threads 40

Feed the result to classify-contigs (the TSV carries a lineage column, so no
nodes.dmp is needed and ranks stay in GTDB's system):
  echidna classify-contigs --bins-dir bins/ \\
    --mmseqs-output contig_taxonomy.mmseqs.tsv --rank class
        """
    )
    parser.add_argument('--bins-dir', type=Path, required=True,
                        help='Directory containing bin FASTA files (contigs are concatenated)')
    parser.add_argument('--db', type=Path, required=True,
                        help='MMseqs2 taxonomy sequence DB (mmseqs databases / createtaxdb)')
    parser.add_argument('-o', '--output', type=Path, required=True,
                        help='Output per-contig taxonomy TSV')
    parser.add_argument('-T', '--threads', type=int, default=1,
                        help='Threads (default: 1)')
    parser.add_argument('--tmp-dir', type=Path,
                        help='Working directory for MMseqs2 DBs/tmp '
                             '(default: <output>.mmseqs_work)')
    parser.add_argument('--extensions', default='fa,fasta,fna',
                        help='Comma-separated bin file extensions (default: fa,fasta,fna)')
    parser.add_argument('--mmseqs-exe', default=MMSEQS_EXE,
                        help='MMseqs2 executable name/path (default: mmseqs)')
    parser.add_argument('--keep-tmp', action='store_true',
                        help='Keep the MMseqs2 working directory (default: remove it)')
    parser.add_argument('--mmseqs-args', default='',
                        help='Extra arguments passed to `mmseqs taxonomy` (quoted)')
    parser.add_argument('--dry-run', action='store_true',
                        help='Print the MMseqs2 commands without running them')
    parser.set_defaults(func=run)
    return parser


def run(args):
    print_header("Echidna - Run MMseqs2 taxonomy on contigs")

    validate_dir_exists(args.bins_dir, "Bins directory")

    if not args.dry_run:
        print_section("Checking dependencies")
        if not check_mmseqs(args.mmseqs_exe):
            sys.exit(1)

    extensions = [e.strip() for e in args.extensions.split(',')]
    tmp_dir = args.tmp_dir or Path(str(args.output) + '.mmseqs_work')
    extra_args = shlex.split(args.mmseqs_args) if args.mmseqs_args else []
    query_fasta = tmp_dir / "contigs.fa"

    print_section("Concatenating bin contigs")
    if args.dry_run:
        log_info(f"(dry-run) would concatenate {args.bins_dir}/*.{{{','.join(extensions)}}} "
                 f"-> {query_fasta}")
    else:
        tmp_dir.mkdir(parents=True, exist_ok=True)
        n = concat_bins(args.bins_dir, extensions, query_fasta)
        if n == 0:
            log_error(f"No contigs found in {args.bins_dir} with extensions: {extensions}")
            sys.exit(1)
        log_info(f"Wrote {n:,} contigs to {query_fasta}")

    steps = build_commands(args.mmseqs_exe, query_fasta, args.db, args.output,
                           tmp_dir, args.threads, extra_args)

    print_section("Running MMseqs2")
    for label, cmd in steps:
        printable = " ".join(shlex.quote(c) for c in cmd)
        log_info(f"\n[{label}] {printable}")
        if args.dry_run:
            continue
        try:
            subprocess.run(cmd, check=True, stdin=subprocess.DEVNULL)
        except subprocess.CalledProcessError as e:
            log_error(f"mmseqs {label} failed with exit code {e.returncode}")
            sys.exit(e.returncode)
        except FileNotFoundError:
            log_error(f"{args.mmseqs_exe} not found. Run inside the mmseqs2 env.")
            sys.exit(127)
        except KeyboardInterrupt:
            log_warning("\nInterrupted by user")
            sys.exit(130)

    if args.dry_run:
        log_warning("\nDRY-RUN: nothing executed")
        return

    if not args.keep_tmp:
        # keep the final TSV (written to args.output, outside tmp_dir); drop DBs
        shutil.rmtree(tmp_dir, ignore_errors=True)
        log_info(f"Removed working dir {tmp_dir} (use --keep-tmp to retain)")

    log_info(f"\nDone. Per-contig taxonomy: {args.output}")
    log_info("Feed it to: echidna classify-contigs --mmseqs-output "
             f"{args.output} --nodes-dmp <nodes.dmp> --rank order")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    register_parser(parser.add_subparsers(dest='command'))
    run(parser.parse_args())
