#!/usr/bin/env python3
"""
MaMISA - run-anvi-misassemblies command

Thin wrapper around anvi'o's `anvi-script-find-misassemblies`, which detects
soft-clipping positions (misassembly signal) from long reads mapped to the
assembly they were built from (Trigodet et al. 2025). Produces
`<output-prefix>-clipping.txt`, the input for process-large-contigs,
classify-clipping and filter-misassemblies.

Runs inside the anvi'o env (like run-checkm2/run-gunc run inside theirs):
`conda run -n anvio-9 mamisa run-anvi-misassemblies ...`.

PacBio HiFi only: the method needs long reads. Map with
`minimap2 -ax map-hifi`, keep supplementary alignments, sort and index.
"""

import sys
import shlex
import subprocess
import argparse
from pathlib import Path
from typing import Optional

from ..utils.validation import check_dependencies
from ..utils.logging import log_info, log_error, log_warning, print_header, print_section

ANVI_EXE = 'anvi-script-find-misassemblies'


def check_anvi() -> bool:
    deps = check_dependencies([ANVI_EXE])
    if not deps[ANVI_EXE]:
        log_error(f"{ANVI_EXE} not found in PATH. Run inside the anvi'o env, e.g. "
                  f"`conda run -n anvio-9 mamisa run-anvi-misassemblies ...`")
        return False
    log_info(f"Found {ANVI_EXE}: {deps[ANVI_EXE]}")
    return True


def run_anvi(bam: Path, output_prefix: Path, threads: int,
             clipping_ratio: float, min_dist_to_end: int,
             just_do_it: bool, extra_args: Optional[list] = None) -> int:
    """Run anvi-script-find-misassemblies. Returns exit code."""
    cmd = [
        ANVI_EXE,
        '-b', str(bam),
        '-o', str(output_prefix),
        '-T', str(threads),
        '--clipping-ratio', str(clipping_ratio),
        '--min-dist-to-end', str(min_dist_to_end),
    ]
    # anvi'o refuses to overwrite / proceed without this in non-interactive runs.
    if just_do_it:
        cmd.append('--just-do-it')
    if extra_args:
        cmd.extend(extra_args)

    log_info("\nRunning command:")
    log_info("  " + " ".join(shlex.quote(c) for c in cmd))

    try:
        result = subprocess.run(cmd, check=True, stdin=subprocess.DEVNULL)
        return result.returncode
    except subprocess.CalledProcessError as e:
        log_error(f"{ANVI_EXE} failed with exit code {e.returncode}")
        return e.returncode
    except FileNotFoundError:
        log_error(f"{ANVI_EXE} not found. Run inside the anvi'o env.")
        return 127
    except KeyboardInterrupt:
        log_warning("\nInterrupted by user")
        return 130


def register_parser(subparsers):
    parser = subparsers.add_parser(
        'run-anvi-misassemblies',
        help="Run anvi'o anvi-script-find-misassemblies (soft-clipping detection)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
PacBio HiFi only. Map HiFi reads to their own assembly with
`minimap2 -ax map-hifi`, KEEP supplementary alignments (do not -F 2048), sort and
index, then:

  conda run -n anvio-9 mamisa run-anvi-misassemblies \\
    --bam hifi_mapping.bam \\
    --output-prefix 01_misassemblies/MisAsm \\
    --threads 20

Produces 01_misassemblies/MisAsm-clipping.txt. --just-do-it is on by default
(anvi'o refuses to proceed without it non-interactively); pass --no-just-do-it to
turn it off.
        """
    )
    parser.add_argument('-b', '--bam', type=Path, required=True,
                        help='Sorted, indexed BAM of long reads mapped to their assembly')
    parser.add_argument('-o', '--output-prefix', type=Path, required=True,
                        help='Output prefix; writes <prefix>-clipping.txt '
                             '(its parent directory is the misassemblies dir)')
    parser.add_argument('-T', '--threads', type=int, default=1,
                        help='Threads (default: 1)')
    parser.add_argument('--clipping-ratio', type=float, default=0.8,
                        help='anvi clipping ratio threshold (default: 0.8)')
    parser.add_argument('--min-dist-to-end', type=int, default=100,
                        help='anvi minimum distance to contig end (default: 100)')
    parser.add_argument('--no-just-do-it', dest='just_do_it', action='store_false',
                        help="Do not pass --just-do-it to anvi'o")
    parser.set_defaults(just_do_it=True)
    parser.add_argument('--anvi-args', default='',
                        help='Extra arguments passed to anvi-script-find-misassemblies (quoted)')
    parser.set_defaults(func=run)
    return parser


def run(args):
    print_header("MaMISA - Run anvi-script-find-misassemblies")

    print_section("Checking Dependencies")
    if not check_anvi():
        sys.exit(1)

    if not args.bam.exists():
        log_error(f"BAM not found: {args.bam}")
        sys.exit(1)

    out_dir = args.output_prefix.parent
    if str(out_dir) and not out_dir.exists():
        out_dir.mkdir(parents=True, exist_ok=True)

    extra_args = shlex.split(args.anvi_args) if args.anvi_args else []

    print_section("Running anvi-script-find-misassemblies")
    exit_code = run_anvi(
        args.bam, args.output_prefix, args.threads,
        args.clipping_ratio, args.min_dist_to_end,
        args.just_do_it, extra_args,
    )

    if exit_code != 0:
        log_error(f"\nFailed (exit {exit_code})")
        sys.exit(exit_code)

    clipping = Path(f"{args.output_prefix}-clipping.txt")
    log_info(f"\nDone. Clipping file: {clipping}")
    log_info("Feed its directory to process-large-contigs -m / classify-clipping -m / "
             "filter-misassemblies -m")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    register_parser(parser.add_subparsers(dest='command'))
    run(parser.parse_args())
