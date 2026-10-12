#!/usr/bin/env python3
"""
Echidna - metagenomic assembly, MAG integrity and chimera screening
Main entry point for package execution
"""

import sys
import argparse
from pathlib import Path

from echidna.commands import filter_misassemblies
from echidna.commands import remove_hq_contigs
from echidna.commands import filter_checkm2
from echidna.commands import run_checkm2
from echidna.commands import run_gtdbtk
from echidna.commands import run_gunc
from echidna.commands import run_mimag_rna
from echidna.commands import process_large_contigs
from echidna.commands import check_chimeras
from echidna.commands import check_read_chimeras
from echidna.commands import classify_contigs
from echidna.commands import run_mmseqs_contigs
from echidna.commands import check_ssu_taxonomy
from echidna.commands import classify_clipping
from echidna.commands import check_zero_coverage
from echidna.commands import setup_workflow
from echidna.commands import fetch_databases
from echidna.commands import check_envs
from echidna.commands import run_anvi_misassemblies
from echidna.commands import misassembly_wf
from echidna.commands import mag_quality_wf
from echidna import __version__


def main():
    parser = argparse.ArgumentParser(
        prog='echidna',
        description='Echidna - metagenomic assembly, MAG integrity and chimera screening',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Available commands:
  filter-misassemblies    Filter assembly based on misassembly detection
  remove-hq-contigs       Remove high-quality genome contigs from assembly
  organize-mags           Split genomes into HQ/MQ/LQ tiers, optionally rename by sample+taxonomy
                          (alias: filter-checkm2, deprecated)
  run-checkm2             Run CheckM2 completeness/contamination prediction
  run-gtdbtk              Run GTDB-Tk taxonomy classification
  run-mimag-rna           Check MIMAG rRNA (5S/16S/23S) + tRNA criterion (barrnap + tRNAscan-SE)
  process-large-contigs   Extract, QC, and filter large contigs intelligently
  run-gunc                Run GUNC gene-level chimerism/contamination detection
  check-chimeras          Detect chimeric MAGs and circular contigs (GC + GUNC)
  check-read-chimeras     Detect chimeric contigs via read-level taxonomy (Kraken2+BAM)
  classify-contigs        Per-contig taxonomy check: flag bins whose contigs span >=2 taxa
  run-mmseqs-contigs      Run MMseqs2 taxonomy on bin contigs (input for classify-contigs)
  check-ssu-taxonomy      Classify each MAG's 16S and check it agrees with GTDB-Tk (blastn)
  classify-clipping       Classify each clipping position with BAM evidence
  check-zero-coverage     Validate assembly regions with no read coverage via BLAST
  setup-workflow          Configure the Snakemake workflow (env yamls + config)
  fetch-databases         Download reference DBs, or register existing ones
  check-envs              Sanity-check existing conda envs (versions + compatibility)

Examples:
  echidna filter-misassemblies --help
  echidna filter-checkm2 --help
  
For more information, visit: https://github.com/l-gallucci/echidna
        """
    )
    
    parser.add_argument('-v', '--version', action='version', 
                        version=f'Echidna version {__version__}')
    
    subparsers = parser.add_subparsers(dest='command', help='Command to run')
    
    # Register all subcommands
    filter_misassemblies.register_parser(subparsers)
    remove_hq_contigs.register_parser(subparsers)
    filter_checkm2.register_parser(subparsers)
    run_checkm2.register_parser(subparsers)
    run_gtdbtk.register_parser(subparsers)
    run_gunc.register_parser(subparsers)
    run_mimag_rna.register_parser(subparsers)
    process_large_contigs.register_parser(subparsers)
    check_chimeras.register_parser(subparsers)
    check_read_chimeras.register_parser(subparsers)
    classify_contigs.register_parser(subparsers)
    run_mmseqs_contigs.register_parser(subparsers)
    check_ssu_taxonomy.register_parser(subparsers)
    classify_clipping.register_parser(subparsers)
    check_zero_coverage.register_parser(subparsers)
    setup_workflow.register_parser(subparsers)
    fetch_databases.register_parser(subparsers)
    check_envs.register_parser(subparsers)
    run_anvi_misassemblies.register_parser(subparsers)
    misassembly_wf.register_parser(subparsers)
    mag_quality_wf.register_parser(subparsers)
    
    # Parse arguments
    args = parser.parse_args()
    
    if not args.command:
        parser.print_help()
        sys.exit(0)
    
    # Execute command
    try:
        args.func(args)
    except KeyboardInterrupt:
        print("\n[!] Interrupted by user", file=sys.stderr)
        sys.exit(130)
    except Exception as e:
        print(f"\n[ERROR] {e}", file=sys.stderr)
        if '--debug' in sys.argv:
            raise
        sys.exit(1)


if __name__ == '__main__':
    main()
