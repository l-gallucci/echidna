"""
Per-MAG final report: join the signals Echidna already produces into one table.

Each source is optional; a field is 'NA' when its source was not given. The
join key is the genome name with a common FASTA extension stripped, so the
CheckM2 'Name', the GTDB-Tk 'user_genome', the bin filenames etc. line up.

A 'concerns' column flags the red signals (high contamination, High chimera
risk, a per-contig taxonomic split, a 16S disagreeing with the genome, an
incomplete MIMAG rRNA/tRNA set) so a reviewer can sort on it; it is a summary,
not a hard pass/fail (that policy stays with the user).

Pure parsers + a pure build_report(); no external calls.
"""

import csv
import math
from pathlib import Path
from typing import Dict, List, Optional, Set

from .checkm2 import assign_tier


_EXTS = ('.fa.gz', '.fasta.gz', '.fna.gz', '.fa', '.fasta', '.fna', '.fa.gz')


def norm_name(name: str) -> str:
    """Strip a trailing FASTA extension so names from different tools match."""
    n = (name or '').strip()
    low = n.lower()
    for ext in _EXTS:
        if low.endswith(ext):
            return n[: -len(ext)]
    return n


def _read_tsv(path: Path) -> List[Dict]:
    with open(path, newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def parse_chimera_report(path: Path) -> Dict[str, str]:
    """check-chimeras chimera_report.tsv -> {genome: chimera_risk}."""
    out: Dict[str, str] = {}
    for row in _read_tsv(path):
        g = norm_name(row.get('bin') or row.get('genome') or '')
        if g:
            out[g] = (row.get('chimera_risk') or 'NA').strip()
    return out


def parse_gunc_pass(path: Path) -> Dict[str, str]:
    """check-chimeras chimera_report.tsv -> {genome: gunc_pass} ('True'/'False'/'N/A').
    Used for the GUNC-vs-MMseqs concordance column."""
    out: Dict[str, str] = {}
    for row in _read_tsv(path):
        g = norm_name(row.get('bin') or row.get('genome') or '')
        if g:
            out[g] = (row.get('gunc_pass') or 'N/A').strip()
    return out


def parse_contig_taxonomy_bins(path: Path) -> Dict[str, Dict]:
    """classify-contigs <out>.bins.tsv -> {genome: {contig_chimera, n_foreign}}."""
    out: Dict[str, Dict] = {}
    for row in _read_tsv(path):
        g = norm_name(row.get('bin') or '')
        if not g:
            continue
        out[g] = {
            'contig_chimera': (str(row.get('chimera_flag', '')).strip().lower()
                               in ('true', '1', 'yes')),
            'n_foreign': row.get('n_foreign_contigs', '0').strip() or '0',
        }
    return out


def parse_ssu_report(path: Path) -> Dict[str, Dict]:
    """check-ssu-taxonomy ssu_taxonomy.tsv -> {genome: {agree, ssu_taxon}}.

    Newer reports may contain multiple SSU rows per genome. Aggregate
    conservatively: any disagreement makes the genome-level result False; if
    there are no disagreements, any agreement makes it True; otherwise NA.
    """
    out: Dict[str, Dict] = {}
    for row in _read_tsv(path):
        g = norm_name(row.get('genome') or '')
        if not g:
            continue
        agree = (row.get('agree') or 'NA').strip()
        taxon = (row.get('ssu_taxon') or 'NA').strip()
        if g not in out:
            out[g] = {'agree': agree, 'ssu_taxon': taxon}
            continue
        prev = out[g]
        if agree == 'False' or prev['agree'] == 'False':
            prev['agree'] = 'False'
        elif agree == 'True' or prev['agree'] == 'True':
            prev['agree'] = 'True'
        else:
            prev['agree'] = 'NA'
        if taxon != 'NA':
            taxa = [t for t in prev['ssu_taxon'].split(';') if t and t != 'NA']
            if taxon not in taxa:
                taxa.append(taxon)
            prev['ssu_taxon'] = ';'.join(taxa) if taxa else 'NA'
    return out


def _fmt(v) -> str:
    if v is None:
        return 'NA'
    if isinstance(v, float):
        return 'NA' if math.isnan(v) else f"{v:.2f}"
    return str(v)


def _concerns(row: Dict) -> str:
    c: List[str] = []
    contam = row.get('contamination')
    if isinstance(contam, float) and not math.isnan(contam) and contam > 10.0:
        c.append('contam>10')
    if row.get('chimera_risk') == 'High':
        c.append('chimera_High')
    elif row.get('chimera_risk') == 'Medium':
        c.append('chimera_Medium')
    if row.get('contig_chimera') is True:
        c.append('contig_chimera')
    if row.get('ssu_agree') == 'False':
        c.append('ssu_disagree')
    if row.get('mimag_rna_ok') is False:
        c.append('mimag_rna_incomplete')
    return ';'.join(c)


REPORT_FIELDS = [
    'genome', 'completeness', 'contamination', 'tier', 'taxonomy',
    'rank_taxon', 'n_16S', 'rrna_complete', 'trna_complete', 'mimag_rna_ok',
    'chimera_risk', 'contig_chimera', 'n_foreign_contigs',
    'ssu_agree', 'ssu_taxon', 'gunc_mmseqs_concordance', 'concerns',
]


def _concordance(gunc_pass: Optional[str], contig_chimera) -> str:
    """GUNC vs MMseqs per-contig agreement on chimeric-or-not: True/False/NA.
    GUNC says chimeric when pass.GUNC is False; MMseqs when contig_chimera."""
    if gunc_pass in (None, '', 'N/A', 'NA') or not isinstance(contig_chimera, bool):
        return 'NA'
    gunc_chimeric = gunc_pass.strip().lower() in ('false', '0', 'no')
    return str(gunc_chimeric == contig_chimera)


def build_report(
    universe: Set[str],
    checkm2: Optional[Dict[str, Dict]] = None,
    taxonomy: Optional[Dict[str, Dict]] = None,   # genome -> {classification, rank_taxon}
    mimag: Optional[Dict[str, Dict]] = None,
    chimera: Optional[Dict[str, str]] = None,
    contig: Optional[Dict[str, Dict]] = None,
    ssu: Optional[Dict[str, Dict]] = None,
    gunc_pass: Optional[Dict[str, str]] = None,
) -> List[Dict]:
    """Join all provided per-genome maps (keyed by normalised name) into rows."""
    checkm2 = checkm2 or {}
    taxonomy = taxonomy or {}
    mimag = mimag or {}
    chimera = chimera or {}
    contig = contig or {}
    ssu = ssu or {}
    gunc_pass = gunc_pass or {}

    rows: List[Dict] = []
    for genome in sorted(universe):
        cm = checkm2.get(genome, {})
        tx = taxonomy.get(genome, {})
        mi = mimag.get(genome, {})
        co = contig.get(genome, {})
        su = ssu.get(genome, {})

        comp = cm.get('completeness')
        contam = cm.get('contamination')
        tier = (assign_tier(comp, contam)
                if isinstance(comp, float) and isinstance(contam, float)
                and not math.isnan(comp) and not math.isnan(contam) else 'NA')

        row = {
            'genome': genome,
            'completeness': comp,
            'contamination': contam,
            'tier': tier,
            'taxonomy': tx.get('classification', 'NA') or 'NA',
            'rank_taxon': tx.get('rank_taxon', 'NA') or 'NA',
            'n_16S': mi.get('n_16S', 'NA') if mi else 'NA',
            'rrna_complete': mi.get('rrna_complete') if mi else 'NA',
            'trna_complete': mi.get('trna_complete') if mi else 'NA',
            'mimag_rna_ok': mi.get('mimag_rna_ok') if mi else 'NA',
            'chimera_risk': chimera.get(genome, 'NA'),
            'contig_chimera': co.get('contig_chimera', 'NA') if co else 'NA',
            'n_foreign_contigs': co.get('n_foreign', 'NA') if co else 'NA',
            'ssu_agree': str(su.get('agree', 'NA')) if su else 'NA',
            'ssu_taxon': su.get('ssu_taxon', 'NA') if su else 'NA',
            'gunc_mmseqs_concordance': _concordance(
                gunc_pass.get(genome),
                co.get('contig_chimera') if co else None),
        }
        row['concerns'] = _concerns(row)
        rows.append(row)
    return rows


def write_report(rows: List[Dict], out_path: Path) -> None:
    with open(out_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=REPORT_FIELDS, delimiter='\t')
        w.writeheader()
        for r in rows:
            w.writerow({k: _fmt(r.get(k)) for k in REPORT_FIELDS})
