"""Reuse sequence-linked contig classifications without rerunning a classifier."""
from collections import Counter, defaultdict
import csv
import math
from pathlib import Path
import re

from .bin_overlap import read_bin, sequence_digest
from .contig_taxonomy import VALID_RANKS


def normalized_rank(rank):
    return 'domain' if rank == 'superkingdom' else rank


def validate_thresholds(min_resolved_frac, min_foreign_frac):
    if not math.isfinite(min_resolved_frac) or not 0 < min_resolved_frac <= 1:
        raise ValueError('Minimum taxonomy-resolved fraction must be in (0, 1]')
    if not math.isfinite(min_foreign_frac) or not 0 <= min_foreign_frac <= 1:
        raise ValueError('Minimum foreign fraction must be in [0, 1]')


class TaxonomyEvidence:
    """A single-rank classify-contigs table, matched by canonical sequence hash.

    Old ID-only tables must be regenerated with classify-contigs and the original
    classifier output; this is a cheap aggregation, not a classifier rerun.
    """
    def __init__(self, path, rank=None, min_resolved_frac=0.5, min_foreign_frac=0.05):
        validate_thresholds(min_resolved_frac, min_foreign_frac)
        self.path = Path(path)
        self.min_resolved_frac = min_resolved_frac
        self.min_foreign_frac = min_foreign_frac
        assignments = defaultdict(set)
        self.lengths = {}
        ranks = set()
        with self.path.open(encoding='utf-8-sig', newline='') as stream:
            reader = csv.DictReader(stream, delimiter='\t')
            required = {'sequence_sha256', 'contig_length', 'taxonomy_rank', 'rank_taxon'}
            if not required.issubset(reader.fieldnames or []):
                raise ValueError('Contig taxonomy needs a sequence-linked classify-contigs TSV '
                                 '(sequence_sha256, contig_length, taxonomy_rank, rank_taxon). '
                                 'Regenerate classify-contigs from the existing classifier output; '
                                 'a per-bin .bins.tsv or legacy ID-only table cannot verify candidate sequences.')
            if len(reader.fieldnames) != len(set(reader.fieldnames)):
                raise ValueError('Duplicate contig-taxonomy columns')
            for row in reader:
                if None in row or any(v is None for v in row.values()):
                    raise ValueError('Malformed contig-taxonomy row')
                digest = row['sequence_sha256'].strip().lower()
                if not re.fullmatch('[0-9a-f]{64}', digest):
                    raise ValueError('Invalid sequence_sha256 in contig taxonomy')
                try:
                    length = int(row['contig_length'])
                except ValueError:
                    raise ValueError('Invalid contig_length in contig taxonomy') from None
                if length < 1 or (digest in self.lengths and self.lengths[digest] != length):
                    raise ValueError('Invalid or inconsistent contig_length in contig taxonomy')
                self.lengths[digest] = length
                row_rank = row['taxonomy_rank'].strip()
                if row_rank not in VALID_RANKS:
                    raise ValueError('Invalid taxonomy_rank in contig taxonomy')
                ranks.add(normalized_rank(row_rank))
                taxon = row['rank_taxon'].strip()
                if taxon.lower() not in ('', 'n/a', 'unresolved', 'unclassified', 'unknown', 'none'):
                    assignments[digest].add(taxon)
        if len(ranks) != 1:
            raise ValueError('Contig taxonomy must contain exactly one taxonomic rank and at least one row')
        self.rank = ranks.pop()
        if rank is not None and normalized_rank(rank) != self.rank:
            raise ValueError(f'Contig-taxonomy rank is {self.rank}, requested {rank}; '
                             'regenerate classify-contigs at the requested rank')
        self.conflicts = {digest for digest, taxa in assignments.items() if len(taxa) > 1}
        self.taxa = {digest: next(iter(taxa)) for digest, taxa in assignments.items() if len(taxa) == 1}

    def taxon(self, digest, length):
        if digest in self.lengths and self.lengths[digest] != length:
            raise ValueError('Taxonomy sequence length does not match its fingerprint')
        return self.taxa.get(digest)

    def assess_records(self, records, keep=None):
        weights = Counter()
        total_bp = resolved_bp = missing_bp = conflict_bp = n_contigs = 0
        for _, digest, length in records:
            if keep is not None and digest not in keep:
                continue
            n_contigs += 1
            total_bp += length
            taxon = self.taxon(digest, length)
            if digest not in self.lengths:
                missing_bp += length
            if digest in self.conflicts:
                conflict_bp += length
            if taxon is not None:
                weights[taxon] += length
                resolved_bp += length
        majority = min(weights, key=lambda taxon: (-weights[taxon], taxon)) if weights else ''
        fraction = resolved_bp / total_bp if total_bp else 0.0
        foreign = (resolved_bp - weights[majority]) / resolved_bp if resolved_bp else None
        status = 'coherent'
        if not resolved_bp or fraction < self.min_resolved_frac:
            status = 'insufficient_taxonomy'
        elif len(weights) > 1 and foreign >= self.min_foreign_frac:
            status = 'incoherent_taxonomy'
        return {'taxonomy_status': status, 'taxonomy_rank': self.rank,
                'taxonomy_resolved_fraction_bp': fraction, 'taxonomy_foreign_fraction': foreign,
                'taxonomy_majority': majority, 'taxonomy_n_taxa': len(weights),
                'taxonomy_total_bp': total_bp, 'taxonomy_resolved_bp': resolved_bp,
                'taxonomy_missing_bp': missing_bp, 'taxonomy_conflicting_bp': conflict_bp,
                'taxonomy_n_contigs': n_contigs}

    @staticmethod
    def records(path):
        return [(name, sequence_digest(seq), len(seq)) for name, seq in read_bin(path)]

    def assess_file(self, path):
        return self.assess_records(self.records(path))


class BinTaxonomy:
    """Cache input records; keep every evaluation for the audit trail."""
    def __init__(self, evidence, inventory):
        self.evidence, self.inventory = evidence, inventory
        self.cache, self.checks = {}, []

    def assess(self, key, keep=None, stage='selection'):
        if key not in self.cache:
            self.cache[key] = self.evidence.records(self.inventory.files[key])
        result = self.evidence.assess_records(self.cache[key], keep)
        self.checks.append(dict(genome=key, stage=stage, **result))
        return result
