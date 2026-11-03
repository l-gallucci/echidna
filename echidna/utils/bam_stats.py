"""BAM evidence for clipping classification.

Local and mean depth both count aligned bases from primary plus secondary
alignments, excluding unmapped/supplementary (mask 2052). Local depth is the
mean per-base depth of the clipped-position window. Pair statistics require
paired primary records; unpaired HiFi reads have NA pair statistics. Primary
SA tags provide inter-contig join evidence without counting supplementary
records twice. BAM starts are converted from 1-based SAM to zero-based positions.
"""

import re
import math
import tempfile
import subprocess
import bisect
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Optional, Tuple

from .logging import log_info, log_warning, log_error
from .validation import check_dependencies


# ──────────────────────────────────────────────────────────────────────────────
# Dependency check
# ──────────────────────────────────────────────────────────────────────────────

def check_samtools() -> bool:
    if check_dependencies(['samtools'])['samtools'] is None:
        log_error("samtools not found in PATH — required for BAM stats")
        log_error("Install: conda install -c bioconda samtools")
        return False
    return True


# ──────────────────────────────────────────────────────────────────────────────
# Contig-level mean depth  (samtools coverage)
# ──────────────────────────────────────────────────────────────────────────────

def get_contig_mean_depths(bam_file: Path, min_mapq: int = 20) -> Dict[str, float]:
    """
    Run ``samtools coverage`` and return {contig: mean_depth}.

    The explicit --ff mask matches collect_position_stats; default samtools
    flags would exclude secondary records but not supplementary ones.
    We pass --min-MQ to apply the same MAPQ filter used in the rest of the
    pipeline so all coverage estimates are on a consistent base.
    """
    cmd = ['samtools', 'coverage', '--min-MQ', str(min_mapq), '--ff', '2052', '-d', '0', str(bam_file)]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"samtools coverage failed: {(e.stderr or chr(32))[:300]}") from e

    depths: Dict[str, float] = {}
    for line in result.stdout.splitlines():
        if line.startswith('#') or not line.strip():
            continue
        fields = line.split('\t')
        if len(fields) < 7:
            continue
        contig = fields[0]
        try:
            mean_depth = float(fields[6])   # column 7 = meandepth
        except (ValueError, IndexError):
            continue
        depths[contig] = mean_depth

    log_info(f"samtools coverage: mean depths for {len(depths):,} contigs")
    return depths


# ──────────────────────────────────────────────────────────────────────────────
# Library insert-size estimation
# ──────────────────────────────────────────────────────────────────────────────

def estimate_insert_size(bam_file: Path, n_reads: int = 20_000,
                         min_mapq: int = 20) -> Tuple[float, float]:
    """
    Estimate library mean and std of insert sizes from the first ``n_reads``
    primary, properly-paired, concordant reads.

    Returns (mean, std). Returns (0, 0) for single-end or if estimation fails.
    """
    # -f 2 = properly paired; -F 2308 = not unmapped (4), not mate-unmapped (8),
    #        not secondary (256), not supplementary (2048)
    cmd = ['samtools', 'view', '-f', '2', '-F', '2308',
           '-q', str(min_mapq), str(bam_file)]
    try:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                stderr=subprocess.DEVNULL, text=True)
    except Exception as e:
        log_warning(f"Could not estimate insert size: {e}")
        return 0.0, 0.0

    sizes: List[int] = []
    for line in proc.stdout:
        fields = line.split('\t')
        if len(fields) < 9:
            continue
        try:
            tlen = int(fields[8])
        except ValueError:
            continue
        if tlen > 0:
            sizes.append(tlen)
        if len(sizes) >= n_reads:
            break

    proc.stdout.close()
    proc.wait()

    if not sizes:
        return 0.0, 0.0

    mean = sum(sizes) / len(sizes)
    variance = sum((x - mean) ** 2 for x in sizes) / len(sizes)
    std = math.sqrt(variance)
    log_info(f"Library insert size: mean={mean:.0f} bp, std={std:.0f} bp "
             f"(estimated from {len(sizes):,} pairs)")
    return mean, std


# ──────────────────────────────────────────────────────────────────────────────
# Per-position stats via single BAM stream
# ──────────────────────────────────────────────────────────────────────────────

class _PosAccumulator:
    """Collects read-level evidence for one clipping position."""
    __slots__ = (
        'clip_pos',
        'this_contig',
        # primary + secondary (for depth)
        'n_depth',
        # primary only
        'n_primary',
        'n_proper_pair',
        'n_paired',
        'n_discordant',
        'n_large_insert',
        'n_forward',
        'n_reverse',
        # clipped-base entropy (from reads soft-clipping at exactly this position)
        '_clipped_bases',
        # supplementary-alignment join evidence (Trigodet et al. 2025)
        '_n_clip_reads',       # primary reads soft-clipping near this position
        '_n_clip_with_sa',     # of those, how many carry an SA tag to another contig
        '_sa_targets',         # {other_contig: n_reads whose clip maps there}
    )

    def __init__(self, clip_pos: int, this_contig: str = ''):
        self.clip_pos = clip_pos
        self.this_contig = this_contig
        self.n_depth = 0
        self.n_primary = 0
        self.n_proper_pair = 0
        self.n_paired = 0
        self.n_discordant = 0
        self.n_large_insert = 0
        self.n_forward = 0
        self.n_reverse = 0
        self._clipped_bases: List[str] = []
        self._n_clip_reads = 0
        self._n_clip_with_sa = 0
        self._sa_targets: Dict[str, int] = defaultdict(int)

    def add_read(self, flag: int, tlen: int, seq: str, cigar: str,
                 read_pos: int, insert_mean: float, insert_std: float,
                 sa_partners: Optional[List[str]] = None, depth_contribution: float = 1.0):
        is_secondary = bool(flag & 0x100)
        is_primary = not is_secondary          # (supplementary already excluded upstream)

        # All non-supplementary reads count toward depth
        self.n_depth += depth_contribution

        if not is_primary:
            return

        self.n_primary += 1

        # Proper pair
        if flag & 0x1:
            self.n_paired += 1
            if flag & 0x2:
                self.n_proper_pair += 1
            else:
                self.n_discordant += 1

        # Large insert (only for proper pairs to avoid noise)
        if (flag & 0x2) and insert_mean > 0 and abs(tlen) > insert_mean + 3 * insert_std:
            self.n_large_insert += 1

        # Strand
        if flag & 0x10:
            self.n_reverse += 1
        else:
            self.n_forward += 1

        # Soft-clip boundary near this position?
        # CIGAR 'S' at the start → clip ends at read_pos-1 (position before alignment).
        # CIGAR 'S' at the end  → clip starts at read_pos + aligned_length.
        clipped = ''
        if seq and cigar and cigar != '*':
            clipped = _extract_clipped_bases(cigar, seq, read_pos, self.clip_pos)

        if clipped:
            self._n_clip_reads += 1
            # bound memory: keep at most 200 clipped sequences for entropy
            if len(self._clipped_bases) < 200:
                self._clipped_bases.append(clipped)

            # Supplementary-alignment join signal: where does the clipped part map?
            # A clipped read whose supplementary alignment lands on ONE other contig
            # is evidence that the two contigs are joined at this position.
            if sa_partners:
                partners = {rn for rn in sa_partners if rn and rn != self.this_contig}
                if partners:
                    self._n_clip_with_sa += 1
                    for rn in partners:
                        self._sa_targets[rn] += 1

    def sa_join_stats(self) -> Tuple[str, int, int]:
        """
        Dominant supplementary-alignment partner for clipped reads.
        Returns (partner_contig, support_reads, n_clip_reads).
        """
        if not self._sa_targets:
            return 'N/A', 0, self._n_clip_reads
        partner, count = max(self._sa_targets.items(), key=lambda kv: kv[1])
        return partner, count, self._n_clip_reads

    def clipped_base_entropy(self) -> float:
        """Shannon entropy (bits) of the base composition of clipped sequences."""
        if not self._clipped_bases:
            return float('nan')
        combined = ''.join(self._clipped_bases).upper()
        n = len(combined)
        if n == 0:
            return float('nan')
        counts = {b: combined.count(b) for b in 'ACGTN'}
        entropy = 0.0
        for b, c in counts.items():
            if c > 0:
                p = c / n
                entropy -= p * math.log2(p)
        return entropy

    def to_dict(self, contig: str, contig_length: int,
                contig_mean_depth: float) -> Dict:
        depth_ratio = (self.n_depth / contig_mean_depth
                       if contig_mean_depth > 0 else float('nan'))
        disc_frac = (self.n_discordant / self.n_paired
                     if self.n_paired > 0 else float('nan'))
        large_ins_frac = (self.n_large_insert / self.n_paired
                          if self.n_paired > 0 and self.n_proper_pair > 0 else float('nan'))
        strand_bias = float('nan')
        if self.n_forward + self.n_reverse > 0:
            strand_bias = self.n_forward / (self.n_forward + self.n_reverse)

        near_end = (
            self.clip_pos < 500
            or (contig_length > 0 and self.clip_pos > contig_length - 500)
        )

        sa_partner, sa_support, n_clip = self.sa_join_stats()
        sa_frac = (sa_support / n_clip) if n_clip > 0 else float('nan')

        return {
            'contig': contig,
            'clip_pos': self.clip_pos,
            'contig_length': contig_length,
            'local_depth': self.n_depth,
            'primary_reads_in_window': self.n_primary,
            'contig_mean_depth': round(contig_mean_depth, 2),
            'depth_ratio': round(depth_ratio, 3),
            'discordant_fraction': round(disc_frac, 3) if not math.isnan(disc_frac) else 'N/A',
            'large_insert_fraction': round(large_ins_frac, 3) if not math.isnan(large_ins_frac) else 'N/A',
            'strand_fwd_fraction': round(strand_bias, 3) if not math.isnan(strand_bias) else 'N/A',
            'clipped_base_entropy': round(self.clipped_base_entropy(), 3)
                                    if not math.isnan(self.clipped_base_entropy()) else 'N/A',
            'near_contig_end': near_end,
            'n_clipped_reads': n_clip,
            'sa_partner_contig': sa_partner,
            'sa_partner_support': sa_support,
            'sa_partner_fraction': round(sa_frac, 3) if not math.isnan(sa_frac) else 'N/A',
        }


def _extract_clipped_bases(cigar: str, seq: str,
                            read_pos: int, clip_pos: int,
                            tol: int = 10) -> str:
    """
    If the read has a soft-clip whose boundary is within `tol` bp of `clip_pos`,
    return the clipped bases; otherwise return ''.
    Uses a simple CIGAR parser — handles S, M, D, I, N, H operations.
    """
    import re
    ops = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
    if not ops:
        return ''

    seq_pos = 0       # position in the read sequence
    ref_pos = read_pos

    for length_s, op in ops:
        length = int(length_s)
        if op == 'S':
            # Soft-clip: read_pos before or after alignment
            # The boundary on the reference side is ref_pos
            if abs(ref_pos - clip_pos) <= tol:
                return seq[seq_pos:seq_pos + length]
            seq_pos += length
        elif op in ('M', '=', 'X'):
            seq_pos += length
            ref_pos += length
        elif op == 'D' or op == 'N':
            ref_pos += length
        elif op == 'I':
            seq_pos += length
        elif op == 'H':
            pass   # hard-clip: bases not in seq

    return ''


def collect_position_stats(
    bam_file: Path, clipping_data: Dict[str, List], window: int, min_mapq: int,
    contig_mean_depths: Dict[str, float], contig_lengths: Dict[str, int],
    insert_mean: float, insert_std: float,
) -> Dict[Tuple[str, int], Dict]:
    """Mean aligned-base depth in each clipped position's local window.

    CIGAR reference coordinates are zero-based. M, = and X contribute bases;
    insertions, deletions and skipped reference bases do not. Primary plus
    secondary depth matches samtools coverage --ff 2052. Pair evidence uses
    paired primary records only. Long reads starting outside the window count.
    """
    if window < 0:
        raise ValueError('window must be non-negative')
    positions = {c: sorted({r[0] for r in rs}) for c, rs in clipping_data.items()}
    accs = {(c, p): _PosAccumulator(p, c) for c, ps in positions.items() for p in ps}
    cmd = ['samtools', 'view', '-F', '2052', '-q', str(min_mapq), str(bam_file)]
    with tempfile.TemporaryFile(mode='w+t') as errors:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=errors, text=True)
        try:
            for line in proc.stdout:
                fields = line.rstrip('\n').split('\t')
                if len(fields) < 11 or fields[2] not in positions:
                    continue
                flag, contig = int(fields[1]), fields[2]
                if flag & 2052:
                    continue
                start = int(fields[3]) - 1
                cursor = start
                blocks = []
                for n, op in re.findall(r'(\d+)([MIDNSHP=X])', fields[5]):
                    n = int(n)
                    if op in 'M=X':
                        blocks.append((cursor, cursor + n))
                        cursor += n
                    elif op in 'DN':
                        cursor += n
                if not blocks:
                    continue
                ps = positions[contig]
                lo = bisect.bisect_left(ps, start - window)
                hi = bisect.bisect_right(ps, cursor + window)
                partners = None
                for opt in fields[11:]:
                    if opt.startswith('SA:Z:'):
                        partners = [entry.split(',', 1)[0] for entry in opt[5:].split(';') if entry]
                        break
                for pos in ps[lo:hi]:
                    length = contig_lengths.get(contig, 0)
                    left, right = max(0, pos - window), pos + window + 1
                    if length:
                        right = min(right, length)
                    overlap = sum(max(0, min(b, right) - max(a, left)) for a, b in blocks)
                    if overlap and right > left:
                        accs[(contig, pos)].add_read(
                            flag, int(fields[8]), fields[9], fields[5], start,
                            insert_mean, insert_std, partners,
                            depth_contribution=overlap / (right - left))
        finally:
            proc.stdout.close()
            proc.wait()
        if proc.returncode:
            errors.seek(0)
            raise RuntimeError(f'samtools view failed (exit {proc.returncode}): {errors.read(300)}')
    return {(c, p): a.to_dict(c, contig_lengths.get(c, 0), contig_mean_depths.get(c, 0))
            for (c, p), a in accs.items()}
