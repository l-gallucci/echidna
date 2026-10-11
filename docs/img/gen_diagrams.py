#!/usr/bin/env python3
"""Generate the two MaMISA dataflow diagrams:
  docs/img/misassembly.svg   - Branch A (PacBio HiFi only)
  docs/img/mag_quality.svg   - Branch B (MAG quality + taxonomy)

Node notation:
  filled dot  = step needs an external input file you provide
  ring (empty)= step runs only on outputs of earlier steps
Each incoming/outgoing chip is tagged with its file type.
"""
import html

# (num, title, tool, left_inputs[(type,label)], right_outputs[(type,label)])
BRANCH_A = [
    (1, "process-large-contigs", "mamisa",
     [("FASTA", "assembly"), ("DIR", "misassemblies")],
     [("FASTA", "assembly_for_filtering"), ("DIR", "HQ_extracted")]),
    (2, "anvi-script-find-misassemblies", "anvi'o (external)",
     [("BAM", "HiFi mapping")],
     [("TSV", "*-clipping.txt")]),
    (3, "check-read-chimeras", "mamisa",
     [("BAM", "HiFi mapping"), ("TSV", "kraken2/kaiju")],
     [("TSV", "chimera_read_report + windows")]),
    (4, "classify-clipping", "mamisa",
     [("BAM", "HiFi mapping")],
     [("TSV", "clipping_classified")]),
    (5, "filter-misassemblies", "mamisa",
     [],
     [("FASTA", "assembly_clean")]),
]
BRANCH_B = [
    (1, "binning", "MetaBAT2 / etc (external)",
     [("FASTA", "clean assembly"), ("BAM", "mapping")],
     [("FASTA", "bins/")]),
    (2, "run-checkm2", "CheckM2",
     [],
     [("TSV", "quality_report")]),
    (3, "organize-mags", "mamisa",
     [],
     [("DIR", "HQ / MQ / LQ")]),
    (4, "run-gtdbtk", "GTDB-Tk",
     [("DIR", "GTDB-Tk DB")],
     [("TSV", "gtdbtk.*.summary")]),
    (5, "run-gunc", "GUNC",
     [("DB", "GUNC .dmnd")],
     [("TSV", "maxCSS + contig_assignments")]),
    (6, "run-mimag-rna", "barrnap / tRNAscan-SE",
     [],
     [("TSV", "mimag_rna summary")]),
    (7, "check-chimeras", "mamisa",
     [],
     [("TSV", "chimera_report")]),
    (8, "organize-mags (final)", "mamisa",
     [],
     [("DIR", "final MAGs")]),
]

TYPE_COLOR = {
    "FASTA": "#2d7d46", "BAM": "#a63d8f", "TSV": "#1f6feb",
    "DB": "#b5651d", "DIR": "#6e7781",
}

W = 940
BOX_X, BOX_W, BOX_H = 330, 300, 56
PITCH = 92
LEFT_EDGE = 320
RIGHT_START = 640

def esc(s):
    return html.escape(s, quote=True)

def chip(x, y, typ, label, anchor):
    color = TYPE_COLOR.get(typ, "#6e7781")
    tag_w = 46
    parts = []
    tag_x = x - tag_w if anchor == "end" else x
    parts.append(f'<rect x="{tag_x}" y="{y-10}" width="{tag_w}" height="20" rx="4" fill="{color}"/>')
    parts.append(f'<text x="{tag_x+tag_w/2}" y="{y+4}" text-anchor="middle" font-size="11" '
                 f'fill="#fff" font-weight="600">{esc(typ)}</text>')
    if anchor == "end":
        parts.append(f'<text x="{tag_x-6}" y="{y+4}" text-anchor="end" font-size="11" '
                     f'fill="#444">{esc(label)}</text>')
    else:
        parts.append(f'<text x="{tag_x+tag_w+6}" y="{y+4}" text-anchor="start" font-size="11" '
                     f'fill="#444">{esc(label)}</text>')
    return "\n".join(parts)

def render(steps, title, note=None):
    svg = []
    top = 150 + (26 if note else 0)
    total_h = top + len(steps) * PITCH + 20
    svg.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{total_h}" '
               f'viewBox="0 0 {W} {total_h}" font-family="-apple-system,Segoe UI,Roboto,sans-serif">')
    svg.append(f'<rect width="{W}" height="{total_h}" fill="#ffffff"/>')
    svg.append('<defs><marker id="arr" markerWidth="8" markerHeight="8" refX="6" refY="3" '
               'orient="auto"><path d="M0,0 L6,3 L0,6 Z" fill="#8a8f98"/></marker></defs>')
    svg.append(f'<text x="{W/2}" y="34" text-anchor="middle" font-size="20" font-weight="700" '
               f'fill="#1a1a1a">{esc(title)}</text>')
    # legend
    ly = 64
    svg.append(f'<circle cx="40" cy="{ly}" r="7" fill="#1f6feb" stroke="#1f6feb"/>')
    svg.append(f'<text x="54" y="{ly+4}" font-size="12" fill="#333">needs an external input file you provide</text>')
    svg.append(f'<circle cx="430" cy="{ly}" r="7" fill="#ffffff" stroke="#8a8f98" stroke-width="2"/>')
    svg.append(f'<text x="444" y="{ly+4}" font-size="12" fill="#333">runs only on earlier steps\' output</text>')
    svg.append(f'<text x="40" y="{ly+24}" font-size="11" fill="#666">chips show file type: '
               f'FASTA sequences, BAM alignments, TSV tables, DB diamond .dmnd, DIR a directory</text>')
    if note:
        svg.append(f'<rect x="20" y="{ly+34}" width="{W-40}" height="26" rx="5" '
                   f'fill="#fff4e5" stroke="#e0a458"/>')
        svg.append(f'<text x="34" y="{ly+51}" font-size="12" font-weight="600" '
                   f'fill="#8a5a0a">{esc(note)}</text>')

    for i, (num, title_s, tool, lefts, rights) in enumerate(steps):
        y = top + i * PITCH
        cy = y + BOX_H / 2
        external = len(lefts) > 0
        svg.append(f'<rect x="{BOX_X}" y="{y}" width="{BOX_W}" height="{BOX_H}" rx="8" '
                   f'fill="#f6f8fa" stroke="#c9d1d9"/>')
        if external:
            svg.append(f'<circle cx="{BOX_X+18}" cy="{cy}" r="7" fill="#1f6feb" stroke="#1f6feb"/>')
        else:
            svg.append(f'<circle cx="{BOX_X+18}" cy="{cy}" r="7" fill="#ffffff" '
                       f'stroke="#8a8f98" stroke-width="2"/>')
        svg.append(f'<text x="{BOX_X+34}" y="{cy-4}" font-size="13.5" font-weight="700" '
                   f'fill="#1a1a1a">[{num}] {esc(title_s)}</text>')
        svg.append(f'<text x="{BOX_X+34}" y="{cy+13}" font-size="11" fill="#6e7781">{esc(tool)}</text>')
        n = len(lefts)
        for j, (typ, label) in enumerate(lefts):
            yy = cy if n == 1 else cy - 14 + j * 28
            svg.append(chip(LEFT_EDGE, yy, typ, label, "end"))
            svg.append(f'<line x1="{LEFT_EDGE+2}" y1="{yy}" x2="{BOX_X-2}" y2="{cy}" '
                       f'stroke="#8a8f98" stroke-width="1.3" marker-end="url(#arr)"/>')
        m = len(rights)
        for j, (typ, label) in enumerate(rights):
            yy = cy if m == 1 else cy - 14 + j * 28
            svg.append(f'<line x1="{BOX_X+BOX_W+2}" y1="{cy}" x2="{RIGHT_START-2}" y2="{yy}" '
                       f'stroke="#8a8f98" stroke-width="1.3" marker-end="url(#arr)"/>')
            svg.append(chip(RIGHT_START, yy, typ, label, "start"))
        if i < len(steps) - 1:
            svg.append(f'<line x1="{BOX_X+BOX_W/2}" y1="{y+BOX_H}" x2="{BOX_X+BOX_W/2}" '
                       f'y2="{y+PITCH-2}" stroke="#444" stroke-width="2" marker-end="url(#arr)"/>')
    svg.append('</svg>')
    return "\n".join(svg)

if __name__ == "__main__":
    import os, sys
    out_dir = sys.argv[1] if len(sys.argv) > 1 else "."
    os.makedirs(out_dir, exist_ok=True)
    a = render(BRANCH_A, "MaMISA - Misassembly cleaning (Branch A)",
               note="PacBio HiFi only: the soft-clipping method needs long reads "
                    "mapped to their own assembly (Trigodet et al. 2025).")
    with open(os.path.join(out_dir, "misassembly.svg"), "w") as f:
        f.write(a)
    b = render(BRANCH_B, "MaMISA - MAG quality and taxonomy (Branch B)")
    with open(os.path.join(out_dir, "mag_quality.svg"), "w") as f:
        f.write(b)
    print("wrote misassembly.svg + mag_quality.svg to", out_dir)
