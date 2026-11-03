"""Draw the README figures (static SVG, light + dark) from results/.

    uv run figures.py            # reads results/, writes figures/*.svg
"""

import csv
import gzip
import re
from collections import Counter
from pathlib import Path

RESULTS, OUT = Path("results"), Path("figures")
W = 760
# SPRY2 worked example: CAT CDS blocks on HG02583 hap2 JBHIIV010000008.1 (minus strand), from the CAT GFF3
EXAMPLE = {"assembly": "HG02583_hap2_hprc_r2_v1.0.1", "transcript": "HG02583_ha_T0328322",
           "blocks": [(69857282, 69857938), (69857940, 69858230)]}

LIGHT = {"bg": "#fcfcfb", "ink": "#0b0b0b", "ink2": "#52514e", "muted": "#898781", "grid": "#e1e0d9",
         "axis": "#c3c2b7", "s1": "#2a78d6", "s2": "#eb6834", "na": "#b9b8b0"}
DARK = {"bg": "#1a1a19", "ink": "#f0efec", "ink2": "#c3c2b7", "muted": "#898781", "grid": "#2c2c2a",
        "axis": "#383835", "s1": "#3987e5", "s2": "#d95926", "na": "#55544f"}


def palette_css(c: dict) -> str:
    return (f".bg,.surf{{fill:{c['bg']}}}text{{fill:{c['ink']}}}.sub{{fill:{c['ink2']}}}.m{{fill:{c['muted']}}}"
            f".s1{{fill:{c['s1']}}}.s2{{fill:{c['s2']}}}.na{{fill:{c['na']}}}.on{{fill:#ffffff}}"
            f".outline{{fill:none;stroke:{c['muted']};stroke-width:1.5}}.grid{{stroke:{c['grid']}}}.axis{{stroke:{c['axis']}}}.stop{{stroke:{c['ink']}}}")


STYLE = ("<style>text{font-family:system-ui,-apple-system,'Segoe UI',Helvetica,Arial,sans-serif;font-size:13px}"
         ".t{font-size:16px;font-weight:600}.m{font-size:12px}.on{font-weight:600;font-size:12px}"
         ".grid,.axis{stroke-width:1}.stop{stroke-width:2}" + palette_css(LIGHT)
         + "@media (prefers-color-scheme:dark){" + palette_css(DARK) + "}</style>")


def esc(s) -> str:
    return str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def text(x, y, s, cls="", anchor="start") -> str:
    return f'<text x="{x:.1f}" y="{y:.1f}" class="{cls}" text-anchor="{anchor}">{esc(s)}</text>'


def rect(x, y, w, h, cls, rx=2) -> str:
    return f'<rect x="{x:.1f}" y="{y:.1f}" width="{max(w, 0):.1f}" height="{h:.1f}" rx="{rx}" class="{cls}"/>'


def svg(name: str, h: int, body: list[str], title: str, desc: str) -> None:
    OUT.mkdir(exist_ok=True)
    (OUT / name).write_text(
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {h}" width="{W}" height="{h}" role="img">'
        f"<title>{esc(title)}</title><desc>{esc(desc)}</desc>{STYLE}"
        f'<rect width="{W}" height="{h}" class="bg"/>' + "".join(body) + "</svg>\n")


def legend(x, y, items) -> list[str]:
    out = []
    for cls, label in items:
        out += [rect(x, y - 10, 12, 12, cls), text(x + 18, y, label, "sub")]
        x += 30 + 7 * len(label)
    return out


def load() -> list[dict]:
    with gzip.open(RESULTS / "r2_cat_cds_gaps.tsv.gz", "rt", newline="") as fh:
        return [r for r in csv.DictReader(fh, delimiter="\t") if r["gaps"]]


# --- 1. mechanism: the SPRY2 example -----------------------------------------------------------

def fig_mechanism(rows):
    r = next(r for r in rows if r["assembly"] == EXAMPLE["assembly"] and r["transcript_id"] == EXAMPLE["transcript"])
    (s1, e1), (s2, e2) = EXAMPLE["blocks"]
    first = e2 - s2 + 1  # minus strand: the higher block comes first in transcript order
    cat, rep = int(r["cat_cds_len"]), int(r["repaired_cds_len"])
    stop_codon, codons = map(float, re.search(r"codon(\d+)/([\d.]+)", r["repaired_orf"]).groups())
    stop_bp = int(stop_codon) * 3
    x0, span, gapw = 190, 430, 8  # the 1-bp gap is drawn 8 px wide so it is visible

    def px(bp):  # transcript-order base -> x; bases after the gap shift by the drawn gap width
        return x0 + bp / rep * (span - gapw) + (gapw if bp > first else 0)

    b = [text(24, 30, f"How CAT drops an insertion: {r['gene_name']} in HG02583 haplotype 2", "t"),
         text(24, 50, f"Single-exon CDS, so the gap cannot be a projected intron. Coding strand 5′→3′; "
                      f"{r['net_gap']}-bp gap drawn wider.", "sub")]
    y = 92
    b += [text(24, y + 4, "Assembly sequence", ""), text(24, y + 20, f"{rep} bp incl. inserted {r['gap_bases']}", "m"),
          rect(x0, y - 9, span, 18, "outline"), rect(px(first) - 1, y - 14, gapw + 2, 28, "s2", rx=1),
          text(px(first) + gapw / 2, y - 20, f"inserted {r['gap_bases']}", "m", "middle")]
    y += 62
    b += [text(24, y + 4, "CAT CDS", ""), text(24, y + 20, f"{cat} bp · proper_orf={r['proper_orf']}", "m"),
          text(24, y + 34, f"frameshift={r['frameshift']}", "m"),
          rect(x0, y - 9, px(first) - x0, 18, "s1"), rect(px(first) + gapw, y - 9, x0 + span - px(first) - gapw, 18, "s1"),
          text(x0 + span + 8, y + 4, "reads as intact", "sub")]
    y += 70
    b += [text(24, y + 4, "Repaired CDS", ""), text(24, y + 20, f"{rep} bp, gap bases put back", "m"),
          rect(x0, y - 9, px(first) + gapw - x0, 18, "s1"),
          rect(px(first) + gapw + 2, y - 9, px(stop_bp) - px(first) - gapw - 2, 18, "s2"),
          rect(px(stop_bp) + 2, y - 4, x0 + span - px(stop_bp) - 2, 8, "na"),
          f'<line x1="{px(stop_bp):.1f}" y1="{y - 16}" x2="{px(stop_bp):.1f}" y2="{y + 16}" class="stop"/>',
          text(px(stop_bp), y + 32, f"stop at codon {int(stop_codon)} of {int(codons)}", "", "middle"),
          text((px(first) + px(stop_bp)) / 2, y - 16, "shifted frame", "m", "middle"),
          text(x0 + span + 8, y + 4, "truncated", "sub")]
    b += legend(x0, y + 66, [("s1", "reference frame"), ("s2", "inserted base / shifted frame"), ("na", "not translated")])
    svg("mechanism.svg", y + 84, b, f"How CAT drops an insertion ({r['gene_name']})",
        f"CAT splits the {r['gene_name']} CDS around a {r['net_gap']}-bp gap that holds an inserted base; "
        f"with the base restored the CDS has a premature stop at codon {int(stop_codon)}.")


# --- 2. CAT's flags ----------------------------------------------------------------------------

def fig_flags(rows):
    fc = [r for r in rows if r["frame_changing"] == "True"]
    n = len(fc)
    missed = sum(r["proper_orf"] == "True" and r["frameshift"] != "True" for r in fc)
    bars = [  # (label, [(cls, count, legend role)])
        ("proper_orf", [("s2", sum(r["proper_orf"] == "True" for r in fc)), ("s1", sum(r["proper_orf"] == "False" for r in fc)),
                        ("na", sum(r["proper_orf"] not in ("True", "False") for r in fc))]),
        ("frameshift", [("s2", sum(r["frameshift"] == "False" for r in fc)), ("s1", sum(r["frameshift"] == "True" for r in fc)),
                        ("na", sum(r["frameshift"] not in ("True", "False") for r in fc))]),
        ("CAT CDS as an ORF", [("s2", sum(r["cat_orf"] == "intact" for r in fc)), ("s1", sum(r["cat_orf"] != "intact" for r in fc))]),
        ("repaired CDS as an ORF", [("s2", sum(r["repaired_orf"] == "intact" for r in fc)),
                                    ("s1", sum(r["repaired_orf"] != "intact" for r in fc))]),
    ]
    x0, span = 200, 520
    b = [text(24, 30, f"CAT's flags pass {100 * missed / n:.0f}% of transcripts with a dropped frame-changing insertion", "t"),
         text(24, 50, f"{n:,} CAT transcripts with a frame-changing 1–29 bp CDS gap, all "
                      f"{len({r['assembly'] for r in rows})} HPRC R2 haplotypes", "sub")]
    b += legend(x0, 80, [("s2", "looks fine"), ("s1", "flags a problem"), ("na", "no value (N/A)")])
    y = 112
    for i, (label, segs) in enumerate(bars):
        if i == 3:
            y += 16
            b.append(f'<line x1="24" y1="{y - 22}" x2="{x0 + span}" y2="{y - 22}" class="grid"/>')
        b.append(text(24, y + 5, label, ""))
        x = x0
        for cls, k in segs:
            w = span * k / n
            b.append(rect(x, y - 11, w - 2, 22, cls))
            if w > 30:
                b.append(text(x + 6, y + 4, f"{100 * k / n:.0f}%", "on" if cls != "na" else "", "start"))
            x += w
        y += 40
    b.append(text(24, y + 4, "Flags: CAT's GFF3 tags. ORF rows: this repo's check (one start-to-stop frame, no internal stop).", "m"))
    svg("flags.svg", y + 20, b, "CAT flags on transcripts with a dropped insertion",
        f"Among {n} CAT transcripts with a frame-changing CDS gap, proper_orf is True for most and frameshift is "
        f"often False or N/A; {missed} pass both flags.")


# --- 3. events per haplotype -------------------------------------------------------------------

def fig_per_haplotype(rows):
    events = {(r["assembly"], r["seqid"], r["gaps"]) for r in rows if r["frame_changing"] == "True"}
    per = Counter(a for a, _, _ in events)
    vals = sorted(per.values())
    lo, hi, bw = (min(vals) // 5) * 5, (max(vals) // 5 + 1) * 5, 5
    bins = Counter((v - lo) // bw for v in vals)
    nb = (hi - lo) // bw
    top = max(bins.values())
    ymax = ((top // 20) + 1) * 20
    x0, y0, pw, ph = 70, 250, 640, 170
    b = [text(24, 30, f"Every haplotype is affected: {vals[len(vals) // 2]} dropped frame-changing insertions per "
                      f"haplotype (median)", "t"),
         text(24, 50, f"Distinct frame-changing CDS gap loci per haplotype, {len(vals)} HPRC R2 haplotypes "
                      f"(range {vals[0]}–{vals[-1]})", "sub")]
    for t in range(0, ymax + 1, 20):
        y = y0 - ph * t / ymax
        b += [f'<line x1="{x0}" y1="{y:.1f}" x2="{x0 + pw}" y2="{y:.1f}" class="grid"/>', text(x0 - 8, y + 4, t, "m", "end")]
    bwpx = pw / nb
    for i in range(nb):
        k = bins.get(i, 0)
        h = ph * k / ymax
        if k:
            b.append(f'<path d="M{x0 + i * bwpx + 1:.1f},{y0} v{-h + 4:.1f} q0,-4 4,-4 h{bwpx - 10:.1f} q4,0 4,4 '
                     f'v{h - 4:.1f} z" class="s1"><title>{lo + i * bw}–{lo + (i + 1) * bw - 1}: {k} haplotypes</title></path>')
        if i % 2 == 0:
            b.append(text(x0 + i * bwpx, y0 + 18, lo + i * bw, "m", "middle"))
    b += [f'<line x1="{x0}" y1="{y0}" x2="{x0 + pw}" y2="{y0}" class="axis"/>',
          text(x0 + pw / 2, y0 + 40, "frame-changing CDS gaps per haplotype", "sub", "middle"),
          text(16, y0 - ph / 2, "haplotypes", "sub", "middle").replace("<text", f'<text transform="rotate(-90 16 {y0 - ph / 2})"')]
    svg("per_haplotype.svg", y0 + 56, b, "Dropped insertions per haplotype",
        f"Histogram of frame-changing CDS gaps per haplotype; median {vals[len(vals) // 2]}, range {vals[0]}-{vals[-1]}.")


# --- 4. gap lengths ----------------------------------------------------------------------------

def fig_gap_lengths(rows):
    gaps = {(r["assembly"], r["seqid"], g) for r in rows for g in r["gaps"].split(",")}
    lengths = Counter(int(g.split("-")[1]) - int(g.split("-")[0]) + 1 for _, _, g in gaps)
    import math
    x0, y0, pw, ph = 70, 285, 640, 180
    decades = range(0, 6)  # 1 .. 100,000
    b = [text(24, 30, "Most dropped gaps are 1–2 bp; nearly all change the frame", "t"),
         text(24, 50, f"Distinct CDS gaps (haplotype × locus) by length, log scale: 1 bp {lengths[1]:,} · "
                      f"2 bp {lengths[2]:,} · all {sum(lengths.values()):,}", "sub")]
    b += legend(x0, 78, [("s1", "frame-changing (length not a multiple of 3)"), ("na", "in-frame (multiple of 3)")])
    for d in decades:
        y = y0 - ph * d / 5
        b += [f'<line x1="{x0}" y1="{y:.1f}" x2="{x0 + pw}" y2="{y:.1f}" class="grid"/>',
              text(x0 - 8, y + 4, f"{10 ** d:,}", "m", "end")]
    bw = pw / 29
    for L in range(1, 30):
        k = lengths.get(L, 0)
        if not k:
            continue
        h = ph * math.log10(k) / 5 if k > 1 else 2
        x = x0 + (L - 1) * bw + 2
        b.append(f'<rect x="{x:.1f}" y="{y0 - h:.1f}" width="{bw - 4:.1f}" height="{h:.1f}" rx="2" '
                 f'class="{"na" if L % 3 == 0 else "s1"}"><title>{L} bp: {k:,} gaps</title></rect>')
    for L in (1, 5, 10, 15, 20, 25, 29):
        b.append(text(x0 + (L - 0.5) * bw, y0 + 18, L, "m", "middle"))
    b += [f'<line x1="{x0}" y1="{y0}" x2="{x0 + pw}" y2="{y0}" class="axis"/>',
          text(x0 + pw / 2, y0 + 40, "gap length (bp)", "sub", "middle"),
          text(x0 + (24.5) * bw, y0 - ph * math.log10(lengths[25]) / 5 - 8, "PITRM1", "m", "middle")]
    svg("gap_lengths.svg", y0 + 56, b, "Dropped CDS gap lengths",
        "Distinct CDS gaps by length on a log scale; 1 and 2 bp dominate and multiples of 3 are rare.")


# --- 5. HiFi spot-check ------------------------------------------------------------------------

def fig_hifi():
    with open(RESULTS / "hifi_spotcheck.tsv", newline="") as fh:
        checks = list(csv.DictReader(fh, delimiter="\t"))
    order = [("supported", "s1", "insertion in ≥20% of reads"), ("ambiguous", "na", "5–20% / low coverage"),
             ("absent", "s2", "insertion in ≤5% of reads")]
    x0, span = 200, 520
    b = [text(24, 30, "HiFi reads confirm the dropped bases", "t"),
         text(24, 50, "Random frame-changing events, one per gene: HPRC R2 HiFi reads at the GRCh38 position", "sub")]
    b += legend(x0, 80, [(c, lab) for _, c, lab in order])
    y = 116
    for role, label in (("carrier", "haplotype with the gap"), ("control", "sample without it")):
        rs = [c for c in checks if c["role"] == role]
        st = Counter("ambiguous" if c["status"] in ("ambiguous", "low_coverage") or c["status"].startswith("error")
                     else c["status"] for c in rs)
        b += [text(24, y + 1, f"{role}s ({len(rs)})", ""), text(24, y + 16, label, "m")]
        x = x0
        for status, cls, _ in order:
            w = span * st[status] / len(rs)
            if w:
                b.append(rect(x, y - 11, w - 2, 22, cls))
                if w > 10:
                    b.append(text(x + 4, y + 4, st[status], "on" if cls != "na" else ""))
            x += w
        y += 44
    b.append(text(24, y, "Reads: HPRC R2 HiFi-to-GRCh38 alignments (the assemblies' own reads, so a consensus check).", "m"))
    svg("hifi.svg", y + 16, b, "HiFi read support for dropped insertions",
        "Most carriers have HiFi read support for the dropped insertion; most controls do not.")


if __name__ == "__main__":
    rows = load()
    fig_mechanism(rows)
    fig_flags(rows)
    fig_per_haplotype(rows)
    fig_gap_lengths(rows)
    fig_hifi()
    print("wrote", ", ".join(sorted(p.name for p in OUT.glob("*.svg"))))
