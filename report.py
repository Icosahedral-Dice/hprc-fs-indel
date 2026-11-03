"""Aggregate per-haplotype `cat_gaps.py scan` outputs and HiFi-check a random sample of events.

    uv run report.py OUT_DIR [--hifi N]

OUT_DIR holds one `{assembly_name}.tsv` per haplotype (scanned with --fasta and --gtf). Writes
results/r2_cat_cds_gaps.tsv.gz, results/summary.md and, with --hifi, results/hifi_spotcheck.tsv.
"""

import argparse
import csv
import gzip
import random
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from cat_gaps import COLUMNS, hifi_support

SUPPORTED, ABSENT, MIN_READS = 0.20, 0.05, 5  # share of spanning reads carrying the indel


def load(out_dir: Path) -> list[dict]:
    rows = []
    for tsv in sorted(out_dir.glob("*.tsv")):
        with open(tsv, newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                rows.append({"assembly": tsv.stem, "sample": tsv.stem.split("_")[0], **r})
    return rows


def pct(n: int, d: int) -> str:
    return f"{n:,} / {d:,} ({100 * n / d:.0f}%)" if d else "0"


def summarize(rows: list[dict], n_hap: int) -> str:
    gap = [r for r in rows if r["gaps"]]
    fc = [r for r in gap if r["frame_changing"] == "True"]
    events = {(r["assembly"], r["seqid"], r["gaps"]) for r in fc}
    genes = {r["gene_name"] for r in fc}
    per_hap = Counter(a for a, _, _ in events)
    lengths = Counter(int(g.split("-")[1]) - int(g.split("-")[0]) + 1 for r in gap for g in r["gaps"].split(","))
    proper = Counter(r["proper_orf"] for r in fc)
    fshift = Counter(r["frameshift"] for r in fc)
    missed = sum(r["proper_orf"] == "True" and r["frameshift"] != "True" for r in fc)
    hidden = sum(r["cat_orf"] == "intact" and r["repaired_orf"] != "intact" for r in fc)
    restored = sum(r["cat_orf"] != "intact" and r["repaired_orf"] == "intact" for r in fc)
    ref_intron = [r for r in rows if int(r["ref_intron_gaps"])]
    unchecked = sum(r["ref_checked"] != "True" for r in gap)
    hap_counts = sorted(per_hap.get(a, 0) for a in {r["assembly"] for r in rows})
    lines = [
        f"Haplotypes scanned: {n_hap}",
        "",
        "| CAT protein-coding transcripts with a 1-29 bp CDS gap not explained by GRCh38 | count |",
        "| --- | ---: |",
        f"| transcripts (all haplotypes) | {len(gap):,} |",
        f"| with a net frame-changing gap | {len(fc):,} |",
        f"| distinct frame-changing gap events (haplotype x locus) | {len(events):,} |",
        f"| genes affected | {len(genes):,} |",
        f"| events per haplotype, median (range) | {hap_counts[len(hap_counts) // 2]} "
        f"({hap_counts[0]}-{hap_counts[-1]}) |",
        "",
        "| Among frame-changing-gap transcripts | |",
        "| --- | ---: |",
        f"| `proper_orf=True` | {pct(proper['True'], len(fc))} |",
        f"| `frameshift=True` | {pct(fshift['True'], len(fc))} |",
        f"| `frameshift=False` | {pct(fshift['False'], len(fc))} |",
        f"| `frameshift=N/A` | {pct(fshift['N/A'], len(fc))} |",
        f"| missed by both flags (`proper_orf=True`, `frameshift` not True) | {pct(missed, len(fc))} |",
        f"| CAT CDS is an intact ORF, repaired CDS is not | {pct(hidden, len(fc))} |",
        f"| CAT CDS is broken, repaired CDS is an intact ORF | {pct(restored, len(fc))} |",
        "",
        "Gap length (bp): " + ", ".join(f"{k}: {v:,}" for k, v in sorted(lengths.items())),
        "",
        f"Short CDS gaps that match a GRCh38 CDS intron (kept as-is): {len(ref_intron):,} transcripts in "
        f"{len({r['gene_name'] for r in ref_intron}):,} genes, e.g. "
        + ", ".join(g for g, _ in Counter(r["gene_name"] for r in ref_intron).most_common(8)),
        f"Gap transcripts whose source transcript is missing from the GTF (unchecked): {unchecked:,}",
    ]
    return "\n".join(lines) + "\n"


def spotcheck(rows: list[dict], n: int, seed: int = 1) -> list[dict]:
    """HiFi support for n random frame-changing insertion events (one per gene) in the carrier,
    and in one random sample whose two haplotypes have no CDS gap in that transcript."""
    rng = random.Random(seed)
    gapped = {(r["sample"], r["source_transcript"]) for r in rows if r["gaps"]}
    samples = sorted({r["sample"] for r in rows})
    pool = [r for r in rows if r["gaps"] and r["frame_changing"] == "True" and "," not in r["gaps"]
            and r["cat_cds_len"] == r["ref_cds_len"] and r["grch38_pos"] not in ("", "NA")
            and r["cat_orf"] == "intact" and r["repaired_orf"] != "intact"]
    rng.shuffle(pool)
    picked, seen = [], set()
    for r in pool:
        if r["gene_name"] not in seen and len(picked) < 2 * n:
            controls = [s for s in samples if (s, r["source_transcript"]) not in gapped]
            if not controls:  # every sample has a gap here
                continue
            seen.add(r["gene_name"])
            picked.append((r, "carrier", r["sample"]))
            picked.append((r, "control", rng.choice(controls)))

    def run(job):
        r, role, sample = job
        try:
            span, carry = hifi_support(sample, r["grch38_pos"], int(r["net_gap"]))
        except (OSError, ValueError) as e:
            return {"role": role, "sample": sample, "gene": r["gene_name"], "status": f"error:{e}"}
        share = carry / span if span else 0
        status = ("low_coverage" if span < MIN_READS else "supported" if share >= SUPPORTED
                  else "absent" if share <= ABSENT else "ambiguous")
        return {"role": role, "sample": sample, "gene": r["gene_name"], "source_transcript": r["source_transcript"],
                "alignment_mode": r["alignment_mode"],
                "grch38_pos": r["grch38_pos"], "net": r["net_gap"], "gap_bases": r["gap_bases"],
                "spanning": span, "carrying": carry, "share": f"{share:.2f}", "status": status}

    with ThreadPoolExecutor(8) as ex:
        return list(ex.map(run, picked))


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("out_dir", type=Path)
    p.add_argument("--hifi", type=int, default=0, help="number of random events to HiFi-check")
    p.add_argument("--results", type=Path, default=Path("results"))
    args = p.parse_args()
    rows = load(args.out_dir)
    args.results.mkdir(exist_ok=True)
    with gzip.open(args.results / "r2_cat_cds_gaps.tsv.gz", "wt", newline="") as fh:
        w = csv.DictWriter(fh, ["assembly", "sample", *COLUMNS], delimiter="\t", extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    text = summarize(rows, len(list(args.out_dir.glob("*.tsv"))))
    if args.hifi:
        checks = spotcheck(rows, args.hifi)
        with open(args.results / "hifi_spotcheck.tsv", "w", newline="") as fh:
            w = csv.DictWriter(fh, ["role", "sample", "gene", "source_transcript", "alignment_mode", "grch38_pos", "net", "gap_bases",
                                    "spanning", "carrying", "share", "status"], delimiter="\t", restval="", lineterminator="\n")
            w.writeheader()
            w.writerows(checks)
        tally = Counter((c["role"], c["status"]) for c in checks)
        text += "\nHiFi spot-check (" + "; ".join(f"{role} {status}: {k}" for (role, status), k in sorted(tally.items())) + ")\n"
    (args.results / "summary.md").write_text(text)
    print(text)


if __name__ == "__main__":
    main()
