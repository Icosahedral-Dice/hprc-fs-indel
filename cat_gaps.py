"""Find and repair frameshifts hidden by assembly bases that HPRC Release 2 CAT annotations leave out of the CDS.

CAT projects a GRCh38 transcript onto a haplotype assembly and keeps the reference frame by splitting
the CDS into two blocks around a short (1-29 bp) gap. The gap holds the bases of a frameshift
insertion, or bases next to a frameshift deletion that make the net change a multiple of 3. The
annotated CDS omits those bases, looks intact, and often keeps proper_orf=True / frameshift=False.

Subcommands:
  download  fetch CAT GFF3 (and optionally the assembly FASTA) for HPRC R2 samples, plus the GRCh38 GTF
  scan      list CAT transcripts with short CDS gaps; with --fasta, put the gap bases back and
            re-check the ORF; with --gtf, report the GRCh38 source intron check and position
  hifi      count HPRC R2 HiFi reads (aligned to GRCh38) carrying an indel at a position
"""

import argparse
import csv
import gzip
import re
import sys
import urllib.request
from collections import defaultdict
from pathlib import Path

S3 = "https://s3-us-west-2.amazonaws.com/human-pangenomics/"
INDEX_URL = ("https://raw.githubusercontent.com/human-pangenomics/hprc_intermediate_assembly/main/"
             "data_tables/assemblies_release2_v1.0.index.csv")
HIFI_BAM = (S3 + "submissions/40399FDD-59DE-43D1-B3A3-DFF0C6E64FAC--YALE_VARIANT_CALLS_R2/samples/{sample}/"
            "alignments/GRCh38_no_alt/hifi_reads/{sample}.GRCh38_no_alt.sorted.bam")
# GRCh38 source models: CAT v1.3 projects GENCODE v48 (= Ensembl 114); every source transcript ID matches.
ENSEMBL_GTF = "https://ftp.ensembl.org/pub/release-114/gtf/homo_sapiens/Homo_sapiens.GRCh38.114.gtf.gz"
GAP_MAX = 30  # bp; a CDS gap shorter than this is treated as dropped assembly sequence, not an intron
STOPS = {"TAA", "TAG", "TGA"}
COLUMNS = ["transcript_id", "gene_name", "source_transcript", "alignment_mode", "seqid", "strand",
           "cds_blocks", "ref_checked", "ref_intron_gaps", "gaps", "net_gap", "frame_changing", "frameshift",
           "proper_orf", "ref_cds_len", "cat_cds_len", "repaired_cds_len", "gap_bases", "cat_orf", "repaired_orf", "grch38_pos"]


# --- download -----------------------------------------------------------------------------------

def assembly_rows(index: Path) -> list[dict]:
    """HPRC R2 sample haplotypes that have a CAT annotation (drops GRCh38, CHM13, HG002)."""
    with open(index, newline="") as fh:
        return [r for r in csv.DictReader(fh) if r["haplotype"] in ("1", "2") and r["sample_id"] != "HG002"]


def cat_url(row: dict) -> str:
    """CAT GFF3 sits next to the assembly in .../assemblies/release2/annotation/cat/ (HG06807, whose
    assembly is elsewhere, uses working/HPRC/)."""
    url = row["assembly"].replace("s3://human-pangenomics/", S3)
    base = (url.split("/assemblies/release2/")[0] if "/assemblies/release2/" in url
            else f"{S3}working/HPRC/{row['sample_id']}")
    return f"{base}/assemblies/release2/annotation/cat/{row['assembly_name']}_cat_v1.3.gff3.gz"


def fetch(url: str, dest: Path) -> None:
    if dest.exists():
        return
    print(f"{url} -> {dest}", file=sys.stderr)
    tmp = dest.with_name(dest.name + ".part")
    for attempt in range(3):
        try:
            urllib.request.urlretrieve(url, tmp)
            break
        except OSError:  # S3 drops a connection now and then
            if attempt == 2:
                raise
    tmp.rename(dest)


def download(args) -> None:
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    index = out / "assemblies_release2_v1.0.index.csv"
    fetch(INDEX_URL, index)
    fetch(ENSEMBL_GTF, out / ENSEMBL_GTF.rsplit("/", 1)[1])
    rows = [r for r in assembly_rows(index) if args.all or r["sample_id"] in args.samples]
    missing = set(args.samples) - {r["sample_id"] for r in rows}
    if missing:
        sys.exit(f"not in the HPRC R2 index: {' '.join(sorted(missing))}")
    for r in rows:
        fetch(cat_url(r), out / cat_url(r).rsplit("/", 1)[1])
        if args.fasta:
            for key in ("assembly", "assembly_fai", "assembly_gzi"):
                url = S3 + r[key].removeprefix("s3://human-pangenomics/")
                fetch(url, out / url.rsplit("/", 1)[1])


# --- scan ---------------------------------------------------------------------------------------

def attrs(field: str) -> dict:
    return dict(kv.split("=", 1) for kv in field.split(";") if "=" in kv)


def read_cat(gff3: str) -> tuple[dict, dict]:
    """(transcript attributes, CDS blocks [(start, end)] per transcript) from a CAT GFF3."""
    meta, blocks = {}, defaultdict(list)
    with (gzip.open if gff3.endswith(".gz") else open)(gff3, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            if f[2] == "transcript":
                a = attrs(f[8])
                meta[a["ID"]] = {"seqid": f[0], "strand": f[6], **a}
            elif f[2] == "CDS":
                blocks[attrs(f[8])["Parent"]].append((int(f[3]), int(f[4])))
    return meta, blocks


def short_gaps(blocks: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """(start, end) of every gap shorter than GAP_MAX between genomically consecutive CDS blocks."""
    return [(a[1] + 1, b[0] - 1) for a, b in zip(blocks, blocks[1:]) if 0 < b[0] - a[1] - 1 < GAP_MAX]


def revcomp(s: str) -> str:
    return s.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def cds_sequences(fasta, seqid: str, strand: str, blocks: list[tuple[int, int]],
                  fill: list[tuple[int, int]]) -> tuple[str, str, list[str]]:
    """(CAT CDS, repaired CDS with the `fill` gaps put back, their bases) in transcript orientation."""
    cat = [fasta.fetch(seqid, s - 1, e).upper() for s, e in blocks]
    gap_bases = {g: fasta.fetch(seqid, g[0] - 1, g[1]).upper() for g in fill}
    repaired = "".join(seq + gap_bases.get((e + 1, nxt[0] - 1), "")
                       for seq, (s, e), nxt in zip(cat, blocks, blocks[1:] + [(0, 0)]))
    cat, bases = "".join(cat), [gap_bases[g] for g in sorted(fill)]
    if strand == "-":
        cat, repaired, bases = revcomp(cat), revcomp(repaired), [revcomp(b) for b in bases[::-1]]
    return cat, repaired, bases


def orf_status(cds: str) -> str:
    """'intact', or the first reason the CDS is not a single complete ORF."""
    codons = [cds[i:i + 3] for i in range(0, len(cds) - len(cds) % 3, 3)]
    for i, c in enumerate(codons[:-1]):
        if c in STOPS:
            return f"premature_stop@codon{i + 1}/{len(cds) / 3:.1f}"
    if len(cds) % 3:
        return f"length_not_multiple_of_3({len(cds)})"
    if not codons or codons[-1] not in STOPS:
        return "no_terminal_stop"
    return "intact"


def read_gtf(gtf: str) -> dict:
    """GRCh38 CDS blocks per unversioned transcript ID: (chrom, strand, [(start, end)]). Ensembl CDS
    features exclude the stop codon, CAT's include it, so stop_codon features are added."""
    out = {}
    with gzip.open(gtf, "rt") if gtf.endswith(".gz") else open(gtf) as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) < 9 or f[2] not in ("CDS", "stop_codon"):
                continue
            tid = re.search(r'transcript_id "([^".]+)', f[8]).group(1)
            out.setdefault(tid, (f[0], f[6], []))[2].append((int(f[3]), int(f[4])))
    for _, _, b in out.values():
        b.sort()
    return out


def gap_offsets(blocks: list[tuple[int, int]], strand: str, gaps: list[tuple[int, int]]) -> list[int]:
    """Transcript-order CDS offset of each gap: the number of CDS bases before it."""
    return [sum(e - s + 1 for s, e in blocks if (e < gs) == (strand == "+")) for gs, _ in gaps]


def ref_position(ref: tuple, offset: int) -> str:
    """GRCh38 position of the CDS base at 0-based transcript offset (assumes no upstream indel)."""
    chrom, strand, blocks = ref
    for s, e in (blocks if strand == "+" else blocks[::-1]):
        if offset < e - s + 1:
            pos = s + offset if strand == "+" else e - offset
            return f"{chrom if chrom.startswith('chr') else 'chr' + chrom}:{pos}"
        offset -= e - s + 1
    return "NA"


def ref_explains(ref: tuple, length: int, offset: int) -> bool:
    """True if the GRCh38 source model has a CDS intron of this length near this CDS offset, i.e.
    the gap is a projected reference feature (programmed frameshift, reference-error fix)."""
    chrom, strand, blocks = ref
    gaps = short_gaps(blocks)
    return any(e - s + 1 == length and abs(o - offset) <= GAP_MAX
               for (s, e), o in zip(gaps, gap_offsets(blocks, strand, gaps)))


def scan_rows(meta: dict, blocks: dict, fasta=None, ref: dict | None = None):
    names = {n.rsplit("#", 1)[-1]: n for n in fasta.references} if fasta else {}  # CAT seqid -> PanSN
    for tid, b in blocks.items():
        b.sort()
        gaps = short_gaps(b)
        if not gaps or tid not in meta:
            continue
        m = meta[tid]
        offsets = gap_offsets(b, m["strand"], gaps)
        r = ref.get(m.get("source_transcript", "").split(".")[0]) if ref is not None else None
        explained = [bool(r) and ref_explains(r, e - s + 1, o) for (s, e), o in zip(gaps, offsets)]
        dropped = [g for g, x in zip(gaps, explained) if not x]
        net = sum(e - s + 1 for s, e in dropped)
        row = {"transcript_id": tid, "gene_name": m.get("gene_name", ""),
               "source_transcript": m.get("source_transcript", ""), "alignment_mode": m.get("alignment_mode", ""),
               "seqid": m["seqid"], "strand": m["strand"], "cds_blocks": len(b),
               "gaps": ",".join(f"{s}-{e}" for s, e in dropped), "ref_intron_gaps": sum(explained),
               "net_gap": net, "frame_changing": net % 3 != 0,
               "frameshift": m.get("frameshift", ""), "proper_orf": m.get("proper_orf", ""),
               "cat_cds_len": sum(e - s + 1 for s, e in b), "repaired_cds_len": sum(e - s + 1 for s, e in b) + net,
               "ref_checked": ref is not None and bool(r)}
        if r:
            row["ref_cds_len"] = sum(e - s + 1 for s, e in r[2])
            row["grch38_pos"] = ",".join(ref_position(r, o) for o, x in zip(offsets, explained) if not x)
        if fasta:
            cat, repaired, bases = cds_sequences(fasta, names.get(m["seqid"], m["seqid"]), m["strand"], b, dropped)
            row |= {"gap_bases": ",".join(bases), "cat_orf": orf_status(cat), "repaired_orf": orf_status(repaired)}
        yield row


def scan(args) -> None:
    fasta = None
    if args.fasta:
        import pysam
        fasta = pysam.FastaFile(args.fasta)
    ref = read_gtf(args.gtf) if args.gtf else None
    meta, blocks = read_cat(args.gff3)
    out = open(args.out, "w") if args.out else sys.stdout
    w = csv.DictWriter(out, COLUMNS, delimiter="\t", restval="", extrasaction="ignore", lineterminator="\n")
    w.writeheader()
    n = fc = po = fs = 0
    for row in scan_rows(meta, blocks, fasta, ref):
        w.writerow(row)
        n += bool(row["gaps"])
        if row["frame_changing"]:
            fc += 1
            po += row["proper_orf"] == "True"
            fs += row["frameshift"] == "True"
    print(f"{args.gff3}: {len(blocks)} CDS transcripts, {n} with a <{GAP_MAX} bp CDS gap{' not explained by GRCh38' if ref is not None else ''}, "
          f"{fc} frame-changing; of those proper_orf=True {po}, frameshift=True {fs}", file=sys.stderr)


# --- hifi ---------------------------------------------------------------------------------------

def read_indels(read, lo: int, hi: int) -> int:
    """Net indel length (insertions positive) of one alignment within GRCh38 [lo, hi]."""
    net, ref = 0, read.reference_start + 1
    for op, n in read.cigartuples:
        if op == 1 and lo <= ref <= hi:
            net += n
        elif op == 2 and lo <= ref <= hi:
            net -= n
        if op in (0, 2, 3, 7, 8):
            ref += n
    return net


def hifi_support(sample: str, region: str, net: int, window: int = 10) -> tuple[int, int]:
    """(spanning reads, reads with an indel of this net length within +-window) in the sample's
    HPRC R2 HiFi-to-GRCh38 alignment. htslib caches the .bai in the working directory."""
    import pysam
    chrom, pos = region.split(":")
    lo, hi = int(pos) - window, int(pos) + window
    span = carry = 0
    with pysam.AlignmentFile(HIFI_BAM.format(sample=sample)) as bam:
        for read in bam.fetch(chrom, lo - 1, hi):
            if read.is_secondary or read.is_supplementary or read.reference_start + 1 > lo or read.reference_end < hi:
                continue
            span += 1
            carry += read_indels(read, lo, hi) == net
    return span, carry


def hifi(args) -> None:
    span, carry = hifi_support(args.sample, args.region, args.net, args.window)
    share = carry / span if span else 0
    print(f"{args.sample}\t{args.region}\tnet={args.net:+d}\tspanning={span}\tcarrying={carry}\tshare={share:.2f}")


def main(argv=None) -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(required=True)
    d = sub.add_parser("download", help="fetch CAT GFF3 (+ assembly FASTA) for HPRC R2 samples")
    d.add_argument("samples", nargs="*")
    d.add_argument("--all", action="store_true", help="every R2 sample haplotype")
    d.add_argument("--fasta", action="store_true", help="also fetch the bgzipped assembly + .fai/.gzi (~1 GB each)")
    d.add_argument("--out", default="data")
    d.set_defaults(func=download)
    s = sub.add_parser("scan", help="list CAT transcripts with short CDS gaps")
    s.add_argument("gff3")
    s.add_argument("--fasta", help="bgzipped assembly (path or https URL) to restore gap bases")
    s.add_argument("--gtf", help="GRCh38 Ensembl 114 / GENCODE v48 GTF: per-gap source-intron check, GRCh38 position")
    s.add_argument("--out", help="TSV path (default stdout)")
    s.set_defaults(func=scan)
    h = sub.add_parser("hifi", help="HiFi read support for an indel at a GRCh38 position")
    h.add_argument("sample")
    h.add_argument("region", help="chrN:pos (GRCh38)")
    h.add_argument("--net", type=int, required=True, help="net indel length, insertions positive")
    h.add_argument("--window", type=int, default=10)
    h.set_defaults(func=hifi)
    args = p.parse_args(argv)
    if getattr(args, "func", None) is download and not (args.samples or args.all):
        p.error("give sample IDs or --all")
    args.func(args)


if __name__ == "__main__":
    main()
