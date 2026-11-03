# hprc-fs-ins

**HPRC Release 2 CAT gene annotations leave frame-changing insertions out of the CDS.**

This repository demonstrates the problem on public HPRC Release 2 data and
provides a small Python tool (`cat_gaps.py`) that finds the affected
transcripts, puts the missing assembly bases back, and re-checks the ORF.

## The problem

When a haplotype assembly carries an insertion that shifts the reading frame,
CAT (Comparative Annotation Toolkit) projects the GRCh38 transcript onto the
assembly and keeps the *reference* frame. To do that it splits the CDS into
two blocks with a short gap (1–29 bp) between them, and the inserted bases sit
in that gap. The annotated CDS therefore **does not contain the inserted
bases**:

```
assembly   ...ATG GCC AAG T CAG CTG...       (T = inserted base)
CAT CDS    [ATG GCC AAG] [CAG CTG...]        two CDS blocks, 1-bp gap, no GRCh38 intron here
           -> stays in frame, often exactly the reference length
true CDS   ATG GCC AAG TCA GCT G...          -> frame-shifted, premature stop
```

The CAT CDS stays in frame, often matches the reference length exactly, and
usually keeps `proper_orf=True`. The `frameshift` transcript tag catches only
some of these cases. Anyone who takes the CAT CDS, or these flags, at face
value will call the copy intact: **frameshift insertions are undercounted**,
while frameshift deletions (which do not need a gap) are mostly unaffected.

The same mechanism shows up in two less common forms:

- **Next to a deletion.** CAT sometimes leaves 1–4 assembly bases next to a
  frameshift deletion out of the CDS, so that the combined change is a
  multiple of 3 and the frame is kept. Putting the bases back restores the
  frameshift deletion.
- **Restoring the frame.** Where the GRCh38 model is itself frame-disrupted
  relative to most haplotypes, the dropped insertion can be the one that
  *restores* the frame. GPATCH4 (`ENST00000368232`) is the clearest case: in
  458 of 462 haplotypes the CAT CDS has a premature stop, and with the 2
  dropped bases (`GT`) put back it is an intact ORF.

`cat_gaps.py` handles all three, because it simply puts the gap bases back and
re-reads the ORF.

## Worked example: SPRY2 in HG02583

SPRY2 has a single-exon CDS (948 bp in GRCh38, `ENST00000377104`), so a
gap inside it cannot be a projected intron.

HG02583 haplotype 2 (`HG02583_hap2_hprc_r2_v1.0.1`), contig
`JBHIIV010000008.1`, minus strand. CAT has three models of this locus:

| CAT transcript | source | mode | CDS blocks | CDS length | `frameshift` | `proper_orf` |
| --- | --- | --- | --- | ---: | --- | --- |
| `HG02583_ha_T0328322` | SPRY2-201 `ENST00000377102.5` | transMap_pairwise | 69857282–69857938, 69857940–69858230 | 948 | **False** | **True** |
| `HG02583_ha_T0328323` | SPRY2-202 `ENST00000377104.4` | transMap | 69857282–69857938, 69857940–69858230 | 948 | True | **True** |
| `HG02583_ha_T0328324` | SPRY2-202 `ENST00000377104.4` | liftoff | 69857282–69858230 | 949 | **False** | **True** |

- The two transMap models skip assembly base 69857939 (a `C`; `G` on the
  coding strand). Their CDS is 948 bp, the reference length, and reads as an
  intact ORF.
- With the base put back, the CDS is 949 bp and hits a premature stop at
  codon 132 of 316.
- The liftoff model does include the base, but is still tagged
  `frameshift=False, proper_orf=True` although its CDS length (949) is not a
  multiple of 3.
- Haplotype 1 (`HG02583_hap1_hprc_r2_v1.0.1`) has a single 948-bp CDS block
  with no gap.
- HPRC R2 HiFi reads aligned to GRCh38: 30 of 73 reads spanning
  chr13:80337414 carry a +1 insertion (42 carry none), as expected for a
  heterozygous insertion present only on haplotype 2.

Reproduce (about 290 MB download plus a minute of remote assembly reads):

```bash
uv sync
uv run cat_gaps.py download HG02583            # CAT GFF3 for both haplotypes + Ensembl 114 GTF -> data/
uv run cat_gaps.py scan data/HG02583_hap2_hprc_r2_v1.0.1_cat_v1.3.gff3.gz \
    --fasta https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC/HG02583/assemblies/release2/HG02583_hap2_hprc_r2_v1.0.1.fa.gz \
    --gtf data/Homo_sapiens.GRCh38.114.gtf.gz --out hap2.tsv
grep -P '\tSPRY2\t' hap2.tsv
uv run cat_gaps.py hifi HG02583 chr13:80337414 --net 1
```

This SPRY2 insertion is rare: across all 462 haplotypes the scan finds it
only in HG02583 haplotype 2.

## Genome-wide: all 462 HPRC Release 2 haplotypes

`cat_gaps.py scan` was run on every R2 CAT v1.3 GFF3 (about 113,000
protein-coding CAT transcripts per haplotype), with the assembly and the
Ensembl 114 GTF (`results/summary.md`, full table in
`results/r2_cat_cds_gaps.tsv.gz`):

| CAT transcripts with a 1–29 bp CDS gap not explained by GRCh38 | count |
| --- | ---: |
| transcripts, all haplotypes | 94,123 |
| — with a net frame-changing gap | 90,379 |
| distinct frame-changing gap events (haplotype × locus) | 55,401 |
| genes affected | 1,795 |
| events per haplotype, median (range) | 119 (93–171) |

| Among the 90,379 frame-changing-gap transcripts | |
| --- | ---: |
| `proper_orf=True` | 68,739 (76%) |
| `frameshift=True` | 33,843 (37%) |
| `frameshift=False` / `N/A` | 31,056 (34%) / 25,480 (28%) |
| **missed by both flags** (`proper_orf=True`, `frameshift` not True) | **41,962 (46%)** |
| CAT CDS is an intact ORF, repaired CDS is not | 62,073 (69%) |
| CAT CDS is broken, repaired CDS is an intact ORF | 7,517 (8%) |

Gap lengths: 1 bp 54,271; 2 bp 17,100; 4 bp 6,160; 5 bp 4,140; 25 bp 4,901
(almost all PITRM1, see Limits); the rest spread over 3–29 bp.

Short gaps in another 18,491 CAT transcripts (15 genes, mostly OAZ1, OAZ3, SHANK3,
SLC37A4, OAZ2, TCEAL6, ABO and PEG10) match a CDS "intron" of the GRCh38
source model and are left alone (see "The fix").

**HiFi spot-check** (`results/hifi_spotcheck.tsv`, `report.py --hifi 40`).
Picks random frame-changing events, one per gene, from transcripts whose CAT
CDS is the reference length and reads as intact. For each event it queries
the carrier's HPRC R2 HiFi-to-GRCh38 alignment at the event's GRCh38
position, plus one random sample with no gap in that transcript on either
haplotype. "Supported" means ≥20% of spanning reads carry an indel of the
same net length within ±10 bp; "absent" means ≤5%.

| | supported | absent | ambiguous |
| --- | ---: | ---: | ---: |
| carriers (40 events) | **37** | 3 | 0 |
| controls (40 samples) | 1 | 32 | 7 |

Supported events range from 1-bp insertions (most) to insertions of up to
29 bp (`GGGGC` in ZFYVE19, 16 bp in CIDEA, 29 bp in ASAH1). Of the
three carriers without read support, one is the 25-bp PITRM1 gap discussed
under Limits: no +25 bp indel at that GRCh38 position fits a shifted splice
site rather than an insertion. The other two are a 1-bp gap (ITGB1BP1) and
a 10-bp repeat (`AGGAGGAGGA`, SLC12A5). All 8 controls that are not
"absent" sit at 1–2 bp A/T homopolymer insertions (5–20% of reads), and 7 of
the 8 events are common (in 88–394 haplotypes), which is consistent with HiFi
homopolymer noise or an unflagged carrier. These are the reads the assemblies
were built from, so this checks that the dropped bases are real assembly
sequence. It is not an independent technology.

## The fix

`cat_gaps.py scan` works on one CAT GFF3 at a time:

1. For each CAT transcript, sort its CDS blocks and find every gap shorter
   than 30 bp between consecutive blocks.
2. With `--gtf`, look up the GRCh38 source transcript (CAT v1.3 projects
   GENCODE v48 = Ensembl 114; every source ID matches). A gap is kept as-is if
   the source model has a CDS intron of the same length at about the same CDS
   offset. A few GRCh38 models do contain such tiny "introns", for example
   the programmed ribosomal frameshifts of the antizymes (OAZ1/2/3) and
   PEG10, and 1–2 bp gaps in SLC37A4, SHANK3 and ABO. Those are projected
   reference features, not dropped insertions, and are not repaired. Every
   other short gap is treated as assembly sequence that CAT left out.
3. With `--fasta`, read the CDS blocks and the dropped gaps from the
   bgzipped assembly (local path or HTTPS URL; CAT seqids are matched to the
   PanSN FASTA names), build the CAT CDS and the repaired CDS, and report
   whether each is an intact ORF.
4. With `--gtf`, also report the approximate GRCh38 position of each gap
   (assuming no upstream indel in the copy) so it can be checked against
   HiFi reads (`cat_gaps.py hifi`) or a variant call set.

Output columns:

| column | meaning |
| --- | --- |
| `transcript_id`, `gene_name`, `source_transcript`, `alignment_mode` | from the CAT GFF3 |
| `seqid`, `strand`, `cds_blocks` | CAT model on the assembly |
| `ref_checked` | source transcript found in the GTF |
| `ref_intron_gaps` | short gaps that match a GRCh38 CDS intron (left alone) |
| `gaps`, `net_gap`, `frame_changing` | dropped gaps (assembly coordinates), their total length, `net_gap % 3 != 0` |
| `frameshift`, `proper_orf` | CAT's own tags |
| `ref_cds_len`, `cat_cds_len`, `repaired_cds_len` | CDS length in GRCh38, in CAT, after repair |
| `gap_bases` | the dropped bases, coding strand |
| `cat_orf`, `repaired_orf` | `intact`, or why not (`premature_stop@codonN/M`, `length_not_multiple_of_3`, `no_terminal_stop`) |
| `grch38_pos` | approximate GRCh38 position of each dropped gap |

Downstream use: replace the CAT CDS with the repaired CDS before calling
loss of function, pseudogenes or ORF integrity, and do not rely on
`proper_orf` or `frameshift` to flag these transcripts.

## Files

- `cat_gaps.py`: `download` (CAT GFF3, optional assembly FASTA, GRCh38 GTF from
  the public HPRC bucket and Ensembl), `scan`, `hifi`.
- `report.py`: aggregates per-haplotype scans into `results/` and HiFi-checks
  a random sample of events.
- `test_cat_gaps.py`: synthetic checks (`uv run python test_cat_gaps.py`).
- `results/`: genome-wide output for all HPRC R2 haplotypes.

Full genome-wide run (≈ 51 GB of GFF3; each scan reads its assembly over
HTTPS, so no FASTA download is needed):

```bash
uv run cat_gaps.py download --all --out data
# for each haplotype:
uv run cat_gaps.py scan data/{assembly}_cat_v1.3.gff3.gz --fasta {assembly URL from the index} \
    --gtf data/Homo_sapiens.GRCh38.114.gtf.gz --out out/{assembly}.tsv
uv run report.py out --hifi 40
```

## Limits

- Gaps of 30 bp or more are not restored; they are indistinguishable from
  introns by length alone.
- A short gap is not always a coding insertion. Most gaps are 1–2 bp, and
  many longer ones are repeat-length changes (for example `GCCCC` in GNRH2,
  `CTGCTGCTGCT` in MICA). A few longer gaps carry splice-signal-like
  sequence: the 25-bp PITRM1 gap seen in about two-thirds of haplotypes,
  `CATTCCTCTCTCTTCCCTAG|CGTGC`, ends in an acceptor-like `...TTCCCTAG`, so it
  may reflect a shifted splice site rather than translated sequence. The
  repaired CDS is the assembly sequence between the CAT blocks; it is not a
  splicing prediction.
- Events present in most haplotypes usually mean that the GRCh38 model
  carries the minor allele or is an unusual isoform. The CAT CDS is still
  wrong for that assembly, but "repaired ORF broken" then says more about
  the reference model than about loss of function.
- The GRCh38 position is approximate when the copy has other indels upstream
  of the gap.
- The HiFi check counts reads with an indel of the same net length within
  ±10 bp; it does not compare the inserted sequence. These are the reads the
  assemblies were built from, so it is a consensus check, not an independent
  technology.
- Scans are per CAT transcript. CAT often has several models (transMap,
  liftoff, alternative source transcripts) for one locus, so transcript
  counts exceed locus counts; `report.py` also counts distinct loci.
