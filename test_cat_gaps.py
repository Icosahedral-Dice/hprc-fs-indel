"""Synthetic check: a CAT model split around a 1-bp insertion is found and repaired on both strands."""

import tempfile

from cat_gaps import read_gtf, revcomp, scan_rows


class Fasta:
    def __init__(self, seqs):
        self.seqs, self.references = seqs, [f"S#1#{k}" for k in seqs]

    def fetch(self, name, start, end):
        return self.seqs[name.rsplit("#", 1)[-1]][start:end]


def test_insertion_repair():
    # Reference ORF ATG AAA CCC GGG TAA; the assembly carries an extra T after ATGAAA, which makes
    # the repaired frame ATG AAA TCC CGG GTA A -> not a multiple of 3 and no terminal stop.
    asm = "NNNN" + "ATGAAA" + "T" + "CCCGGGTAA" + "NNNN"
    blocks = [(5, 10), (12, 20)]  # CAT skips assembly base 11 (the T)
    meta = {"t": {"seqid": "c", "strand": "+", "gene_name": "G", "source_transcript": "ENST1.1",
                  "proper_orf": "True", "frameshift": "False"}}
    ref = {"ENST1": ("7", "+", [(101, 115)])}
    for strand, seq in (("+", asm), ("-", revcomp(asm))):
        n = len(seq)
        b = blocks if strand == "+" else [(n - e + 1, n - s + 1) for s, e in reversed(blocks)]
        meta["t"]["strand"] = strand
        ref["ENST1"] = ("7", strand, [(101, 115)])
        (row,) = scan_rows(meta, {"t": list(b)}, Fasta({"c": seq}), ref)
        assert row["net_gap"] == 1 and row["frame_changing"], row
        assert row["ref_cds_len"] == row["cat_cds_len"] == 15 and row["repaired_cds_len"] == 16, row
        assert row["gap_bases"] == "T", row
        assert row["cat_orf"] == "intact", row
        assert row["repaired_orf"] != "intact", row
        assert row["ref_checked"] and row["ref_intron_gaps"] == 0, row
        assert row["grch38_pos"] == ("chr7:107" if strand == "+" else "chr7:109"), row


def test_reference_intron_is_kept():
    # GRCh38 itself has a 1-bp CDS "intron" here (e.g. a programmed frameshift): not a dropped insertion.
    meta = {"t": {"seqid": "c", "strand": "+", "source_transcript": "ENST2.3"}}
    ref = {"ENST2": ("1", "+", [(100, 105), (107, 115)])}
    (row,) = scan_rows(meta, {"t": [(5, 10), (12, 20)]}, ref=ref)
    assert row["ref_intron_gaps"] == 1 and row["gaps"] == "" and not row["frame_changing"], row


def test_gtf_includes_stop_codon():
    with tempfile.NamedTemporaryFile("w", suffix=".gtf") as fh:
        for kind, s, e in (("CDS", 101, 106), ("CDS", 110, 112), ("stop_codon", 113, 115)):
            fh.write(f'7\tensembl\t{kind}\t{s}\t{e}\t.\t+\t0\ttranscript_id "ENST3"; transcript_version "1";\n')
        fh.flush()
        assert read_gtf(fh.name) == {"ENST3": ("7", "+", [(101, 106), (110, 112), (113, 115)])}


def test_long_gap_is_intron():
    blocks = {"t": [(1, 30), (61, 90)]}  # 30-bp gap: kept as an intron
    meta = {"t": {"seqid": "c", "strand": "+"}}
    assert list(scan_rows(meta, blocks)) == []


if __name__ == "__main__":
    test_insertion_repair()
    test_reference_intron_is_kept()
    test_gtf_includes_stop_codon()
    test_long_gap_is_intron()
    print("ok")
