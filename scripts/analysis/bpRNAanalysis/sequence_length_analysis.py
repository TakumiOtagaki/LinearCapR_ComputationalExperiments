#!/usr/bin/env python3
# multifasta の配列長を boxplot にするだけの極小スクリプト
import sys
from pathlib import Path
from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt

def main():
    fasta = Path("data/processed/bprna/bpRNA_multifasta.fasta")
    out = Path("result/bprna/sequence_len_dist.png")

    lengths = [len(r.seq) for r in SeqIO.parse(fasta.open(), "fasta")]
    if not lengths:
        sys.exit("No sequences found.")

    out.parent.mkdir(parents=True, exist_ok=True)

    # boxplot + histogram（同一PNGに2枚並べる）
    fig, (ax_box, ax_hist) = plt.subplots(1, 2, figsize=(8, 4))

    # boxplot
    ax_box.boxplot(lengths, vert=True, widths=0.5)
    ax_box.set_ylabel("Sequence length (nt)")
    ax_box.set_yscale("log")
    ax_box.set_title(f"{fasta.name} — boxplot")
    ax_box.set_xticks([])

    # histogram
    ax_hist.hist(lengths, bins="auto")
    ax_hist.set_xlabel("Sequence length (nt)")
    # x を log scale にする
    ax_hist.set_xscale("log")
    ax_hist.set_ylabel("Count")
    ax_hist.set_title("histogram")

    fig.tight_layout()
    fig.savefig(out, dpi=300)


if __name__ == "__main__":
    main()