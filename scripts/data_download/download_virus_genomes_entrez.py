#!/usr/bin/env python
"""
Download viral GenBank records from NCBI nuccore using Biopython + Entrez.

- Reads a CSV file with an 'accession' column (e.g. NC_045512.2).
- Fetches each record in GenBank format and saves as data/raw/virus/genbank/{accession}.gbff
- Skips accessions that are already downloaded.

Usage:
    python scripts/data_download/download_virus_genomes_entrez.py \
        --input data/raw/virus/metadata/virus_list.csv \
        --outdir data/raw/virus/genbank \
        --email YOUR_EMAIL \
        [--api-key YOUR_NCBI_API_KEY]
"""

import argparse
import csv
import os
import sys
import time
from pathlib import Path

from Bio import Entrez


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        required=True,
        help="Path to CSV file containing at least an 'accession' column.",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Directory to save downloaded GenBank files.",
    )
    parser.add_argument(
        "--email",
        required=True,
        help="Email address for NCBI Entrez (required by NCBI policy).",
    )
    parser.add_argument(
        "--api-key",
        default=None,
        help="NCBI API key (optional but recommended).",
    )
    parser.add_argument(
        "--sleep",
        type=float,
        default=0.4,
        help="Sleep seconds between requests (to respect NCBI rate limits).",
    )
    return parser.parse_args()


def read_accessions(csv_path: str) -> list[str]:
    accessions: list[str] = []
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        if "accession" not in reader.fieldnames:
            raise ValueError("CSV must contain an 'accession' column.")
        for row in reader:
            acc = row["accession"].strip()
            if acc:
                accessions.append(acc)
    return accessions


def fetch_genbank(accession: str) -> str:
    """
    Fetch a GenBank record for given accession from NCBI nuccore.

    Returns the raw GenBank text.
    """
    handle = Entrez.efetch(
        db="nuccore",
        id=accession,
        rettype="gbwithparts",  # includes features; gb or gbwithparts are both ok
        retmode="text",
    )
    text = handle.read()
    handle.close()
    return text


def main():
    args = parse_args()

    Entrez.email = args.email
    if args.api_key:
        Entrez.api_key = args.api_key

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    accessions = read_accessions(args.input)
    print(f"Found {len(accessions)} accessions in {args.input}")

    for i, acc in enumerate(accessions, start=1):
        out_path = outdir / f"{acc}.gbff"
        if out_path.exists():
            print(f"[{i}/{len(accessions)}] {acc}: already exists, skipping")
            continue

        print(f"[{i}/{len(accessions)}] Fetching {acc} ...", end="", flush=True)
        try:
            gb_text = fetch_genbank(acc)
        except Exception as e:
            print(f" FAILED: {e}", file=sys.stderr)
            continue

        # very simple sanity check: GenBank records start with 'LOCUS'
        if not gb_text.lstrip().startswith("LOCUS"):
            print(f" WARNING: unexpected content for {acc}, not starting with LOCUS", file=sys.stderr)

        with open(out_path, "w") as f:
            f.write(gb_text)

        print(f" saved to {out_path}")
        time.sleep(args.sleep)


if __name__ == "__main__":
    main()