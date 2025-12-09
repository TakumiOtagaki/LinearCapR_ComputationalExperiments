#!/usr/bin/env python
"""
Build a virus accession list CSV by querying NCBI nuccore (Entrez),
fetching GenBank records, and heuristically deciding whether each
genome is "annotated" or "unannotated" based on feature types.

Output: data/raw/virus/metadata/virus_list.csv

Usage:
    uv run python scripts/data_download/build_virus_list_entrez.py \
        --out data/raw/virus/metadata/virus_list.csv \
        --email $NCBI_EMAIL \
        --api-key $NCBI_API_KEY \
        --max-records 50
"""

import argparse
import csv
import time
from pathlib import Path
from typing import Iterable

from Bio import Entrez, SeqIO


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--out",
        required=True,
        help="Path to output CSV (virus_list.csv).",
    )
    parser.add_argument(
        "--email",
        required=True,
        help="Email for NCBI Entrez.",
    )
    parser.add_argument(
        "--api-key",
        default=None,
        help="NCBI API key (optional but recommended).",
    )
    parser.add_argument(
        "--max-records",
        type=int,
        default=50,
        help="Maximum number of records to fetch from Entrez search.",
    )
    parser.add_argument(
        "--sleep",
        type=float,
        default=0.5,
        help="Sleep seconds between Entrez requests.",
    )
    parser.add_argument(
        "--min_len",
        type=int,
        default=7000,
        help="Minimum genome length to keep (nt).",
    )
    parser.add_argument(
        "--max_len",
        type=int,
        default=50000,
        help="Maximum genome length to keep (nt).",
    )
    return parser.parse_args()


def entrez_search_ids(max_records: int) -> list[str]:
    """
    Search nuccore for candidate viral RNA complete genomes and return a list of IDs.
    You can tweak the search term as needed.
    """
    # Virus (txid2559587) AND "complete genome" in title.
    # Riboviria = RNAウイルスを含む realm
    # You can refine this query later to specific families, hosts, etc.
    term = 'txid2559587[Organism:exp] AND "complete genome"[Title] NOT "Severe acute respiratory syndrome coronavirus 2"[Organism]'
    handle = Entrez.esearch(
        db="nuccore",
        term=term,
        retmax=max_records,
    )
    result = Entrez.read(handle)
    handle.close()
    ids = result.get("IdList", [])
    return ids


def classify_record(record) -> dict:
    """
    Given a Biopython SeqRecord for a viral genome, return a dict with:
    - accession
    - virus_name
    - length
    - is_refseq
    - n_features
    - has_cds
    - has_rbp_like_feature
    - set_type ("annotated" or "unannotated")
    """
    accession = record.id
    desc = record.description
    length = len(record.seq)
    is_refseq = accession.startswith("NC_") or accession.startswith("NG_")

    feature_types = [f.type for f in record.features]
    n_features = len(feature_types)

    has_cds = "CDS" in feature_types

    rbp_like_types = {
        "protein_bind",
        "misc_binding",
    }

    struct_like_types = {
        "stem_loop",
        "misc_structure",
        "regulatory",
        "misc_feature",
        "ncRNA",
        "misc_RNA",
    }

    has_rbp_like = any(ft in rbp_like_types for ft in feature_types)
    has_struct_like = any(ft in struct_like_types for ft in feature_types)

    # 「CDS を持っていて、構造っぽい feature もある」= 怪しい候補
    has_candidate = has_cds and has_struct_like

    if has_cds and has_rbp_like:
        # 明示的に RBP binding (protein_bind / misc_binding) がある
        set_type = "annotated"
    elif has_candidate:
        # RBP は明示されていないが、構造エレメントとして怪しい
        set_type = "candidate"
    elif has_cds:
        # 普通の CDS だけ（構造 annotation ほぼ無し）
        set_type = "unannotated"
    else:
        # ゲノムとして扱いづらい（断片など）
        set_type = "unusable"

    # Virus name: take the first part of description before the first comma or bracket.
    virus_name = desc
    for sep in [",", "["]:
        if sep in virus_name:
            virus_name = virus_name.split(sep)[0].strip()

    print(accession, has_cds, has_rbp_like, has_struct_like, n_features, sorted(set(feature_types)))

    return {
        "accession": accession,
        "virus_name": virus_name,
        "length": length,
        "is_refseq": "yes" if is_refseq else "no",
        "n_features": n_features,
        "has_cds": "yes" if has_cds else "no",
        "has_rbp_like": "yes" if has_rbp_like else "no",
        "set_type": set_type,
    }


def fetch_genbank_records(id_list: Iterable[str], sleep: float) -> Iterable:
    """
    Fetch GenBank records for given nuccore IDs and yield SeqRecord objects.
    """
    for idx, nid in enumerate(id_list, start=1):
        print(f"[{idx}/{len(id_list)}] Fetching nuccore ID {nid} ...", flush=True)
        handle = Entrez.efetch(
            db="nuccore",
            id=nid,
            rettype="gbwithparts",
            retmode="text",
        )
        # Some records may contain multiple entries; we take them one by one.
        for record in SeqIO.parse(handle, "gb"):
            yield record
        handle.close()
        time.sleep(sleep)


def main():
    args = parse_args()

    Entrez.email = args.email
    if args.api_key:
        Entrez.api_key = args.api_key

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    ids = entrez_search_ids(args.max_records)
    print(f"Entrez search returned {len(ids)} IDs")

    rows = []
    for record in fetch_genbank_records(ids, sleep=args.sleep):
        length = len(record.seq)
        if not (args.min_len <= length <= args.max_len):
            # skip too short / too long genomes (likely fragments or huge DNA viruses)
            continue

        info = classify_record(record)
        rows.append(info)

    if not rows:
        print("No records passed filters; nothing to write.")
        return

    fieldnames = [
        "accession",
        "virus_name",
        "length",
        "is_refseq",
        "n_features",
        "has_cds",
        "has_rbp_like",
        "set_type",
    ]

    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    print(f"Wrote {len(rows)} entries to {out_path}")


if __name__ == "__main__":
    main()