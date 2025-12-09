# 結局 7.1 GB の全データをダウンロードしたのでこのスクリプトは使わない

import argparse
import random
import requests
import math
import time

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bins', type=int, default=20)
    parser.add_argument('--min', type=int, default=4000)
    parser.add_argument('--max', type=int, default=985945)
    parser.add_argument('--out_fasta', required=True)
    parser.add_argument('--seed', type=int, default=42)
    args = parser.parse_args()

    random.seed(args.seed)
    # 対数スケールでbin分割
    bin_edges = [int(math.exp(math.log(args.min) + (math.log(args.max) - math.log(args.min)) * i / args.bins)) for i in range(args.bins + 1)]
    records = []
    for i in range(args.bins):
        print(f"Processing bin {i + 1}/{args.bins}: length {bin_edges[i]} to {bin_edges[i + 1] - 1}")
        min_len = bin_edges[i]
        max_len = bin_edges[i + 1] - 1
        url = f"https://rnacentral.org/api/v1/rna/?min_length={min_len}&max_length={max_len}&page_size=100"
        resp = requests.get(url)
        time.sleep(0.5)
        if resp.status_code != 200:
            continue
        results = resp.json().get('results', [])
        if not results:
            continue
        rec = random.choice(results)
        records.append(rec)

    with open(args.out_fasta, 'w') as fasta_out:
        for rec in records:
            seq = rec.get('sequence', '')
            acc = rec.get('rnacentral_id', 'unknown')
            fasta_out.write(f">{acc}\n{seq}\n")

if __name__ == '__main__':
    main()
