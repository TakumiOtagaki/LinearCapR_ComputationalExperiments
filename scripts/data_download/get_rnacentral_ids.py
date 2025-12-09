# 結局 7.1 GB の全データをダウンロードしたのでこのスクリプトは使わない

import argparse
import random
import requests
import math

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bins', type=int, default=20)
    parser.add_argument('--min', type=int, default=4000)
    parser.add_argument('--max', type=int, default=985945)
    parser.add_argument('--out_ids', required=True)
    parser.add_argument('--seed', type=int, default=42)
    args = parser.parse_args()

    random.seed(args.seed)
    bin_edges = [int(math.exp(math.log(args.min) + (math.log(args.max) - math.log(args.min)) * i / args.bins)) for i in range(args.bins + 1)]
    ids = []
    for i in range(args.bins):
        min_len = bin_edges[i]
        max_len = bin_edges[i + 1] - 1
        url = f"https://rnacentral.org/api/v1/rna/?min_length={min_len}&max_length={max_len}&page_size=100"
        resp = requests.get(url)
        if resp.status_code != 200:
            continue
        results = resp.json().get('results', [])
        if not results:
            continue
        rec = random.choice(results)
        acc = rec.get('rnacentral_id', None)
        if acc:
            ids.append(acc)

    with open(args.out_ids, 'w') as out:
        for acc in ids:
            out.write(acc + '\n')

if __name__ == '__main__':
    main()
