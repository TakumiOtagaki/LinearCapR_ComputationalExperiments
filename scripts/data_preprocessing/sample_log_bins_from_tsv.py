import math
import random

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--tsv', required=True)
    parser.add_argument('--bins', type=int, default=20)
    parser.add_argument('--k', type=int, default=1)
    parser.add_argument('--out', required=True)
    parser.add_argument('--seed', type=int, default=42)
    args = parser.parse_args()

    random.seed(args.seed)
    # TSVからname, lengthを読み込む
    records = []
    with open(args.tsv) as f:
        next(f)  # skip header
        for line in f:
            name, length = line.strip().split('\t')
            records.append((name, int(length)))

    if not records:
        print("No records found.")
        return

    min_len = min(length for _, length in records)
    max_len = max(length for _, length in records)
    bin_edges = [int(math.exp(math.log(min_len) + (math.log(max_len) - math.log(min_len)) * i / args.bins)) for i in range(args.bins + 1)]

    selected = []
    for i in range(args.bins):
        bin_min = bin_edges[i]
        bin_max = bin_edges[i + 1] - 1
        bin_records = [name for name, length in records if bin_min <= length <= bin_max]
        if bin_records:
            sample_n = min(args.k, len(bin_records))
            selected.extend(random.sample(bin_records, sample_n))

    with open(args.out, 'w') as f:
        for name in selected:
            f.write(name + '\n')

if __name__ == '__main__':
    main()
