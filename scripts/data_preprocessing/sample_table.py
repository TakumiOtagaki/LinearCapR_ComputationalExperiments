import csv
import argparse
import random

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True)
    parser.add_argument('--n', type=int, required=True)
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--fasta', required=True)
    parser.add_argument('--table', required=True)
    args = parser.parse_args()

    random.seed(args.seed)
    with open(args.input, newline='') as infile:
        reader = list(csv.DictReader(infile, delimiter='\t'))
        # 長さフィルタ
        filtered = [row for row in reader if 19 <= int(row['Length']) <= 4021]
        sampled = random.sample(filtered, min(args.n, len(filtered)))

    # FASTA出力
    with open(args.fasta, 'w') as fasta_out:
        for row in sampled:
            fasta_out.write(f">{row['ID']}\n{row['Sequence']}\n")

    # TSV出力
    with open(args.table, 'w', newline='') as table_out:
        writer = csv.DictWriter(table_out, fieldnames=reader[0].keys(), delimiter='\t')
        writer.writeheader()
        for row in sampled:
            writer.writerow(row)

if __name__ == '__main__':
    main()
