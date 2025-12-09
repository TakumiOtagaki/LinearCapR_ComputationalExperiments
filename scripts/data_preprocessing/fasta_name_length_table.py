import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--in_fasta', required=True)
    parser.add_argument('--out_tsv', required=True)
    args = parser.parse_args()

    with open(args.in_fasta, 'r') as fin, open(args.out_tsv, 'w') as fout:
        fout.write("name\tlength\n")
        name = None
        seq = []
        count = 0
        for line in fin:
            if count % 1000000 == 0:
                print(f"Processed {count} lines ({count * 100.0 / 490549792 :.2f} %)")
            if line.startswith(">"):
                if name is not None:
                    fout.write(f"{name}\t{len(''.join(seq))}\n")
                name = line[1:].strip().split()[0]
                seq = []
            else:
                seq.append(line.strip())
            count += 1
        if name is not None:
            fout.write(f"{name}\t{len(''.join(seq))}\n")

if __name__ == '__main__':
    main()
