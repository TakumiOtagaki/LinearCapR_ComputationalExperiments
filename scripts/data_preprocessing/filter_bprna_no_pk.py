import csv
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()

    with open(args.input, newline='') as infile, open(args.output, 'w', newline='') as outfile:
        reader = csv.DictReader(infile, delimiter='\t')
        writer = csv.DictWriter(outfile, fieldnames=reader.fieldnames, delimiter='\t')
        writer.writeheader()
        for row in reader:
            # PageNumber列が1ならpseudoknotなし
            if 'PageNumber' in row and row['PageNumber'] == '1':
                writer.writerow(row)

if __name__ == '__main__':
    main()
