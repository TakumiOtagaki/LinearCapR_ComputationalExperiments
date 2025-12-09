import requests

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--ids', required=True)
    parser.add_argument('--out', required=True)
    args = parser.parse_args()

    with open(args.ids) as f:
        ids = [line.strip() for line in f if line.strip()]

    with open(args.out, 'w') as out_f:
        for acc in ids:
            url = f"https://rnacentral.org/api/v1/rna/{acc}?format=fasta"
            resp = requests.get(url)
            if resp.status_code == 200:
                out_f.write(resp.text)
            else:
                print(f"Failed to fetch {acc}: {resp.status_code}")

if __name__ == '__main__':
    main()
