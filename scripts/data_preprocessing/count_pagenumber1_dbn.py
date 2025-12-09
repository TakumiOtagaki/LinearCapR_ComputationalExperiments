import os

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--indir', required=True, help='dbnファイルのディレクトリ')
    args = parser.parse_args()

    count = 0
    for fname in os.listdir(args.indir):
        if fname.endswith('.dbn'):
            path = os.path.join(args.indir, fname)
            with open(path, encoding='utf-8') as f:
                for line in f:
                    if line.strip() == '#PageNumber: 1':
                        count += 1
                        break
    print(count)

if __name__ == '__main__':
    main()
