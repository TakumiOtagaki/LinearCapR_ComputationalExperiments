import os
import shutil

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--indir', required=True, help='dbnファイルのディレクトリ')
    parser.add_argument('--outdir', required=True, help='抽出先ディレクトリ')
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    for fname in os.listdir(args.indir):
        if fname.endswith('.dbn'):
            path = os.path.join(args.indir, fname)
            with open(path, encoding='utf-8') as f:
                for line in f:
                    if line.strip() == '#PageNumber: 1':
                        shutil.copy2(path, os.path.join(args.outdir, fname))
                        break

if __name__ == '__main__':
    main()
