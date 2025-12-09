# python3 slice_csv.py 10 20 < input.csv > output.csv

#!/usr/bin/env python3
import sys

def slice_csv_lines(lines, start_idx, end_idx):
    """
    lines: 入力ファイルの各行（リスト）
    start_idx, end_idx: 0-origin の開始・終了インデックス（inclusive）
    """
    # 1行目はそのまま出力
    header = lines[0].rstrip("\n")
    print(header)

    for line in lines[1:]:
        parts = line.strip().split()
        if not parts:
            continue
        name = parts[0]
        values = parts[1:]
        # start_idx–end_idx までを含めてスライス
        sliced = values[start_idx:end_idx + 1]
        # 名前とスライスした値をスペース区切りで出力
        print(name, *sliced)

def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} start_idx end_idx", file=sys.stderr)
        sys.exit(1)

    try:
        start_idx = int(sys.argv[1])
        end_idx = int(sys.argv[2])
    except ValueError:
        print("start_idx と end_idx は整数を指定してください", file=sys.stderr)
        sys.exit(1)

    # stdin から全行読み込み
    lines = sys.stdin.readlines()
    if len(lines) < 2:
        print("データ行が存在しません", file=sys.stderr)
        sys.exit(1)

    slice_csv_lines(lines, start_idx, end_idx)

if __name__ == "__main__":
    main()