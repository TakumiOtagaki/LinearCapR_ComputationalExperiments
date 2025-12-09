import subprocess
import argparse
import os
import sys
from tqdm import tqdm

from multiprocessing import Pool

def process_one(args_tuple):
    name, seq, cmd_template = args_tuple
    from tempfile import NamedTemporaryFile
    import os
    seq_len = len(seq)
    with NamedTemporaryFile('w', delete=False, suffix='.fa', dir='result/time_memory/tmp') as tmp_fa:
        tmp_fa.write(f'>{name}\n')
        for i in range(0, len(seq), 60):
            tmp_fa.write(seq[i:i+60] + '\n')
        tmp_fa.flush()
        tmp_fa_name = tmp_fa.name
    cmd = cmd_template.format(fasta=tmp_fa_name)
    # print(f'Running: {cmd}')
    result = subprocess.run(f'/usr/bin/time -v {cmd}', shell=True, stderr=subprocess.PIPE, stdout=subprocess.DEVNULL, text=True, cwd=os.getcwd())
    stderr_path = f"result/time_memory/tmp/{name}.stderr.txt"
    with open(stderr_path, "w") as debug_stderr:
        debug_stderr.write(result.stderr)
    time_sec = None
    max_rss = None
    for line in result.stderr.splitlines():
        if 'Elapsed (wall clock) time' in line:
            t = line.split(':', 1)[1].strip()
            parts = t.split()[-1].split(":")
            # print(f'Parsed time: {t} -> {parts}')
            try:
                if len(parts) == 3:
                    sec = int(parts[0])*3600 + float(parts[1]) *60 + float(parts[2])
                elif len(parts) == 2:
                    sec = int(parts[0])*60 + float(parts[1])
                else:
                    sec = float(parts[0])
                time_sec = sec
            except Exception:
                time_sec = None
        if 'Maximum resident set size' in line:
            try:
                max_rss = int(line.split(':')[1].strip())
            except Exception:
                max_rss = None
    os.remove(tmp_fa_name)
    # print(f'Finished: {name} - time: {time_sec}, max_rss: {max_rss}, length: {seq_len}')
    return (name, time_sec, max_rss, seq_len)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta', required=True, help='multi-fastaファイル')
    parser.add_argument('--cmd', required=True, help='実行コマンド（{fasta}が一時fastaファイル名に置換される）')
    parser.add_argument('--out', required=True, help='出力TSV')
    parser.add_argument('--jobs', type=int, default=1, help='並列実行数')
    parser.add_argument('--progress-interval', type=int, default=500, help='進捗を表示する単位（0で無効）')
    parser.add_argument(
        '--metadata',
        nargs='*',
        default=[],
        help='key=value 形式で追加列を指定（例: dataset=bpRNA mode=LinCapR）',
    )
    args = parser.parse_args()

    from tempfile import NamedTemporaryFile

    def parse_fasta(fp):
        name = None
        seq = []
        for line in fp:
            if line.startswith('>'):
                if name:
                    yield name, ''.join(seq)
                name = line[1:].strip().split()[0]
                seq = []
            else:
                seq.append(line.strip())
        if name:
            yield name, ''.join(seq)

    metadata_items = []
    for entry in args.metadata:
        if '=' not in entry:
            parser.error(f"Invalid metadata entry (expect key=value): {entry}")
        key, value = entry.split('=', 1)
        metadata_items.append((key, value))

    os.makedirs('result/time_memory/tmp', exist_ok=True)
    with open(args.fasta) as fasta_fp, open(args.out, 'w') as out:
        header_fields = ['fasta', 'time_sec', 'max_rss_kb', 'length']
        header_fields.extend(key for key, _ in metadata_items)
        out.write(','.join(header_fields) + '\n')
        fasta_records = [(name, seq, args.cmd) for name, seq in parse_fasta(fasta_fp)]
        results = []
        total_records = len(fasta_records)
        progress_interval = args.progress_interval
        show_progress_messages = progress_interval is not None and progress_interval > 0
        if show_progress_messages:
            progress_interval = max(1, progress_interval)
        with Pool(processes=args.jobs) as pool:
            iterator = tqdm(
                pool.imap(process_one, fasta_records, chunksize=1),
                total=total_records,
                desc="Processing",
                disable=show_progress_messages,
            )
            for idx, result in enumerate(iterator, start=1):
                results.append(result)
                if show_progress_messages and total_records:
                    if idx % progress_interval == 0 or idx == total_records:
                        pct = idx / total_records * 100
                        print(f"[progress] {idx}/{total_records} ({pct:.1f}%) completed", flush=True)
        for name, time_sec, max_rss, seq_len in results:
            row = [name, time_sec, max_rss, seq_len]
            row.extend(value for _, value in metadata_items)
            out.write(','.join(str(item) for item in row) + '\n')

if __name__ == '__main__':
    main()
