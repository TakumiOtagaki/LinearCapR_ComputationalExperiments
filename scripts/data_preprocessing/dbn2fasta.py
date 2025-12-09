import os

input_dir = "data/processed/bprna/dbn_pagenumber1"

for filename in os.listdir(input_dir):
    if filename.endswith(".dbn"):
        dbn_path = os.path.join(input_dir, filename)
        with open(dbn_path, "r") as f:
            lines = f.readlines()
            if len(lines) < 4:
                continue
            # 1行目から sequence name を取得
            if lines[0].startswith("#1Name:"):
                seq_name = lines[0].strip().split(":", 1)[1].strip()
            else:
                seq_name = "unknown"
            # 4行目が sequence
            seq = lines[3].strip()
        fasta_path = os.path.join(input_dir, filename.replace(".dbn", ".fasta"))
        with open(fasta_path, "w") as out:
            out.write(f">{seq_name}\n{seq}\n")
