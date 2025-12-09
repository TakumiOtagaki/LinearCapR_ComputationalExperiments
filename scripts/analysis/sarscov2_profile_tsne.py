import os
import sys
import numpy as np
import pandas as pd
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
from ..common.config import WORKDIR

# 設定
# mode = "covid_only" # or "total"
mode = "total"
beamsize = 100
ID_LIST_PATH = WORKDIR / f"data/raw/sarscov2/{mode}_id_list.txt"
PROFILE_DIR = WORKDIR / "result/profile/sarscov2s"
EXCLUDE_ID = "NC_006213.1"
OUTPUT_DIR = WORKDIR / "result/covid_lincapr_msa"
os.makedirs(OUTPUT_DIR, exist_ok=True)

def load_id_list(id_list_path, exclude_id=None):
    with open(id_list_path) as f:
        ids = [line.strip() for line in f if line.strip()]
    if exclude_id is not None:
        ids = [i for i in ids if i != exclude_id]
    return ids

def load_profile_csv(profile_path):
    # 1行目がヘッダ（>で始まる）、2行目以降が "ラベル 数値 数値 ..." 形式のスペース区切り
    rows = []
    with open(profile_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            parts = line.split()
            # 先頭はラベルなので除外
            nums = []
            for p in parts[1:]:
                try:
                    nums.append(float(p))
                except Exception:
                    continue
            if nums:
                rows.append(nums)
    # 7行N列のデータを1次元ベクトルにflatten
    if len(rows) == 0:
        return np.array([])
    arr = np.array(rows)
    return arr.flatten()

def main():
    ids = load_id_list(ID_LIST_PATH, EXCLUDE_ID)
    profiles = []
    valid_ids = []
    for id_ in ids:
        print(f"Processing ID: {id_}")
        csv_path = os.path.join(PROFILE_DIR, f"{id_}.b{beamsize}.csv")
        if not os.path.exists(csv_path):
            print(f"Warning: {csv_path} not found, skipping.")
            continue
        try:
            vec = load_profile_csv(csv_path)
            if len(vec) == 0:
                print(f"Warning: {csv_path} is empty after filtering, skipping.")
                continue
            profiles.append(vec)
            valid_ids.append(id_)
        except Exception as e:
            print(f"Error loading {csv_path}: {e}")
            continue

    if len(profiles) == 0:
        print("No valid profiles loaded. Exiting.")
        return
    print([f"{id_}: {vec.shape}" for id_, vec in zip(valid_ids, profiles)])
    # 長さが異なる場合はpadding（0埋め）で揃える
    maxlen = max(len(v) for v in profiles)
    X = np.zeros((len(profiles), maxlen))
    for i, v in enumerate(profiles):
        X[i, :len(v)] = v
    # 値が小さすぎるところは 0 に置き換えておく
    X[X < 1e-3] = 0.0

    perplexity = min(30, max(1, len(X) - 1))
    # t-SNE
    print("size", X.shape, "perplexity", perplexity)
    tsne = TSNE(n_components=2, random_state=42, perplexity = perplexity)
    X_embedded = tsne.fit_transform(X)

    # 結果保存
    df_out = pd.DataFrame({
        "id": valid_ids,
        "tsne1": X_embedded[:, 0],
        "tsne2": X_embedded[:, 1]
    })
    out_csv = os.path.join(OUTPUT_DIR, f"sarscov2_profile_tsne.{mode}.b{beamsize}.csv")
    df_out.to_csv(out_csv, index=False)
    print(f"t-SNE結果を {out_csv} に保存しました。")

    # 可視化
    plt.figure(figsize=(8, 6))
    plt.scatter(X_embedded[:, 0], X_embedded[:, 1], alpha=0.7)
    for i, txt in enumerate(valid_ids):
        plt.annotate(txt, (X_embedded[i, 0], X_embedded[i, 1]), fontsize=6, alpha=0.5)
    plt.title("SARS-CoV-2 Profile t-SNE")
    plt.xlabel("t-SNE 1")
    plt.ylabel("t-SNE 2")
    plt.tight_layout()
    out_png = os.path.join(OUTPUT_DIR, f"sarscov2_profile_tsne.{mode}.b{beamsize}.png")
    plt.savefig(out_png, dpi=300)
    print(f"t-SNEプロットを {out_png} に保存しました。")

if __name__ == "__main__":
    main()
