import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def load_and_process_data(file_path):
    # 1行目は使わない
    df = pd.read_csv(file_path, skiprows=1, sep=" ", index_col=0, header=None).T
    return df

def plot_motif_graph(df, graph_path, pos_name):
    plt.figure(figsize=(18, 6))
    for motif in df.columns: #  Bulge, Exterior, Hairpin, Internal, Multibranch, Stem.
        plt.plot(df.index, df[motif], label=motif)
    plt.title(f"{pos_name} Motifs")
    plt.xlabel("Position")
    plt.ylabel("Score")
    plt.grid(True, linestyle='--', alpha=0.5, which='both')
    # x, y の 軸の値を細かく増やす
    # plt.xticks(np.arange(1, len(df.index) + 1, step=10))
    # plt.yticks(np.arange(0, df.max().max() + 0.1, step=0.1))
    # plt.xlim(1, len(df.index))
    plt.legend()
    plt.savefig(graph_path)

def main():
    # ファイルパスを指定
    file_path = Path("result/profile/NC_045512.profile.w500.txt")
    graph_dir = Path("result/NC_045512_analysis")
    
    # データを読み込み、処理
    df = load_and_process_data(file_path)
    
    # 結果を表示
    print("Data loaded and processed successfully:")
    print(df.head())  # 最初の5行を表示

    pos_name_map = {
        "5UTR": (1, 266), # 1 origin (end は含まない)
        "3UTR": (29675, 29903) # 1 origin (end は含まない) 
    }

    for pos_name, (start_pos, end_pos) in pos_name_map.items():
        print(f"\n{pos_name} Position: {start_pos}..{end_pos}")
        df_pos = df.iloc[start_pos - 1:end_pos - 1]
        print("\n5' UTR Data:")
        # print(df_pos)

        # plot
        plot_motif_graph(df_pos, graph_dir / f"{pos_name}_motifs.png", pos_name)
        print("graph saved to:", graph_dir / f"{pos_name}_motifs.png")

        # df filtered: T 未満となっている値は 0 に書き換える
        T = 0.05
        df_filtered = df_pos.copy()
        df_filtered[df_filtered < T] = 0
        print("\nFiltered Data (T < 0.05):")
        # print(df_filtered)
        plot_motif_graph(df_filtered, graph_dir / f"{pos_name}_motifs_filtered.png", pos_name)
        print("Filtered graph saved to:", graph_dir / f"{pos_name}_motifs_filtered.png")
       


    # stem_loop       29609..29644

    # 3'UTR           29675..29903

if __name__ == "__main__":
    main()