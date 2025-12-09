import pandas as pd
import numpy as np

def read_lincapr_profile(file_path):
    """
    LinCapRのプロファイルファイルを読み込んでDataFrameに変換する
    
    Parameters:
    -----------
    file_path : str
        読み込むファイルのパス
    
    Returns:
    --------
    pd.DataFrame
        position, bulge, exterior, hairpin, internal, multibranch, stem の列を持つDataFrame
    """
    
    # ファイルを読み込み（1行目はヘッダーとしてスキップ）
    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    # 1行目（ヘッダー）をスキップ
    data_lines = lines[1:]
    
    # 各行を解析してデータを取得
    structure_types = []
    data_matrix = []
    
    for line in data_lines:
        line = line.strip()
        if not line:  # 空行はスキップ
            continue
        
        parts = line.split()
        structure_type = parts[0]  # 最初の要素が構造タイプ
        values = [float(x) for x in parts[1:]]  # 残りの要素が数値データ
        
        structure_types.append(structure_type)
        data_matrix.append(values)
    
    # DataFrameを作成
    df = pd.DataFrame(data_matrix, index=structure_types)
    
    # 転置して、各列が塩基位置、各行が構造タイプになるようにする
    df_transposed = df.T
    
    # positionカラムを追加（0から始まるインデックス）
    df_transposed.insert(0, 'position', range(len(df_transposed)))
    
    # カラム名を小文字に変換
    df_transposed.columns = ['position'] + [col.lower() for col in df_transposed.columns[1:]]
    
    return df_transposed

# 使用例
if __name__ == "__main__":
    # ファイルパスを指定
    file_path = "/large/otgk/LinCapR_Experiments/result/profile/NC_045512.2_5UTR/localCapR_5UTR_100.tsv"
    
    # データを読み込み
    df = read_lincapr_profile(file_path)
    
    # 結果を表示
    print("データの形状:", df.shape)
    print("\nカラム名:")
    print(df.columns.tolist())
    
    print("\n最初の10行:")
    print(df.head(10))
    
    print("\n最後の5行:")
    print(df.tail(5))
    
    print("\n各構造タイプの統計:")
    structure_columns = ['bulge', 'exterior', 'hairpin', 'internal', 'multibranch', 'stem']
    for col in structure_columns:
        if col in df.columns:
            print(f"{col}: mean={df[col].mean():.6f}, std={df[col].std():.6f}")
