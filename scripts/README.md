# LinCapR Experiments Scripts

## 概要
LinCapR（Linear-time RNA CapR structure prediction）の実験用スクリプト集を機能別に整理したものです。

## フォルダ構成

### 📥 data_download/
データダウンロード関連のスクリプト（一度実行すれば十分）
- `download_rnacentral_fasta_by_ids.py`: RNAcentral APIから配列をダウンロード
- `download_rnacentral_bins.py`: RNAcentralのビンデータをダウンロード（非使用）
- `get_rnacentral_ids.py`: RNAcentral IDリストを生成

### 🔧 data_preprocessing/
データ前処理関連のスクリプト（rawデータを編集して保存）
- `filter_bprna_no_pk.py`: bpRNAから擬似結びを除去
- `extract_pagenumber1_dbn.py`: PageNumber=1のDBN構造を抽出
- `sample_table.py`: テーブルからランダムサンプリング
- `sample_log_bins_from_tsv.py`: 対数スケールでビンサンプリング
- `fasta_name_length_table.py`: FASTA配列の名前と長さのテーブル作成
- その他の抽出・変換スクリプト

### 📊 analysis/
解析関連のスクリプト（よく呼び出される）
- `benchmark_lincapr.py`: LinCapRのベンチマーク実行
- `time_memory/aggregate_results.py`: 時間・メモリ測定CSVの集計
- `time_memory/plot_time_memory.py`: 時間・メモリ測定の散布図作成
- `sarscov2_profile_tsne.py`: SARS-CoV-2プロファイルのt-SNE解析
- `entropy_on_sarscov2.py`: SARS-CoV-2のエントロピー解析
- `NC_045512.2_analysis.py`: 特定配列の解析
- `lincapr_sarscov2.sh`: SARS-CoV-2用LinCapR実行スクリプト
- **`rRNAanalysis/`**: 16S rRNA構造解析モジュール（進行中）

### 🖼️ generate_publication_figs/
論文図表生成パッケージ（`python -m scripts.generate_publication_figs --config graph_config.json`）
- `multiloop_unpaired.py`: bpRNA の multiloop 内連続非ペア長を集計し、C=30 のカットオフを可視化（グラフ + 集計表）

### 🛠️ common/
共通ユーティリティとツール
- `utils.py`: プロファイル読み込み、アライメント処理等の共通関数
- `config.py`: 作業ディレクトリ等の設定

## rRNAanalysis について
16S rRNAデータをLinCapRで解析するためのモジュールです。ROCカーブ作成のための正解ラベル付けを行います。

- `structure_parser.py`: RNA二次構造解析クラス（現在編集中）
- 将来: `StructureContextLabeler`クラスを実装予定

## 使用方法
各スクリプトは以下のように実行できます：

```bash
# データダウンロード例
python scripts/data_download/download_rnacentral_fasta_by_ids.py --ids id_list.txt --out sequences.fa

# 前処理例
python scripts/data_preprocessing/filter_bprna_no_pk.py --input bpRNA.tsv --output filtered.tsv

# 解析例
python scripts/analysis/benchmark_lincapr.py [オプション]
```

## 移行履歴
2025年7月29日: 元の`scripts/`フォルダから機能別に整理・移動
