# ./LinCapR <input_file> <output_file>
# inputs: /large/otgk/LinCapR_Experiments/data/raw/sarscov2/${id}.fasta
    # id は /large/otgk/LinCapR_Experiments/data/raw/sarscov2/total_id_list.txt の各行
# outputs: /large/otgk/LinCapR_Experiments/result/profile/sarscov2s/${id}.csv

ID_LIST="/large/otgk/LinCapR_Experiments/data/raw/sarscov2/total_id_list.txt"
INPUT_DIR="/large/otgk/LinCapR_Experiments/data/raw/sarscov2"
OUTPUT_DIR="/large/otgk/LinCapR_Experiments/result/profile/sarscov2s"
LINCAPR="/large/otgk/LinCapR_Experiments/LinearCapR/LinCapR"
beamsize=100
NUM_CPU=4  # 並列実行数をここで指定

mkdir -p "$OUTPUT_DIR"

process_id() {
    id="$1"
    input_file="${INPUT_DIR}/${id}.fasta"
    output_file="${OUTPUT_DIR}/${id}.b${beamsize}.csv"
    log_file="${OUTPUT_DIR}/${id}.b${beamsize}.log"
    echo "Processing $id...(Input: $input_file, Output: $output_file)"
    /usr/bin/time -v $LINCAPR "$input_file" "$output_file" "$beamsize" \
        > "${log_file}" 2>&1
    # 実行時間と最大メモリ消費量を抽出して記録
    elapsed=$(grep "Elapsed (wall clock) time" "${log_file}" | awk '{print $8}')
    mem=$(grep "Maximum resident set size" "${log_file}" | awk '{print $6}')
    echo -e "${id}\t${elapsed}\t${mem}" >> "${OUTPUT_DIR}/resource_usage.tsv"
}

export -f process_id
export INPUT_DIR OUTPUT_DIR LINCAPR beamsize

# ヘッダーを出力
echo -e "id\telapsed_time\tmax_memory_kb" > "${OUTPUT_DIR}/resource_usage.tsv"

cat "$ID_LIST" | xargs -P $NUM_CPU -I{} bash -c 'process_id "$@"' _ {}
