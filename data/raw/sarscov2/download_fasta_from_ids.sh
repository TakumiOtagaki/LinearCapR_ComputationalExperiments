#!/usr/bin/env bash
# id_list.txt を１行ずつ読み込んで fasta を取得し、<アクセッション>.fasta に保存

# while read -r id; do
#   echo "Downloading $id …"
#   efetch -db nuccore -id $id -format fasta > "${id}.fasta"
# done < id_list.txt


# multi fasta 
# out=/large/otgk/LinCapR_Experiments/result/covid_msa/total.multi.fa
# : > "$out"
# while read -r id; do
#   cat "writing ${id}.fasta to ${out}"
#   cat ${id}.fasta >> $out
# done < /large/otgk/LinCapR_Experiments/data/raw/sarscov2/total_id_list.txt

# multi fasta 
out=/large/otgk/LinCapR_Experiments/result/covid_msa/covid_only.multi.fa
echo -n > "$out"
while read -r id; do
  echo "writing ${id}.fasta to ${out}:"
  cat /large/otgk/LinCapR_Experiments/data/raw/sarscov2/${id}.fasta >> $out
done < /large/otgk/LinCapR_Experiments/data/raw/sarscov2/covid_only_id_list.txt