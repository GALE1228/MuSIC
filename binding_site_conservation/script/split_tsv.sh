#!/bin/bash

# 文件名
input_file="/home/bioinfo/02_project/05_RBP/05_homo/script/anno_multiCol.tsv"

# 获取列数
num_cols=$(awk -F'\t' '{print NF; exit}' $input_file)

# 循环遍历每一列，从第二列开始
for i in $(seq 2 $num_cols)
do
    # 使用awk提取第一列和当前列（i列）
    awk -F'\t' -v col="$i" '{print $1, $col}' OFS='\t' $input_file > "/home/bioinfo/02_project/05_RBP/05_homo/script/anno_multiCol_$i.tsv"
done
