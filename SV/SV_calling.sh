#!/bin/bash

threads=100
ref=HapA.fasta
queries=("HapB.fasta" "HapC.fasta" "HapD.fasta")

for q in "${queries[@]}"
do
    prefix=$(basename $ref .fasta)_map_$(basename $q .fasta)

    # 1. nucmer 全基因组比对
    nucmer --mum -t $threads -p $prefix $ref $q

    # 2. 过滤最佳比对
    delta-filter -1 ${prefix}.delta > ${prefix}.filter.delta

    # 3. 统计比对信息
    dnadiff -d ${prefix}.filter.delta

    # 4. 输出坐标
    show-coords -THrd ${prefix}.filter.delta > ${prefix}.filter.coords

    # 5. Syri 检测结构变异
    syri -k \
        -s /public-supool/home/jytian/miniconda3/bin/show-snps \
        -c ${prefix}.filter.coords \
        -d ${prefix}.filter.delta \
        -r $ref \
        -q $q

done
