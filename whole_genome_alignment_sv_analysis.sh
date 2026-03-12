#mummer 全基因组比对
nucmer --mum -t 100 -p HapA_map_HapB HapA.fasta HapB.fasta

delta-filter  -1 HapA_map_HapB.delta > HapA_map_HapB.filter.delta

dnadiff -d HapA_map_HapB.filter.delta

show-coords -THrd HapA_map_HapB.filter.delta > HapA_map_HapB.filter.coords

#syri 提取结构变异
syri -k -s /public-supool/home/jytian/miniconda3/bin/show-snps -c HapA_map_HapB.filter.coords  -d HapA_map_HapB.filter.delta -r HapA.fasta  -q HapB.fasta