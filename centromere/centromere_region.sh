
# 1. split_genome_to_100kb_windows
bedtools makewindows -g genome.lengths.txt -w 100000 > windows_100kb.bed
# 2. sequence_to_bed
awk 'BEGIN{OFS="\t"} !/^#/ {print $1, $4-1, $5, $9}' trf_results.gff > trf.bed
awk 'BEGIN{OFS="\t"} !/^#/ {print $1, $4-1, $5, $9}' repeatmasker_results.gff > repeatmasker.bed

# 3. get_coverage
bedtools coverage -a windows_100kb.bed -b trf.bed > trf_cov.bed
bedtools coverage -a windows_100kb.bed -b repeatmasker.bed > repeatmasker_cov.bed

# 4. find_target_region
awk 'BEGIN{OFS="\t"} FNR==NR{trf[$1":"$2]=$7;next} 
{key=$1":"$2; if(key in trf && $7<20 && trf[key]>50) print $1,$2,$3,trf[key],$7}' \
trf_cov.bed repeatmasker_cov.bed | bedtools merge -c 4,5 -o distinct,mean > divergent_regions.bed

