
#map_with_other_reads
for i in {001..040} ; do
    for j in {001..040}; do
        file_1="pass.ul.100k.98.part_${i}.fq.gz"
        file_2="pass.ul.100k.98.part_${j}.fq.gz"
        out="first_ont_${i}_${j}_aln.sam"
        out_detail="first_ont_${i}_${j}_aln.detail"
        out_filter_detail="first_ont_${i}_${j}_aln.filter.detail"
        out_filter_detail_filter_inside="first_ont_${i}_${j}_aln.filter_inside.detail"
        base_dir=$(dirname "$file_1")
        job_name=$(basename "${file_1%.fq.gz}")
        job_log="${base_dir}/minimap_${job_name}.ava.log"
        job_err="${base_dir}/minimap_${job_name}.ava.err"
        bsub -J "$job_name" -q fatA -n 15 -o "$job_log" -e "$job_err" \
        "minimap2 -ax map-ont -t 15 '$file_1' '$file_2' > '$out' "
        bsub -J "${job_name}" -q fatA -n 1 -o "${job_log}" -e "${job_err}" \
         "python detail_from_sam.py '$out' '$out_detail'"
        bsub -J "${job_name}" -q fatA -n 1 -o "${job_log}" -e "${job_err}" \
         "python filter_detail_rm_indentity_095.py ont_len.txt '$out_detail' > '$out_filter_detail'"
        grep "head" "$out_filter_detail"  > "$out_filter_detail_filter_inside"
        grep "tail" "$out_filter_detail"  >> "$out_filter_detail_filter_inside"
    done
done

#map_to_genome
minimap2 -ax map-ont --eqx -t 100 genome.fasta pass.ul.100k.98.fq.gz > 100k_98_ONT_correct_map.sam 
samtools view -bS 100k_98_ONT_correct_map.sam 1> 100k_98_ONT_correct_map.bam 2> 100k_98_ONT_correct_map.bam.log
samtools sort -@ 100 100k_98_ONT_correct_map.bam -o 100k_98_ONT_correct_map_sort.bam 2> 100k_98_ONT_correct_map_sort.bam.log
samtools index 100k_98_ONT_correct_map_sort.bam


#gap_region_reads_inf
python  01.find_GAP_side_500kb.py | awk '$5>=0.95&&$6>=0.98&&($10==0||$10==16) {print $0}' | sort -k7,7 -k8,8n | uniq >  01.find_GAP_side_500kb_sort_uniq.detail
#Create AGP files containing gap information
python 02.create_gap_path_agp.py 02.gap.list > 03.GAP.agp 
#Generate paths of reads surrounding gaps
python 03.ctg_side_reads_paths.py --agp 03.GAP.agp --gfa 01.find_GAP_side_500kb_sort_uniq.detail
#Split paths for each gap
mkdir split
split -l 1 chr_paths split/prefix_
#fill gap
for file in $(find . -name "prefix_ad"); do
    base_dir=$(dirname "$file")
    job_name=$(basename "$file")
    job_log="${base_dir}/fa_${job_name}.log"
    job_err="${base_dir}/fa_${job_name}.err"
    bsub -J "${job_name}" -q standardA -n 5 -o "${job_log}" -e "${job_err}"  python fill_gap  "$file"
done

