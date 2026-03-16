

#blast
bsub -q standardA -n 5 -o out.log -e out.err "blastn -query genome.fasta -db 144bp.fasta -outfmt 6 -out MS_map_144bp_cent_seq.out -num_threads 10 -evalue 0.001"
bsub -q standardA -n 5 -o out.log -e out.err "blastn -query genome.fasta -db 149bp.fasta -outfmt 6 -out MS_map_149bp_cent_seq.out -num_threads 10 -evalue 0.001"
bsub -q standardA -n 5 -o out.log -e out.err "blastn -query genome.fasta -db 168bp.fasta -outfmt 6 -out MS_map_168bp_cent_seq.out -num_threads 10 -evalue 0.001"
bsub -q standardA -n 5 -o out.log -e out.err "blastn -query genome.fasta -db 187bp.fasta -outfmt 6 -out MS_map_180bp_cent_seq.out -num_threads 10 -evalue 0.001"
bsub -q standardA -n 5 -o out.log -e out.err "blastn -query genome.fasta -db 545bp.fasta -outfmt 6 -out MS_map_545bp_cent_seq.out -num_threads 10 -evalue 0.001"

#提取对应位置及其方向
for file in $(find . -name "MS_map_*_cent_seq.out"); do
    out="${file%_cent_seq.out}.bed"
    base=$(basename "$file")
    index=$(echo "$base" | awk -F'_' '{print $3}' | sed "s/bp//g")
    if [[ "$index" =~ ^[0-9]+$ ]] && [ "$index" -gt 200 ]; then
        echo "Index is greater than 200: $index"
        # 执行后续命令
        python  00.find_target_monoer_bed_and_forward.py "$file" 100 > "$out"
    else
        python 00.find_target_monoer_bed_and_forward.py "$file" 500 > "$out"
    fi
done

##数目及长度
for dir in $(find . -name "*" -type d);do
    base=$(basename "$dir")
    file="${dir}/MS_map_${base}bp.bed"
    out="${dir}/01.MS_map_${base}bp.num"
    if [[ ! -f "$file" ]]; then
        echo "File not found: $file"
        continue
    fi
       python 01.find_target_num.py "$file" | sort -k1,1 > "$out"
done

##生成簇和簇的数目和位置
for file in $(find . -name "*.bed");do
    basedir=$(dirname "$file")
    base=$(basename "$file")
    index=$(echo "$base" | awk -F'_' '{print $3}' | sed "s/bp.bed//g")
    out_file="${basedir}/02.MS_map_${index}_add_name.bed"
    awk -v name="$index" '{print $1"\t"$2"\t"$3"\t"name"\t"$4}' "$file" > "$out_file"
done

find . -name "02.MS_map_*_add_name.bed" -exec cat {} \; | sort -k1,1 -k2,2n > merged_all.bed

#统计这些簇在染色体上的位置
python 02.find_breakpoint.py merged_all.bed | sort -k1,1 -k3,3n | awk '$6>=10 {print $0}' > 01.break_point.bed

#分割相同单体连续区间的单体序列
##（1）检查是否有错误分割区间
python 03.find_continue_bed.py 01.break_point.bed
##（2）提取对应单体区间的bed文件
###生成标准bed文件
python 04.target_cluster_bed.py 01.break_point.bed > 02.break_point_standard.bed


###bedtools 提取ovlp
while IFS= read -r line; do
    echo "$line"
    file="${line}/MS_map_${line}bp.bed"
    out="${line}/${line}.ovlp.bed"
    bsub -q standardA -n 1 -o out.log -r out.err "bedtools intersect -a 02.break_point_standard.bed -b $file -wa -wb > $out"
done < "00.repeat.list"

###拆分ovlp文件并生成以cluster为最小单位的fasta文件
while IFS= read -r line; do
    echo "$line"
    dir="${line}"
    cd "$dir"
    file="${line}.ovlp.bed"
    python 05.split_chr_cluster_bed.py "$file"
done < "00.repeat.list"



###以染色体为单位的最小单位的fasta文件
#!/bin/bash
# 输入文件路径
input_file="00.repeat.list"
# 检查输入文件是否存在且非空
if [[ ! -f "$input_file" || ! -s "$input_file" ]]; then
    echo "Error: Input file is missing or empty: $input_file"
    exit 1
fi


while IFS= read -r line; do
    if [[ -z "$line" ]]; then
        echo "Warning: Empty line encountered, skipping..."
        continue
    fi
    for i in {1..8}; do
        for j in {1..4}; do
            # 构造目录名
            dir_name="01.chr/${line}/Chr${i}/Chr${i}_${j}"

            if [[ ! -d "$dir_name" ]]; then
                mkdir -p "$dir_name"
                echo "Created directory: $dir_name"
            else
                echo "Directory already exists: $dir_name"
            fi
            out_file="${dir_name}/Chr${i}_${j}.fasta"
            search_dir="${line}/Chr${i}/Chr${i}_${j}"
            fi
            out_file="${dir_name}/Chr${i}_${j}.fasta"
            search_dir="${line}/Chr${i}/Chr${i}_${j}"
            if [[ ! -d "$search_dir" ]]; then
                echo "Warning: Search directory does not exist: $search_dir"
                continue
            fi
            find "$search_dir" -name "*.fasta" -exec cat {} + >> "$out_file"
            if [[ -f "$out_file" && -s "$out_file" ]]; then
                echo "Successfully created: $out_file"
            else
                echo "Warning: Failed to create or empty output file: $out_file"
            fi
        done
    done
done < "$input_file"

#create whole genome cent_fa
while IFS= read -r line; do
    dir="02.whole_genome/${line}"
    out_file="${dir}/${line}_merged.fasta"
    search_dir="01.chr/${line}/"
    find "$search_dir" -name "*.fasta" -exec cat {} + >> "$out_file"
done < "00.repeat.list" 


#identity
base_dir="cent168"

for i in {1..8}; do
    for j in {1..4}; do
        ref="Chr${i}_${j}"
        
        for k in {1..8}; do
            for q in {1..4}; do
                qry="Chr${k}_${q}"
                
                # 创建目录
                map_dir="${base_dir}/map_self/${ref}_map_all_chr/${qry}"
                mkdir -p "$map_dir"
                
                # 文件路径
                ref_fasta="${base_dir}/split_chr/${ref}.fasta"
                qry_fasta="${base_dir}/split_chr/${qry}.fasta"
                out_file="${map_dir}/${ref}_map_${qry}_out_result"
                script_file="${map_dir}/${ref}_vs_${qry}_bsub.sh"
                
                # 生成脚本
                cat > "$script_file" << EOF
#!/bin/bash
#BSUB -J ${ref}_vs_${qry}
#BSUB -q fatA
#BSUB -n 1
#BSUB -o ${map_dir}/job.out
#BSUB -e ${map_dir}/job.err

python map_pair_identity.py \\
    "$ref_fasta" \\
    "$qry_fasta" \\
    "$map_dir" > "$out_file"
EOF
                
                chmod +x "$script_file"
            done
        done
    done
done
