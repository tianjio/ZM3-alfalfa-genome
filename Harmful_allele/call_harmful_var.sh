snpeff_dir="/public-supool/home/jytian/bin/software/snpEff"
dir="Mt_A17_v6.0"
A17_v6_genome="${dir}/Mt_A17_v6_sed_name.fasta"
A17_v6_cds="${dir}/Mt_A17_v6.cds.fasta"
A17_v6_pro="${dir}/Mt_A17_v6.pep.fasta"
A17_v6_gff="${dir}/A17_T2T.v6.0.fa.gff"

cp "$A17_v6_genome" "${snpeff_dir}/data/genomes/A17_v6.fa"
A17_v6_dir="${snpeff_dir}/data/A17_v6"

mkdir -p "$A17_v6_dir"
cp "$A17_v6_genome" "${A17_v6_dir}/sequences.fa"
cp "$A17_v6_cds" "${A17_v6_dir}/cds.fa"
cp "$A17_v6_pro" "${A17_v6_dir}/protein.fa"
cp "$A17_v6_gff" "${A17_v6_dir}/genes.gff"

java -jar /public-supool/home/jytian/bin/software/snpEff/snpEff.jar build -v A17_v6 2>&1 | tee A17_v6.build

echo "HM078.genome:HM078" >>   /public-supool/home/jytian/bin/software/snpEff/snpEff.config
snpeff_dir="/public-supool/home/jytian/bin/software/snpEff" 
dir="/public-supool/home/jytian/Project/jytian/Medicago_project/08_other_genome/Mt_HM078/"
HM078_genome="${dir}/HM078.chr.fasta"
HM078_cds="${dir}/HM078.cds.fasta"
HM078_pro="${dir}/HM078.pep.longest.fasta"
HM078_gff="${dir}/HM078_Chr.gene.gff3"

cp "$HM078_genome" "${snpeff_dir}/data/genomes/HM078.fa"
HM078_dir="${snpeff_dir}/data/HM078"

mkdir -p "$HM078_dir"
cp "$HM078_genome" "${HM078_dir}/sequences.fa"
cp "$HM078_cds" "${HM078_dir}/cds.fa"
cp "$HM078_pro" "${HM078_dir}/protein.fa"
cp "$HM078_gff" "${HM078_dir}/genes.gff"

java -jar /public-supool/home/jytian/bin/software/snpEff/snpEff.jar build -v HM078 2>&1 | tee HM078.build


for i in {1..4};do
    echo "ZM3_${i}_v8.genome:ZM3_${i}_v8" >> /public-supool/home/jytian/bin/software/snpEff/snpEff.config

    snpeff_dir="/public-supool/home/jytian/bin/software/snpEff" 
    dir="/public-supool/home/jytian/Project/jytian/Medicago_project/02_annotation/01_gene/new_v8"
    ZM3_genome="/public-supool/home/jytian/Project/jytian/Medicago_project/01_genome/04.01_final_genome_two_version/v4_final/Hap${i}_genome.fasta"
    ZM3_cds="${dir}/01.seq_file/ZM3_gene_chr_H${i}_update.cds.fasta"
    ZM3_pro="${dir}/01.seq_file/ZM3_gene_chr_H${i}_update.pep.fasta"
    ZM3_gff="${dir}/00.primary_file/gff/ZM3_gene_chr_H${i}_updated.gff3"

    cp "$ZM3_genome" "${snpeff_dir}/data/genomes/ZM3_${i}_v8.fa"
    ZM3_dir="${snpeff_dir}/data/ZM3_${i}_v8"

    mkdir -p "$ZM3_dir"
    cp "$ZM3_genome" "${ZM3_dir}/sequences.fa"
    cp "$ZM3_cds" "${ZM3_dir}/cds.fa"
    cp "$ZM3_pro" "${ZM3_dir}/protein.fa"
    cp "$ZM3_gff" "${ZM3_dir}/genes.gff"

    java -jar /public-supool/home/jytian/bin/software/snpEff/snpEff.jar build -v ZM3_${i}_v8 2>&1 | tee ZM3_${i}_v8.build
done

%%bash
#!/bin/bash

base_dir="/public-supool/home/jytian/Project/jytian/Medicago_project/07_comparative_genome/00_mummer_new_v4_final"
chr_dir="/public-supool/home/jytian/Project/jytian/Medicago_project/01_genome/04.01_final_genome_two_version/v4_final/split_fa"
show_snps="/public-supool/home/jytian/miniconda3/envs/mummer4_new/bin/show-snps"

for k in {1..8}; do
  for i in {1..4}; do
    for j in {1..4}; do
      if [ "$i" -ne "$j" ]; then

        mummer_dir="${base_dir}/chr_${k}/${i}_vs_${j}"
        mkdir -p "$mummer_dir"

        pushd "$mummer_dir" > /dev/null 
        pwd
        prefix="Chr${k}_${i}_map_Chr${k}_${j}"
        bsub_file="${mummer_dir}/mummer.sh"
        echo "$bsub_file"
        rm "$bsub_file"
        cd $mummer_dir
        cat << EOF > "$bsub_file"
#BSUB -J mummer_chr${k}_${i}_${j}
#BSUB -q standardA
#BSUB -n 2
#BSUB -o ${mummer_dir}/mummer.log
#BSUB -e ${mummer_dir}/mummer.err

#nucmer --mum -t 20 -p $prefix \
#  ${chr_dir}/Chr${k}_${i}.fasta \
#  ${chr_dir}/Chr${k}_${j}.fasta

#delta-filter -1 ${prefix}.delta > ${prefix}.filter.delta
#dnadiff -d ${prefix}.filter.delta
#show-coords -THrd ${prefix}.filter.delta > ${prefix}.filter.coords

syri -k \
  -s $show_snps \
  -c ${prefix}.filter.coords \
  -d ${prefix}.filter.delta \
  -r ${chr_dir}/Chr${k}_${i}.fasta \
  -q ${chr_dir}/Chr${k}_${j}.fasta
EOF

        bsub < "$bsub_file"
        popd > /dev/null 
        sleep 3
      fi
    done
  done
done
%%bash
#!/bin/bash
base_dir="/public-supool/home/jytian/Project/jytian/Medicago_project/08_other_genome/00.mummer/HM078_vs_A17_v6"
HM078_chr="/public-supool/home/jytian/Project/jytian/Medicago_project/08_other_genome/Mt_HM078/split_fa"
A17_chr="/public-supool/home/jytian/Project/jytian/Medicago_project/08_other_genome/Mt_A17_v6.0/split_fa"
show_snps="/public-supool/home/jytian/miniconda3/envs/mummer4_new/bin/show-snps"

for i in {1..8};do
    mummer_dir="${base_dir}/chr_${i}/"
    mkdir -p "$mummer_dir"

    pushd "$mummer_dir" > /dev/null
    pwd
    prefix="HM078_Chr${i}_map_A17_Chr${i}"
    bsub_file="${mummer_dir}/mummer.sh"
    rm "$bsub_file"
    cat << EOF > "$bsub_file"
#BSUB -J mummer_ch_${i}
#BSUB -q standardA
#BSUB -n 2
#BSUB -o ${mummer_dir}/mummer.log
#BSUB -e ${mummer_dir}/mummer.err

#nucmer --mum -t 20 -p $prefix \
#  ${HM078_chr}/Chr${i}.fasta \
#  ${A17_chr}/Chr${i}.fasta

#delta-filter -1 ${prefix}.delta > ${prefix}.filter.delta
#dnadiff -d ${prefix}.filter.delta
#show-coords -THrd ${prefix}.filter.delta > ${prefix}.filter.coords

#syri -k \
#  -s $show_snps \
#  -c ${prefix}.filter.coords \
#  -d ${prefix}.filter.delta \
#  -r ${HM078_chr}/Chr${i}.fasta \
#  -q ${A17_chr}/Chr${i}.fasta
EOF

    bsub < "$bsub_file"
    popd > /dev/null 
    sleep 2
done

%%bash
#!/bin/bash
base_dir="/public-supool/home/jytian/Project/jytian/Medicago_project/08_other_genome/00.mummer/HM078_vs_A17_v6"
HM078_chr="/public-supool/home/jytian/Project/jytian/Medicago_project/08_other_genome/Mt_HM078/split_fa"
A17_chr="/public-supool/home/jytian/Project/jytian/Medicago_project/08_other_genome/Mt_A17_v6.0/split_fa"
show_snps="/public-supool/home/jytian/miniconda3/envs/mummer4_new/bin/show-snps"

i=8
j=4
    mummer_dir="${base_dir}/chr_${i}_vs_${j}/"
    mkdir -p "$mummer_dir"

    pushd "$mummer_dir" > /dev/null
    pwd
    prefix="HM078_Chr${i}_map_A17_Chr${j}"
    bsub_file="${mummer_dir}/mummer.sh"
    rm "$bsub_file"
    cat << EOF > "$bsub_file"
#BSUB -J mummer_ch_${i}_map_${j}
#BSUB -q standardA
#BSUB -n 2
#BSUB -o ${mummer_dir}/mummer.log
#BSUB -e ${mummer_dir}/mummer.err

#nucmer --mum -t 20 -p $prefix \
#  ${HM078_chr}/Chr${i}.fasta \
#  ${A17_chr}/Chr${j}.fasta

#delta-filter -1 ${prefix}.delta > ${prefix}.filter.delta
#dnadiff -d ${prefix}.filter.delta
#show-coords -THrd ${prefix}.filter.delta > ${prefix}.filter.coords

#syri -k \
#  -s $show_snps \
#  -c ${prefix}.filter.coords \
#  -d ${prefix}.filter.delta \
#  -r ${HM078_chr}/Chr${i}.fasta \
#  -q ${A17_chr}/Chr${j}.fasta
EOF

    bsub < "$bsub_file"
    popd > /dev/null 

%%bash
base_dir="/public-supool/home/jytian/Project/jytian/Medicago_project/07_comparative_genome/00_mummer_new_v4_final"
snpeff_dir="${base_dir}/snpeff"
mkdir -p "$snpeff_dir"

for i in {1..4};do
    for j in {1..4};do
        if [ "$i" -ne "$j" ]; then
            outfile="${base_dir}/syri/Hap${i}_vs_Hap${j}_syri.vcf"
            if [ -f "$outfile" ]; then
                rm "$outfile"
            fi
            for k in {1..8};do
                tar_dir="${base_dir}/chr_${k}/${i}_vs_${j}"
                file="${tar_dir}/syri.vcf"
                cat "$file" >> "$outfile"
            done

            outfile_2="${base_dir}/syri/Hap${i}_vs_Hap${j}_syri.out"
            filter_2="${base_dir}/syri/Hap${i}_vs_Hap${j}_syri_filter.out"
            if [ -f "$outfile_2" ]; then
                rm "$outfile_2"
            fi 

            for k in {1..8};do
                tar_dir="${base_dir}/chr_${k}/${i}_vs_${j}"
                file="${tar_dir}/syri.out"
                cat "$file" >> "$outfile_2"
            done

            awk '$11!="SNP"&&$11!="INS"&&$11!="DEL" {print $0}' "$outfile_2" > "$filter_2"


            #bsub -q fatA -m b005 -n 1 -o out.log -e out.err "java -jar /public-supool/home/jytian/bin/software/snpEff/snpEff.jar  eff \
            #    -c /public-supool/home/jytian/bin/software/snpEff/snpEff.config \
            #    ZM3_${i}_v8 $outfile > ${snpeff_dir}/Hap${i}_vs_Hap${j}_syri.snpeff.vcf \
            #    -csvStats ${snpeff_dir}/positive_Hap${i}_vs_Hap${j}_syri.snpeff.csv \
            #    -stats ${snpeff_dir}/positive_Hap${i}_vs_Hap${j}_syri.snpeff.html"
        fi
    done
done

        %%bash
#!/bin/bash
#A17_ref
base_dir="/public-supool/home/jytian/Project/jytian/Medicago_project/08_other_genome/00.mummer/A17_v6_vs_HM078"
outfile="${base_dir}/merged_syri.vcf"
snpeff_dir="${base_dir}/snpeff"
mkdir -p "$snpeff_dir"

if [ -f "$outfile" ]; then
    rm "$outfile"
fi
for i in {1..8};do
    file="${base_dir}/chr_${i}/syri.vcf"
    if [ -f "$file" ];then
        cat "$file" >> "$outfile"
        echo "$file"
    fi
done
cat ${base_dir}/chr_4_vs_8/syri.vcf >> "$outfile"
cat ${base_dir}/chr_8_vs_4/syri.vcf >> "$outfile"

bsub -q fatA -m b005 -n 1 -o ${base_dir}/out.log -e ${base_dir}/out.err "java -jar /public-supool/home/jytian/bin/software/snpEff/snpEff.jar  eff \
                -c /public-supool/home/jytian/bin/software/snpEff/snpEff.config \
                A17_v6 $outfile > ${snpeff_dir}/syri.snpeff.vcf \
                -csvStats ${snpeff_dir}/positive_syri.snpeff.csv \
               -stats ${snpeff_dir}/positive_syri.snpeff.html"

    
    
    
    

    

