# (1) Align Hi-C data to the assembly, remove PCR duplicates and filter out secondary and supplementary alignments
#bwa index  ONT_l0.asm.bp.p_ctg.fa
bwa mem -t 80 -5SP  ONT_l0.asm.bp.p_ctg.fa /supool-old/Multi-omics_DATA/Medicago/Medicago_sativa/HiC_fromLiangLab/M056/HIC-1-1A_R1.fq.gz  /supool-old/Multi-omics_DATA/Medicago/Medicago_sativa/HiC_fromLiangLab/M056/HIC-1-1A_R2.fq.gz | samblaster | samtools view - -@ 60 -S -h -b -F 3340 -o M056-1A.bam
bwa mem -t 80 -5SP ONT_l0.asm.bp.p_ctg.fa /supool-old/Multi-omics_DATA/Medicago/Medicago_sativa/HiC_fromLiangLab/M056/HIC-1-2A_R1.fq.gz /supool-old/Multi-omics_DATA/Medicago/Medicago_sativa/HiC_fromLiangLab/M056/HIC-1-2A_R2.fq.gz  | samblaster | samtools view - -@ 80 -S -h -b -F 3340 -o M056-2A.bam

# (2) Filter the alignments with MAPQ 1 (mapping quality ≥ 1) and NM 3 (edit distance < 3)
/public-supool/home/jytian/bin/software/HapHiC/utils/filter_bam.py M056-1A.bam 1 --NM 3 --threads 80 | samtools view - -b -@ 60 -o M056-1A.filtered.bam
/public-supool/home/jytian/bin/software/HapHiC/utils/filter_bam.py M056-2A.bam 1 --NM 3 --threads 80 | samtools view - -b -@ 80 -o M056-2A.filtered.bam
# (3) Merged all bam
samtools merge merged.filterd.bam M056-1A.filtered.bam M056-2A.filtered.bam
# (4) Run HapHiC pipeline
/public-supool/home/jytian/bin/software/HapHiC/haphic pipeline  ONT_l0.asm.bp.p_ctg.fa  merged.filterd.bam 32 --RE "GATC" --correct_nrounds 2 --remove_allelic_links 4 --threads 20 --processes 20
# (5) Run jucerbox
ln -s /public-as5300/gap/jytian/Medicago_project/01_genome/01_HIC/00_hifiasm_l0/ONT_l0.asm.bp.p_ctg.fa .
samtools faidx ONT_l0.asm.bp.p_ctg.fa
/public-supool/home/jytian/bin/software/HapHiC/scripts/../utils/juicer pre -a -q 1 -o out_JBAT /public-as5300/gap/jytian/Medicago_project/01_genome/01_HIC/00_hifiasm_l0/merged.filterd.bam scaffolds.raw.agp ONT_l0.asm.bp.p_ctg.fa.fai >out_JBAT.log 2>&1
(java -jar -Xmx32G /public-supool/home/jytian/bin/software/HapHiC/scripts/../utils/juicer_tools.1.9.9_jcuda.0.8.jar pre out_JBAT.txt out_JBAT.hic.part <(cat out_JBAT.log | grep PRE_C_SIZE | awk '{print $2" "$3}')) && (mv out_JBAT.hic.part out_JBAT.hic)
# (6) Plot
/public-supool/home/jytian/bin/software/HapHiC/haphic plot scaffolds.agp  ../merged.filterd.bam --separate_plots
