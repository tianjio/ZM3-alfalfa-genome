pybarrnap genome.fasta -k euk -t 10 > genome_rDNA.gff
python 00.split_5s_and_45s.py genome_rDNA.gff 5S_rDNA.gff 45S_rDNA.gff 
grep -v "partial" 45S_rDNA.gff > 45S_rDNA_filter.gff 
python  01.stats_45s_step1_v1.py   45S_rDNA_filter.gff >  45S_rDNA_filter_sub.list
python 01.stats_45s_step2_v1.py 45S_rDNA_filter_sub.list > 45S_rDNA_filter_units.bed

fasta="genome.fasta"
bedtools getfasta -fi $fasta -bed 45S_rDNA_filter_units.bed -fo 45S_rDNA_filter_units.fasta -s

bsub -q fatA -n 1 -o out.log -e out.err "trf 45S_rDNA_filter_units.fasta 2 7 7 80 10 50 500 -d -h"
trf2gff -i 45S_rDNA_filter_units.fasta.2.7.7.80.10.50.500.dat

python setk_sub_repeat_from_gff_v2.py 45S_rDNA_filter_units.fasta.2.7.7.80.10.50.500.gff3 > 45S_rDNA_filter_units_repeat.bed
awk '$2>=310&&$2<=325 {print $0}' 45S_rDNA_filter_units_repeat.bed > 45S_rDNA_filter_units_repeat_filter.bed
awk '{if ($3 == 3.4) $4 = "Type1";else if ($3 == 4.4) $4 = "Type2";else if ($3 == 5.4) $4 = "Type3";print $0}' 45S_rDNA_filter_units_repeat_filter.bed > 45S_rDNA_filter_units_repeat_filter_retype.bed
python /public-supool/home/jytian/Project/jytian/Medicago_project/03_rDNA/00_rDNA_reads/HiFi/filter_retype.py 45S_rDNA_filter_units_repeat_filter_retype.bed > 45S_rDNA_filter_units_repeat_filter_retype_twice.bed
awk '{if ($4 == "") $4 = "Other";count[$4]++;total++;} END {for (type in count) {printf "%s\t%d\t%.2f%%\n", type, count[type], (count[type]/total)*100;}}' 45S_rDNA_filter_units_repeat_filter_retype_twice.bed > 45S_rDNA_filter_units_repeat_filter_retype_twice.stats
