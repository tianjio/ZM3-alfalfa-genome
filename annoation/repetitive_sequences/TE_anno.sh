#RepeatModeler
/public-supool/home/jytian/miniconda3/envs/EDTA2/share/RepeatModeler/BuildDatabase  -name genome_lib  genome.fasta
/public-supool/home/jytian/miniconda3/envs/EDTA2/share/RepeatModeler/RepeatModeler -pa 5 -database genome_lib -rmblast_dir /public-supool/home/jytian/miniconda3/envs/repeatmasker/bin -trf_prgm /public-supool/software/bin/trf -recon_dir /public-supool/software/RECON-1.08/bin -rscout_dir /public-supool/home/jytian/miniconda3/envs/EDTA2/bin 

#EDTA
perl /public-supool/home/jytian/bin/software/EDTA/EDTA.pl --genome genome.fasta -species others -cds CDS_longest.fa  -step all -t 50  --force 1

#LTR_retriever
ltr_finder -D 15000 -d 1000 -L 7000 -l 100 -p 20 -C -M 0.9 genome.fasta  > genome.finder.scn 
gt suffixerator -db  genome.fasta -indexname genome -tis -suf -lcp -des -ssp -sds -dna
gt ltrharvest -index genome -similar 90 -vic 10 -seed 20 -seqids yes -minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6 -motif TGCA -motifmis 1  > genome.harvest.scn

#filter_family
blat ltr_genome.fa.mod.LTRlib.fa modeler_ZM3-families.fa -t=dna -q=dna > modler_map_retriever.psl
python filter_result.py > lib.fasta

#RepeatMasker
RepeatMasker -lib lib.fasta -pa 24 -xsmall genome.fasta -html -gff -dir genome_TE

