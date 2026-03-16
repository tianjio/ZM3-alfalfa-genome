#hifiasm ont+hifi
hifiasm -o hifiasm_hifi_ont --hg-size 3.2g  -t 70  --write-ec --write-paf  -l 0 --n-hap 4  --ul-rate 0.02 -ul /public-supool/home/jytian/M056_ont/merged_four_150k_98correct.fq.gz   /public-supool/home/jytian/M056hifi/M056.merged.Q30.reads.fq  2> hifiasm.log
#hifiasm ont
hifiasm -t 220 --ont -l 0 --hg-size 3.2g  --n-hap 4  -o hifiasm_ont_only  /public-supool/home/jytian/M056_ont/merged_four_150k_98correct.fq.gz  2>hifiasm.log
#verkko
verkko -d verkko --hifi /public-supool/home/jytian/M056hifi/M056.merged.Q30.reads.fq --nano /public-supool/home/jytian/M056_ont/merged_four_150k_98correct.fq --threads 70 2>02.verkko1.sh1.o
