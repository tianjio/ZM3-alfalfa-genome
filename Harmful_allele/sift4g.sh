for i in {1..4};do
  perl make-SIFT-db-all.pl -config ZM3_${i}.txt
done

/public-supool/home/jytian/bin/software/sift4g/bin/sift4g -t 50 -d /public-supool/home/jytian/bin/software/uniprot_sprot.fasta -q /public-supool/home/jytian/bin/software/scripts_to_build_SIFT_db/test_files/homo_sapiens_small/all_prot.fasta --subst /public-supool/home/jytian/bin/software/scripts_to_build_SIFT_db/test_files/homo_sapiens_small/test_subst --out /public-supool/home/jytian/bin/software/scripts_to_build_SIFT_db/test_files/homo_sapiens_small/SIFT_predictions --sub-results
