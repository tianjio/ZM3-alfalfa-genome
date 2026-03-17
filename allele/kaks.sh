python split_pair.py /three_allele.tsv three_allele_split_pair.tsv
wgd ksd --n_threads 100 three_allele_split_pair.tsv ZM3_gene_update.cds.longest.rename.fasta --tmp_dir three_allele_kaks_result
