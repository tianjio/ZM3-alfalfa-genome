base_dir="new_v9"
target_dir="allele_new"

for i in {1..4};do
        gff_ref_file="${base_dir}/00.primary_file/gff/ZM3_gene_chr_H${i}_updated.gff3"
        for j in {1..4};do
                if [ $i -lt $j ];then
                        gff_qry_file="${base_dir}/00.primary_file/gff/ZM3_gene_chr_H${j}_updated.gff3"
                        jcvi_dir="${base_dir}/04.jcvi/ZM3_${i}_vs_ZM3_${j}"
                        temp_dir="${target_dir}/jcvi/ZM3_${i}_vs_ZM3_${j}"
                        #mkdir -p "$temp_dir"
                        #echo "$jcvi_dir"
                        anchors="${jcvi_dir}/ZM3_${i}.ZM3_${j}.anchors"
                        out="${temp_dir}/00.filter_homo.anchors_1"
                        out_1="${temp_dir}/01.find_gene_pos_1"
                        multi="${temp_dir}/02.00.multi_genes_pos_1"
                        uniq_multi="${temp_dir}/02.00.multi_genes_pos_1_uniq"
                        single="${temp_dir}/02.01.single_gene_pos_1"
                        uniq_single="${temp_dir}/02.01.single_gene_pos_1_uniq"
                        best_pos="${temp_dir}/03.all_best_match_pos_1_1"
                        #python ${target_dir}/00.script/jcvi/00.filter_homo_gene_2.py $gff_ref_file $gff_qry_file $anchors > $out
                        #python ${target_dir}/00.script/jcvi/01.find_gene_pos_v1.py $gff_ref_file $gff_qry_file  $out > $out_1
                        #python ${target_dir}/00.script/jcvi/02.filter_multi_map_pos_anchor.py $out_1 $multi $single
                        #sort -k1,1 -k5,5 $single | uniq > $uniq_single
                        #sort -k1,1 -k5,5 $multi | uniq > $uniq_multi
                        #awk '{print $2"\t"$3"\t"$4"\t"$5"\t"$7"\t"$8"\t"$9"\t"$10}' $uniq_single > $best_pos
                fi
        done
done
