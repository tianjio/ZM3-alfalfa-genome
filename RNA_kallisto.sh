# 1. 参考注释建立索引
~/tools/kallisto_linux_v0.46.1/kallisto/kallisto index --index=ZM3_mRNA_251221_kallisto_v0.46_index M056_mRNA_cleanname_251221.fa

# 1.软件输出out
[build] loading fasta file M056_mRNA_cleanname_251221.fa
[build] k-mer length: 31
[build] warning: clipped off poly-A tail (longer than 10)
        from 1 target sequences
[build] counting k-mers ... done.
[build] building target de Bruijn graph ...  done 
[build] creating equivalence classes ...  done
[build] target de Bruijn graph has 8169672 contigs and contains 278687153 k-mers 

# 2. kallisto运行 下面是脚本循环多个样本提交slurm脚本
#!/bin/bash

# 样本名列表
samples=("leaf_1" "leaf_2" "leaf_3" \
         "stem_1" "stem_2" "stem_3" \
         "root_1" "root_2" "root_3")

# Kallisto 索引
index="/public1/home/shuairun/02_alfalfa/02_M056/00_genome/00.gene_annotation_251221/ZM3_mRNA_251221_kallisto_v0.46_index"

# 数据目录
data_dir="/public1/home/shuairun/02_alfalfa/09_M056_root_leaf_stem_TA/01_data/CleanData"

# kallisto 路径
kallisto_bin="/public1/home/shuairun/tools/kallisto_linux_v0.46.1/kallisto/kallisto"

for sample in "${samples[@]}"; do

cat > kallisto_${sample}.slurm <<EOF
#!/bin/bash
#SBATCH -J kallisto_${sample}
#SBATCH -N 1
#SBATCH -c 20
#SBATCH -p amd9654
#SBATCH --mem=50G
#SBATCH -o kallisto_${sample}_%j.out

echo "Start ${sample} at \$(date)"

${kallisto_bin} quant \
    --index ${index} \
    --output-dir ${sample} \
    --threads 20 \
    -b 100 \
    ${data_dir}/${sample}/${sample}_1.fq.gz \
    ${data_dir}/${sample}/${sample}_2.fq.gz

echo "Finished ${sample} at \$(date)"
EOF

    # 提交任务
    sbatch kallisto_${sample}.slurm
done

# 3. R整理表达量表格
# allele TPM
library(dplyr)
library(data.table)

# 读入数据并重命名 TPM 列
leaf_1 <- read.table("../05_kallisto_251222/leaf_1/abundance.tsv", 
                     header = TRUE)[, c("target_id", "tpm")]
colnames(leaf_1)[2] <- "leaf_1"

leaf_2 <- read.table("../05_kallisto_251222/leaf_2/abundance.tsv", 
                     header = TRUE)[, c("target_id", "tpm")]
colnames(leaf_2)[2] <- "leaf_2"

leaf_3 <- read.table("../05_kallisto_251222/leaf_3/abundance.tsv", 
                     header = TRUE)[, c("target_id", "tpm")]
colnames(leaf_3)[2] <- "leaf_3"

root_1 <- read.table("../05_kallisto_251222/root_1/abundance.tsv", 
                     header = TRUE)[, c("target_id", "tpm")]
colnames(root_1)[2] <- "root_1"

root_2 <- read.table("../05_kallisto_251222/root_2/abundance.tsv", 
                     header = TRUE)[, c("target_id", "tpm")]
colnames(root_2)[2] <- "root_2"

root_3 <- read.table("../05_kallisto_251222/root_3/abundance.tsv", 
                     header = TRUE)[, c("target_id", "tpm")]
colnames(root_3)[2] <- "root_3"

stem_1 <- read.table("../05_kallisto_251222/stem_1/abundance.tsv", 
                     header = TRUE)[, c("target_id", "tpm")]
colnames(stem_1)[2] <- "stem_1"

stem_2 <- read.table("../05_kallisto_251222/stem_2/abundance.tsv", 
                     header = TRUE)[, c("target_id", "tpm")]
colnames(stem_2)[2] <- "stem_2"

stem_3 <- read.table("../05_kallisto_251222/stem_3/abundance.tsv", 
                     header = TRUE)[, c("target_id", "tpm")]
colnames(stem_3)[2] <- "stem_3"

# 合并成一个表达矩阵
all_samples_trans_tpm <- list(leaf_1, leaf_2, leaf_3,
                    root_1, root_2, root_3,
                    stem_1, stem_2, stem_3)

merged_data <- Reduce(function(x, y) merge(x, y, by = "target_id", all = TRUE), 
                      all_samples_trans_tpm)

# 可选：保存结果
write.table(merged_data, file = "transcript_tpms_all_samples_251222.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

merged_data$gene_id <- sub("\\.t.*", "", merged_data$target_id)

# 计算基因层级的 TPM 加和
# 转换为 data.table
dt <- as.data.table(merged_data)

# 去掉非数值列（如 target_id），只保留 gene_id + 数值列
# 注意：这里假设 target_id 是非数值列，其它都是样本
gene_tpm <- dt[, lapply(.SD, sum, na.rm = TRUE), by = gene_id, .SDcols = is.numeric]

write.table(gene_tpm, file = "gene_tpms_all_samples_251222.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)


