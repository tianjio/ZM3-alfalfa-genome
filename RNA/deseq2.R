
#!/usr/bin/env Rscript

# ========== 0. 加载库 ==========
suppressMessages({
  library(tximport)
  library(readr)
  library(dplyr)
  library(DESeq2)
  library(tibble)
  library(pheatmap)
  library(rhdf5)   # 读取 abundance.h5 必需
})

# ========== 1. 命令行参数解析 ==========
args <- commandArgs(trailingOnly = TRUE)

if(length(args) < 6){
  stop("Usage: Rscript pipeline_deseq2_h5.R <output_dir> <output_prefix> <allele_pair_file> <base_dir> <sample_list_file> <tx2gene_file>")
}

output_dir <- args[1]          # 输出目录
output_prefix <- args[2]       # 输出文件前缀
allele_file <- args[3]         # 等位基因对文件
base_dir <- args[4]            # abundance.h5 所在基本目录
sample_list_file <- args[5]    # 样本名列表，每行一个样本
tx2gene_file <- args[6]        # 转录本 -> 基因映射文件

# 如果输出目录不存在，创建
if(!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
}

# 输出文件名
result_file <- file.path(output_dir, paste0(output_prefix, "_DESeq2.csv"))
matrix_file <- file.path(output_dir, paste0(output_prefix, "_counts.csv"))
up_file <- file.path(output_dir, paste0(output_prefix, "_upregulated.csv"))
down_file <- file.path(output_dir, paste0(output_prefix, "_downregulated.csv"))
heatmap_file <- file.path(output_dir, paste0(output_prefix, "_heatmap.pdf"))

cat("将使用以下输出文件：\n")
cat("DESeq2结果:", result_file, "\n")
cat("表达矩阵:", matrix_file, "\n")
cat("上调基因:", up_file, "\n")
cat("下调基因:", down_file, "\n")
cat("热图文件:", heatmap_file, "\n")

# ========== 2. 读取样本名 ==========
samples <- readLines(sample_list_file)
cat("样本列表：", paste(samples, collapse=", "), "\n")

# abundance.h5 文件路径
files <- file.path(base_dir, samples, "abundance.h5")
names(files) <- samples

# ========== 3. 加载转录本 -> 基因映射 ==========
tx2gene <- read.table(
  tx2gene_file,
  header = FALSE, col.names = c("TXNAME","GENEID"), sep="\t", stringsAsFactors=FALSE
)

# ========== 4. 读取等位基因对文件 ==========
allele_pairs <- read.table(allele_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(allele_pairs) <- c("gene_A", "gene_B")

# ========== 5. 使用 tximport 读取 abundance.h5 ==========
txi <- tximport(files,
                type = "kallisto",
                tx2gene = tx2gene,
                countsFromAbundance = "lengthScaledTPM",
                importer = "h5")  # 指定读取 h5

# 转换为数据框
gene_counts <- as.data.frame(txi$counts) %>%
  rownames_to_column("gene_id")

# ========== 6. 过滤有效等位基因对 ==========
filtered_pairs <- allele_pairs %>%
  filter(gene_A %in% gene_counts$gene_id &
         gene_B %in% gene_counts$gene_id)

# ========== 7. 构建等位基因表达矩阵 ==========
allele_matrix <- matrix(0,
                        nrow = nrow(filtered_pairs),
                        ncol = length(samples)*2,
                        dimnames = list(
                          paste(filtered_pairs$gene_A, filtered_pairs$gene_B, sep="_vs_"),
                          paste0(rep(samples, each=2), "_", rep(c("A","B"), times=length(samples)))
                        ))

for(i in seq_len(nrow(filtered_pairs))){
  geneA <- filtered_pairs$gene_A[i]
  geneB <- filtered_pairs$gene_B[i]
  
  allele_matrix[i, paste0(samples,"_A")] <- as.numeric(gene_counts[gene_counts$gene_id==geneA, samples])
  allele_matrix[i, paste0(samples,"_B")] <- as.numeric(gene_counts[gene_counts$gene_id==geneB, samples])
}

# ========== 8. 样本信息 ==========
sample_info <- data.frame(
  sample = rep(samples, each=2),
  allele = rep(c("A","B"), times=length(samples)),
  row.names = paste0(rep(samples, each=2), "_", rep(c("A","B"), times=length(samples)))
)

# ========== 9. DESeq2 ==========
dds <- DESeqDataSetFromMatrix(
  countData = round(allele_matrix),
  colData = sample_info,
  design = ~ allele
)

dds <- DESeq(dds)
res <- results(dds, contrast = c("allele","A","B"))

# ========== 10. 保存结果 ==========
res_df <- as.data.frame(res) %>%
  rownames_to_column("allele_pair")
write.csv(res_df, result_file, row.names = FALSE, quote = FALSE)
write.csv(allele_matrix, matrix_file, quote = FALSE)

# ========== 11. 筛选显著差异基因 ==========
sig_genes_df <- subset(res_df, abs(log2FoldChange) > 1 & padj < 0.05)
up_genes <- subset(sig_genes_df, log2FoldChange > 0)
down_genes <- subset(sig_genes_df, log2FoldChange < 0)

write.csv(up_genes, up_file, row.names = FALSE)
write.csv(down_genes, down_file, row.names = FALSE)

# ========== 12. 热图 ==========
if(nrow(sig_genes_df) > 0){
  vsd <- vst(dds, blind=FALSE)
  mat <- assay(vsd)[rownames(sig_genes_df), ]
  
  pheatmap(mat,
           annotation_col = sample_info["allele"],
#           scale="row",
           show_rownames=FALSE,
           filename = heatmap_file)
} else {
  message("No significant genes found for heatmap.")
}

cat("✅ DESeq2 analysis using abundance.h5 完成\n")
