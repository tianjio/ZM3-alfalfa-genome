#!/usr/bin/env python3
import re,os
from collections import defaultdict

# =========================
# 1. 输入文件
# =========================

gff_file  = "merged.gff"
longest_cds_file="ZM3_gene_all_longest_mrna_map_gene.txt"
base_dir="05.frameshift/"


# =========================
# 1.1解析最长的转录本
# =========================

longest_cds = set()
with open(longest_cds_file) as file:
    for line in file:
        cols = line.strip().split()
        longest_cds.add(cols[0])

# =========================
# 2. 读取 CDS exon 结构
# =========================
# gene -> list of (chr, start, end, strand)
cds_exons = defaultdict(list)

with open(gff_file) as f:
    for line in f:
        if line.startswith("#"):
            continue
        c = line.strip().split("\t")
        if c[2] != "CDS":
            continue

        chr_, start, end, strand, attr = c[0], int(c[3]), int(c[4]), c[6], c[8]

        m = re.search(r'Parent=([^;]+)', attr)
        if not m:
            continue

        transcript = m.group(1)
        if transcript  in longest_cds:    
            gene = transcript.split(".")[0]

            cds_exons[gene].append((chr_, start, end, strand))

# 按转录方向排序 CDS exon
for gene in cds_exons:
    strand = cds_exons[gene][0][3]
    cds_exons[gene].sort(
        key=lambda x: x[1],
        reverse=(strand == "-")
    )

# =========================
# 3. 解析 CDS alignment
# =========================
def read_fasta(fp):
    seqs = {}
    name = None
    seq = []
    with open(fp) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(seq)
                name = line[1:]
                seq = []
            else:
                seq.append(line)
        if name:
            seqs[name] = "".join(seq)
    return seqs

def read_name(fp):
    name_list = []
    with open(fp) as f:
        for line in f:
            line.strip()
            if line.startswith(">"):
                name_list.append(line[1:])
    return name_list

def extract_number_after_g(text):
    """提取G后面的数字并转换为整数"""
    match = re.search(r'G(\d{2})', text)
    if match:
        Hap = ((int(match.group(1)) -1) % 4) + 1
        return Hap
    return None
# =========================
# 4. 定位 frameshift 在 CDS 中的位置
# =========================
def get_frameshift_cds_positions(aln_seq):
    """
    返回所有 ! 对应的 CDS 坐标（1-based）
    """
    cds_pos = 0
    positions = []

    for i, base in enumerate(aln_seq, 1):
        if base != "-":
            cds_pos += 1
        if base == "!":
            positions.append(cds_pos)

    return positions

# =========================
# 5. CDS → 基因组坐标
# =========================
def cds_to_genome(gene, cds_pos):
    """
    cds_pos: 1-based
    """
    exons = cds_exons[gene]
    offset = cds_pos

    for chr_, start, end, strand in exons:
        length = end - start + 1

        if offset <= length:
            if strand == "+":
                return chr_, start + offset - 1
            else:
                return chr_, end - offset + 1
        offset -= length

    return None

# =========================
# 6. 主程序
# =========================

outfile=os.path.join(base_dir,"merged_all_frameshift_pos.txt")

with open(outfile,'w') as out_file:
    out_file.write(f"Ref_Hap\tQuery_Hap\tGene\tCDS_pos\tChr\tGenome_pos\n")
    #print("Ref_Hap\tQuery_Hap\tGene\tCDS_pos\tChr\tGenome_pos")

    for i in range(1,5):
        have_frameshift_allele = set()
        summary_file=os.path.join(base_dir,f"ZM3_{i}","comprehensive_analysis_stats_add.csv")
        with open(summary_file,'r') as f:
            for line in f:
                cols =  line.strip().split(',')
                if "True" in cols:
                    allele = cols[0]
                    have_frameshift_allele.add(allele)
    
        for allele in have_frameshift_allele:
            cds_fasta=os.path.join(base_dir,f"ZM3_{i}",f"{allele}_macse_nt.fasta")
            seqs = read_fasta(cds_fasta)
            name_list = read_name(cds_fasta)
            another_hap = None
            for name in name_list:
                hap = extract_number_after_g(name)
                if hap != i:
                    another_hap = hap
            for gene, aln_seq in seqs.items():
                base_gene = gene.split(".")[0]

                if base_gene not in cds_exons:
                    continue

                fs_positions = get_frameshift_cds_positions(aln_seq)

                for cds_pos in fs_positions:
                    res = cds_to_genome(base_gene, cds_pos)
                    if res:
                        chr_, gpos = res
                        out_file.write(f"{i}\t{another_hap}\t{gene}\t{cds_pos}\t{chr_}\t{gpos}\n")
