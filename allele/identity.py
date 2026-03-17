
#!/usr/bin/env python3

import os
import sys
import argparse
import subprocess
from multiprocessing import Pool, cpu_count
from Bio import SeqIO
from Bio import AlignIO

# --------------------- 输入读取函数 ---------------------

def read_allele_list(allele_file):
    """读取等位基因列表文件"""
    allele_groups = []
    with open(allele_file, 'r') as f:
        for line in f:
            genes = line.strip().split()
            if len(genes) == 4:
                filtered_genes = []
                for gene in genes:
                    if gene != "NA":
                        filtered_genes.append(gene)
                if len(filtered_genes) >= 2:
                    allele_groups.append(filtered_genes)
    return allele_groups

def read_fasta_sequences(fasta_file):
    """读取FASTA序列为字典"""
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        gene_id = record.id.split('.')[0]
        #gene_id = record.id.replace('_upstream(+)', '').replace('_upstream(-)', '')
        sequences[gene_id] = str(record.seq)
    return sequences

# -------------------- 核心功能函数 ----------------------

def run_muscle_alignment(input_fasta, output_fasta):
    """运行 MUSCLE 进行序列比对"""
    cmd = ['muscle', '-align', input_fasta, '-output', output_fasta]
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        return True
    except subprocess.CalledProcessError as e:
        print(f"[错误] MUSCLE 比对失败: {input_fasta}")
        print(f"返回码: {e.returncode}\n错误信息: {e.stderr}")
        return False
    except FileNotFoundError:
        print("错误：找不到 muscle 命令，请检查是否已安装并加入 PATH。")
        return False

def calculate_pairwise_identity(alignment_file):
    """从MUSCLE比对结果计算所有序列对的Identity"""
    try:
        alignment = AlignIO.read(alignment_file, "fasta")
        num_seqs = len(alignment)
        pairwise_results = []

        if num_seqs < 2:
            return []

        seq_names = [record.id for record in alignment]

        for i in range(num_seqs):
            seq1 = alignment[i]
            for j in range(i+1, num_seqs):
                seq2 = alignment[j]
                matches, total = 0, 0
                for k in range(alignment.get_alignment_length()):
                    a, b = seq1[k].upper(), seq2[k].upper()
                    if a != '-' and b != '-':
                        total += 1
                        if a == b:
                            matches += 1
                identity = matches / total if total else 0
                gene1 = seq_names[i].split('.')[0]
                gene2 = seq_names[j].split('.')[0]
                pairwise_results.append((gene1, gene2, identity))
        return pairwise_results
    except Exception as e:
        print(f"计算 Identity 出错: {alignment_file}，错误：{e}")
        return []

def safe_remove_file(file_path):
    try:
        if os.path.exists(file_path):
            os.remove(file_path)
    except Exception as e:
        print(f"[警告] 删除文件失败 {file_path}: {e}")

# -------------------- 并行处理函数 -----------------------

def process_allele_group(args):
    """并行处理一组等位基因序列比对与 Identity 计算"""
    group_index, gene_group, sequences_dict, output_dir = args
    group_name = f"group_{group_index}"
    fasta_path = os.path.join(output_dir, f"{group_name}.fasta")
    align_path = fasta_path.replace('.fasta', '.align.fasta')

    try:
        with open(fasta_path, 'w') as f:
            for gene_id in gene_group:
                if gene_id in sequences_dict:
                    f.write(f">{gene_id}\n{sequences_dict[gene_id]}\n")
                else:
                    print(f"[跳过] {gene_id} 缺失于FASTA序列中")
                    return []

        if not run_muscle_alignment(fasta_path, align_path):
            return []

        result = calculate_pairwise_identity(align_path)
        safe_remove_file(fasta_path)
        safe_remove_file(align_path)
        return result

    except Exception as e:
        print(f"[错误] 处理组 {group_name} 失败: {e}")
        return []

# ------------------------ 主函数 -------------------------

def main():
    parser = argparse.ArgumentParser(description='等位基因下游序列 Identity 分析')
    parser.add_argument('--allele_list', required=True, help='等位基因列表文件')
    parser.add_argument('--gene_fasta', required=True, help='下游序列FASTA文件')
    parser.add_argument('--output_dir', default='./output', help='输出目录')
    parser.add_argument('--output_file', default='all_pairs_identity.tsv', help='最终输出文件名')

    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    final_output = os.path.join(args.output_dir, args.output_file)

    print("=" * 60)
    print("等位基因下游序列 Identity 分析（并行版）")
    print("=" * 60)

    # 步骤1：读取数据
    print("读取输入文件...")
    allele_groups = read_allele_list(args.allele_list)
    sequences_dict = read_fasta_sequences(args.gene_fasta)
    print(f"共 {len(allele_groups)} 组等位基因，{len(sequences_dict)} 条序列")

    # 步骤2：并行处理每组
    print(f"\n启动并行处理（线程数: {cpu_count()}）...")
    args_list = [(i + 1, group, sequences_dict, args.output_dir)
                 for i, group in enumerate(allele_groups)]

    all_results = []
    with Pool(processes=cpu_count()) as pool:
        all_results = pool.map(process_allele_group, args_list)

    # 步骤3：汇总结果
    print(f"\n汇总结果到文件：{final_output}")
    total_pairs = 0
    with open(final_output, 'w') as out_f:
        out_f.write("Gene_A\tGene_B\tIdentity\n")
        for group_results in all_results:
            for gene1, gene2, identity in group_results:
                out_f.write(f"{gene1}\t{gene2}\t{identity:.4f}\n")
                total_pairs += 1

    print("\n分析完成！")
    print(f"共计算 {total_pairs} 对基因的 Identity")
    print(f"结果文件保存在：{final_output}")
    print("=" * 60)

if __name__ == "__main__":
    main()
