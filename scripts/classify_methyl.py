import csv

# 读取参考基因组序列
def read_fasta(fasta_file):
    fasta_dict = {}
    with open(fasta_file, 'r') as f:
        chrom = ""
        seq = ""
        for line in f:
            if line.startswith(">"):
                if chrom != "":
                    fasta_dict[chrom] = seq
                chrom = line[1:].strip()
                seq = ""
            else:
                seq += line.strip().upper()
        fasta_dict[chrom] = seq
    return fasta_dict

# 分类CG、CHG、CHH甲基化位点，并标注甲基化类型
def classify_methylation(txt_file, fasta_dict, output_file):
    with open(txt_file, 'r') as f, \
         open(output_file, 'w', newline='') as out_file:

        reader = csv.reader(f, delimiter='\t')
        writer = csv.writer(out_file, delimiter='\t')

        # 写入新文件的表头，添加 "C_type" 列
        writer.writerow(['chrom', 'start', 'end', 'methylation_level', 'methylated_reads', 'unmethylated_reads', 'C_type'])

        for row in reader:
            chrom, start, end, methylation_level, meth_reads, unmeth_reads = row
            start = int(start) - 1  # 0-based indexing for Python
            if chrom in fasta_dict:
                seq = fasta_dict[chrom]
                if start >= 0 and start + 2 < len(seq):  # 确保context长度足够长
                    context = seq[start:start + 3]  # 获取上下文，至少3个碱基
                else:
                    context = ""

                # 根据上下文分类并标注
                if context.startswith("CG"):
                    c_type = "CG"
                elif len(context) == 3 and context[1] != "G" and context[2] == "G":
                    c_type = "CHG"
                else:
                    c_type = "CHH"

                # 写入新的一行，包含甲基化类型
                writer.writerow([chrom, start + 1, end, methylation_level, meth_reads, unmeth_reads, c_type])
# 主函数
if __name__ == "__main__":
    fasta_file = "../ngctrl1/ngctrl1.fasta"  # 参考基因组的文件
    txt_file = "suptongctrl.txt"  # 原始的 txt 文件
    output_file = "suptongctrl_with_methylation_types.txt"  # 输出文件

    fasta_dict = read_fasta(fasta_file)  # 读取参考基因组
    classify_methylation(txt_file, fasta_dict, output_file)  # 分类并写入新的 txt 文件