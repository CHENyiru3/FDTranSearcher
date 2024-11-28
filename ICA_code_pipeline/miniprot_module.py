import os
import subprocess
from typing import List
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import tempfile
from data.scripts.ICA_code_pipeline.GFF_entry_class import GFFEntry
from dataclasses import dataclass

def run_miniprot_analysis(
    miniprot_path: str,
    Tpase: str,
    transposable_elements: List[SeqRecord],
    threads: int = 16
) -> List[GFFEntry]:
    try:
        miniprot_results = []

        # 创建一个临时文件来存储所有转座子序列
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fasta') as all_te_fasta:
            SeqIO.write(transposable_elements, all_te_fasta, "fasta")
            all_te_fasta_path = all_te_fasta.name

        # 创建一个临时文件来存储 GFF 输出
        miniprot_output = tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.gff').name

        # 运行 miniprot，使用 --gff 参数
        miniprot_cmd = f"{miniprot_path} -ut {threads} --gff {all_te_fasta_path} {Tpase} > {miniprot_output}"

        # 使用subprocess.run执行命令
        result = subprocess.run(miniprot_cmd, shell=True, capture_output=True, text=True)

        if result.returncode != 0:
            print(f"miniprot执行错误。返回码: {result.returncode}")
            print(f"标准错误输出: {result.stderr}")
            return []

        # 解析 miniprot 输出的 GFF 文件
        with open(miniprot_output, 'r') as miniprot_file:
            for line in miniprot_file:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue
                seqid, source, feature, start, end, score, strand, phase, attributes = fields
                if feature == 'CDS':
                    # 查找与该序列相关联的转座子记录
                    te_record = next((te for te in transposable_elements if te.id == seqid), None)

                    if te_record is None:
                        print(f"警告：未找到序列ID {seqid} 对应的转座子记录")
                        continue

                    # 获取 BLAST 的链方向
                    blast_strand = te_record.annotations['blast_match']['strand']

                    # 如果 BLAST 匹配的是负链，反转 miniprot 输出的链方向
                    if blast_strand == '-':
                        if strand == '+':
                            strand = '-'
                        elif strand == '-':
                            strand = '+'

                    gff_entry = GFFEntry(
                        seqid=seqid,  # 使用 GFF 文件中的 seqid
                        source='miniprot',
                        type='protein_coding_region',
                        start=int(start),
                        end=int(end),
                        score=float(score) if score != '.' else 0.0,
                        strand=strand,  # 使用调整后的链方向
                        phase=phase,
                        attributes=parse_gff_attributes(attributes)
                    )
                    miniprot_results.append(gff_entry)

        # 清理临时文件
        os.remove(all_te_fasta_path)
        os.remove(miniprot_output)

        return miniprot_results

    except Exception as e:
        print(f"miniprot过程中发生错误: {str(e)}")
        return []
