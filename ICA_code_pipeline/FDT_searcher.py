import os
import subprocess
from dataclasses import dataclass
from typing import List, Dict, Tuple
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline, NcbimakeblastdbCommandline
from io import StringIO
import tempfile
import datetime
import csv
import math
import argparse
from data.scripts.ICA_code_pipeline.GFF_entry_class import GFFEntry
from miniprot_module import run_miniprot_analysis
from blast_module import run_blastn_search
from architecture_indentification import analyze_structures
from reporter_module import generate_report, generate_csv_report, generate_gff3_report



def run_analysis(
    genome_fasta_path: str,
    target_sequence_file: str,
    miniprot_path: str,
    Tpase: str,
    output_report_file: str,
    output_csv_file: str,
    output_gff3_file: str,
    evalue_threshold: float = 1e-5,
    identity_threshold: float = 60.0,
    alignment_length_threshold: int = 30,
    max_target_seqs: int = 1000,
    max_hsps: int = 10,
    threads: int = 16
):
    try:
        # 第一步：使用BLASTN识别转座子
        print("步骤1: 使用BLASTN识别转座子...")
        gff_entries = run_blastn_search(
            target_sequence_file,
            genome_fasta_path,
            evalue_threshold,
            identity_threshold,
            alignment_length_threshold,
            max_target_seqs,
            max_hsps,
            threads
        )
        print(f"找到了 {len(gff_entries)} 个转座子匹配。")

        # 第二步：提取转座子区域的序列
        print("步骤2: 提取转座子区域的序列...")
        transposable_elements = extract_transposable_elements(
            gff_entries,
            genome_fasta_path
        )
        print(f"提取了 {len(transposable_elements)} 个转座子序列。")

        # 第三步：使用miniprot识别转座酶蛋白区域
        print("步骤3: 使用miniprot识别转座酶蛋白区域...")
        miniprot_results = run_miniprot_analysis(
            miniprot_path,
            Tpase,
            transposable_elements,
            threads
        )
        print(f"miniprot找到了 {len(miniprot_results)} 个蛋白质编码区域。")

        print("步骤4: 分析转座子结构...")
        analysis_results = []
        for te_record in transposable_elements:
            # 直接从 annotations 获取位置信息
            blast_match = te_record.annotations['blast_match']

            # 使用更灵活的匹配方式
            protein_regions = []
            for miniprot_entry in miniprot_results:
                # 修改匹配逻辑
                if (miniprot_entry.seqid == te_record.id or
                    te_record.id in miniprot_entry.seqid or
                    miniprot_entry.seqid.split('_')[-1] in te_record.id):
                    protein_regions.append((miniprot_entry.start, miniprot_entry.end))

            # 只保留找到转座酶蛋白的转座子
            if protein_regions:
                try:
                    structures = analyze_structures(
                        str(te_record.seq),
                        protein_regions
                    )

                    # 检查结构是否为空
                    if not structures:
                        print(f"Warning: No structures found for TE {te_record.id}")
                        continue

                    # 计算转座子结构的精确基因组位置
                    te_structure_start = min(
                        structures[0]['TSD1'][0],
                        structures[0]['TIR1'][0],
                        structures[0]['TIR2'][0],
                        structures[0]['TSD2'][0]
                    )
                    te_structure_end = max(
                        structures[0]['TSD1'][1],
                        structures[0]['TIR1'][1],
                        structures[0]['TIR2'][1],
                        structures[0]['TSD2'][1]
                    )
                    achitechture_score = structures[0]['scores']['total_score']
                    # 计算在参考基因组中的精确位置
                    blast_position = f"{blast_match['chromosome']}:{blast_match['start']}-{blast_match['end']}"
                    element_position = f"{blast_match['chromosome']}:{blast_match['start'] + te_structure_start}-{blast_match['start'] + te_structure_end}"

                    analysis_results.append({
                        'te_id': te_record.id,
                        'blast_position': blast_position,
                        'element_position': element_position,
                        'blast_match': blast_match,
                        'te_structure_position': {
                            'start': te_structure_start,
                            'end': te_structure_end
                        },
                        'protein_regions': protein_regions,
                        'structures': structures,
                        'source': 'FDT_pipeline',
                        'architechture_score': achitechture_score,
                        'type': 'functional_transposable_element',
                        'strand': blast_match['strand'],  # 使用 BLAST 的 strand 信息
                        'score': 0.0,
                        'miniprot_strand': miniprot_entry.strand , # 可选地将 miniprot 的 strand 记录在 attributes 中
                        'phase': '.',
                    })

                except Exception as e:
                    print(f"Error processing TE {te_record.id}: {str(e)}")
                    continue

        print(f"成功分析的转座子数量: {len(analysis_results)}")

        # 步骤 5: 生成各种报告
        print("步骤5: 生成报告...")
        generate_report(analysis_results, output_report_file)
        print(f"分析报告已保存到 {output_report_file}。")

        generate_csv_report(analysis_results, output_csv_file)
        print(f"详细的CSV报告已保存到 {output_csv_file}。")

        generate_gff3_report(analysis_results, output_gff3_file)
        print(f"GFF3格式报告已保存到 {output_gff3_file}。")

    except Exception as e:
        print(f"在执行过程中发生错误: {str(e)}")
    finally:
        # 清理临时文件
        for file in [f for f in os.listdir() if f.endswith('.tmp')]:
            os.remove(file)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Transposable Element Analysis Tool")

    # 添加缩写形式的参数，并设置为可选的关键字参数
    parser.add_argument("-g", "--genome_fasta_path", help="Path to the genome FASTA file", metavar="g", nargs='?')
    parser.add_argument("-ref", "--target_sequence_file", help="Path to the target sequence file", metavar="ref", nargs='?')
    parser.add_argument("-miniprot", "--miniprot_path", help="Path to the miniprot executable", metavar="miniprot", nargs='?')
    parser.add_argument("-tp", "--Tpase", help="Path to the hAT sequences file", metavar="tp", nargs='?')
    parser.add_argument("-or", "--output_report_file", help="Path to the output report file", metavar="or", nargs='?')
    parser.add_argument("-ot", "--output_csv_file", help="Path to the output CSV file", metavar="ot", nargs='?')
    parser.add_argument("-gff", "--output_gff3_file", help="Path to the output GFF3 file", metavar="gff", nargs='?')
    parser.add_argument("--evalue_threshold", type=float, default=1e-5, help="E-value threshold for BLASTN", metavar="et", nargs='?')
    parser.add_argument("--identity_threshold", type=float, default=60.0, help="Identity threshold for BLASTN", metavar="it", nargs='?')
    parser.add_argument("--alignment_length_threshold", type=int, default=30, help="Alignment length_threshold for BLASTN", metavar="alt", nargs='?')
    parser.add_argument("--max_target_seqs", type=int, default=1000, help="Maximum number of target sequences for BLASTN", metavar="mts", nargs='?')
    parser.add_argument("--max_hsps", type=int, default=10, help="Maximum number of HSPs for BLASTN", metavar="mhs", nargs='?')
    parser.add_argument("--threads", type=int, default=16, help="Number of threads to use", metavar="t", nargs='?')

    args = parser.parse_args()

    # 检查必需参数是否都已提供
    required_args = ['genome_fasta_path', 'target_sequence_file','miniprot_path', 'Tpase', 'output_report_file', 'output_csv_file', 'output_gff3_file']
    for arg in required_args:
        if getattr(args, arg) is None:
            parser.error(f"{arg} is a required argument.")

    run_analysis(
        args.genome_fasta_path,
        args.target_sequence_file,
        args.miniprot_path,
        args.Tpase,
        args.output_report_file,
        args.output_csv_file,
        args.output_gff3_file,
        args.evalue_threshold,
        args.identity_threshold,
        args.alignment_length_threshold,
        args.max_target_seqs,
        args.max_hsps,
        args.threads
    )