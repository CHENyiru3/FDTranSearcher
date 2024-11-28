from dataclasses import dataclass
from typing import List, Dict
from GFF_entry import GFFEntry
import datetime
import csv


def generate_report(analysis_results: List[Dict], output_file: str):
    with open(output_file, 'w') as f:
        f.write("Transposable Element Analysis Report\n")
        f.write(f"Analysis Date: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Total Transposable Elements with Protein: {len(analysis_results)}\n\n")

        for result in analysis_results:
            f.write(f"TE ID: {result['te_id']}\n")

            # 位置信息
            f.write("Genomic Positions:\n")
            f.write(f"  BLAST Match Position: {result['blast_position']}\n")
            f.write(f"  Element Structure Position: {result['element_position']}\n")

            # 详细的蛋白质区域报告
            f.write("Protein Regions:\n")
            for i, (start, end) in enumerate(result['protein_regions'], 1):
                f.write(f"  Region {i}: {start}-{end}\n")

            # 结构详情
            f.write("\nTransposable Element Structures:\n")
            for i, structure in enumerate(result['structures'][:3], 1):
                f.write(f"  Structure {i}:\n")
                for key in ['TSD1', 'TIR1', 'TIR2', 'TSD2']:
                    start, end = structure[key]
                    f.write(f"    {key}: {start}-{end} ({end - start + 1} bp)\n")
                    f.write(f"    {key} Sequence: {structure[f'{key.lower()}_seq']}\n")

                f.write("\n")

            f.write("-" * 50 + "\n\n")


def generate_csv_report(analysis_results: List[Dict], output_file: str):
    with open(output_file, 'w', newline='') as csvfile:
        fieldnames = [
            'te_id',
            'blast_chromosome',
            'blast_start',
            'blast_end',
            'element_chromosome',
            'element_start',
            'element_end',
            'protein_region_start',
            'protein_region_end'
        ]

        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for result in analysis_results:
            for protein_region in result['protein_regions']:
                writer.writerow({
                    'te_id': result['te_id'],
                    'blast_chromosome': result['blast_match'].get('chromosome', 'N/A'),
                    'blast_start': result['blast_match'].get('start', 'N/A'),
                    'blast_end': result['blast_match'].get('end', 'N/A'),
                    'element_chromosome': result['blast_match'].get('chromosome', 'N/A'),
                    'element_start': result['blast_match'].get('start', 'N/A') + result['te_structure_position']['start'],
                    'element_end': result['blast_match'].get('start', 'N/A') + result['te_structure_position']['end'],
                    'protein_region_start': protein_region[0],
                    'protein_region_end': protein_region[1]
                })


def generate_gff3_report(analysis_results: List[Dict], output_file: str):
    """
    生成规范化的 GFF3 格式报告。
    对于每个转座子，仅保留得分最高的结构预测。
    """
    with open(output_file, 'w') as f:
        # 写入 GFF3 头部
        f.write("##gff-version 3\n")

        for result in analysis_results:
            # 获取转座子基本信息
            te_id = result['te_id']
            blast_match = result['blast_match']
            te_structure_position = result['te_structure_position']
            protein_regions = result['protein_regions']
            structures = result['structures']
            architecture_score = result['architechture_score']

            # 如果没有结构预测，跳过该转座子
            if not structures:
                print(f"Warning: No structures found for TE {te_id}")
                continue

            # 筛选最佳结构（假设每个结构有一个 'score' 字段）
            best_structure = max(structures, key=lambda x: x['scores']['total_score'])

            # 提取最佳结构中的 TSD 和 TIR 信息
            attributes = []
            attributes.append(f"ID={te_id}")
            attributes.append(f"BLAST_Match={blast_match['chromosome']}:{blast_match['start']}-{blast_match['end']}")

            # 元素结构位置信息
            if te_structure_position:
                attributes.append(f"Element_Structure={blast_match['chromosome']}:{blast_match['start'] + te_structure_position['start']}-{blast_match['start'] + te_structure_position['end']}")

            # 添加蛋白区域信息
            for i, (start, end) in enumerate(protein_regions, 1):
                attributes.append(f"Protein_Region_{i}={start}-{end}")

            # 添加结构信息
            for key in ['TSD1', 'TIR1', 'TIR2', 'TSD2']:
                if key in best_structure:
                    start, end = best_structure[key]
                    seq = best_structure.get(f"{key.lower()}_seq", "")
                    attributes.append(f"{key}={start}-{end};{key}_Sequence={seq}")

            # 写入 GFF3 内容
            f.write(f"{blast_match['chromosome']}\t{result['source']}\t{result['type']}\t"
                    f"{blast_match['start']}\t{blast_match['end']}\t{result['score']}\t"
                    f"{result['strand']}\t{result['phase']}\t"
                    f"{';'.join(attributes)}\n")
