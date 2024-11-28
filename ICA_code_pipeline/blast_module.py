import os
from dataclasses import dataclass
from typing import List, Tuple
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbiblastnCommandline, NcbimakeblastdbCommandline
import tempfile
from data.scripts.ICA_code_pipeline.GFF_entry_class import GFFEntry
from dataclasses import dataclass
# 解析GFF属性
def parse_gff_attributes(attr_string: str) -> dict:
    attrs = {}
    for attr in attr_string.split(';'):
        if '=' in attr:
            key, value = attr.split('=', 1)
            attrs[key] = value
    return attrs

# 运行BLASTN并返回匹配的GFF条目
def run_blastn_search(
    target_sequence_file: str,
    genome_fasta_path: str,
    evalue_threshold: float = 1e-5,
    identity_threshold: float = 60.0,
    alignment_length_threshold: int = 30,
    max_target_seqs: int = 10000,
    max_hsps: int = 10,
    threads: int = 16
) -> List[GFFEntry]:

    # 定义BLAST数据库路径
    genome_blast_db_path = os.path.splitext(genome_fasta_path)[0]

    # 构建BLAST数据库
    if not os.path.exists(genome_blast_db_path + ".nin"):
        makeblastdb_cmd = NcbimakeblastdbCommandline(
            dbtype='nucl',
            input_file=genome_fasta_path,
            out=genome_blast_db_path
        )
        makeblastdb_cmd()

    # 运行BLASTN
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.outfmt6') as temp_output:
        blast_output_path = temp_output.name

    blastn_cmd = NcbiblastnCommandline(
        query=target_sequence_file,
        db=genome_blast_db_path,
        evalue=evalue_threshold,
        outfmt='6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore',
        max_target_seqs=max_target_seqs,
        max_hsps=max_hsps,
        num_threads=threads,
        out=blast_output_path,
        perc_identity=identity_threshold
    )
    blastn_cmd()

    # 解析BLAST结果
    gff_results = []
    with open(blast_output_path, 'r') as blast_file:
        for line in blast_file:
            fields = line.strip().split('\t')
            qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = fields
            if (float(evalue) <= evalue_threshold and
                float(pident) >= identity_threshold and
                int(length) >= alignment_length_threshold):
                
                # 根据 sstart 和 send 的相对大小判断链的方向
                sstart = int(sstart)
                send = int(send)
                if sstart < send:
                    strand = '+'
                else:
                    strand = '-'
                
                gff_entry = GFFEntry(
                    seqid=sseqid,
                    source='BLASTN',
                    type='transposable_element',
                    start=min(sstart, send),  # 确保 start 是较小值
                    end=max(sstart, send),    # 确保 end 是较大值
                    score=float(evalue),
                    strand=strand,  # 动态设置链方向
                    phase='.',
                    attributes={'Name': 'Transposable_Element'}
                )
                gff_results.append(gff_entry)

        return gff_results

def extract_transposable_elements(
    gff_entries: List[GFFEntry],
    genome_fasta_path: str
) -> List[SeqRecord]:
    genome_dict = SeqIO.to_dict(SeqIO.parse(genome_fasta_path, "fasta"))
    transposable_elements = []

    for entry in gff_entries:
        chrom = entry.seqid
        start = entry.start - 1500
        end = entry.end + 1500
        # safety checks
        if start < 0:
            start = 0
        if end > len(genome_dict[chrom].seq):
            end = len(genome_dict[chrom].seq)

        te_sequence = genome_dict[chrom].seq[start-1:end]

        # 修改 ID 和 description
        te_id = f"TE_{chrom}_{start}_{end}"
        te_description = f"{chrom}:{start}-{end}"

        te_record = SeqRecord(
            te_sequence,
            id=te_id,
            description=te_description
        )

        # 添加位置信息到注释中
        te_record.annotations['blast_match'] = {
            'chromosome': chrom,
            'start': start,
            'end': end,
            'original_match_start': entry.start,
            'original_match_end': entry.end,
            'strand': entry.strand  # 保存 BLAST 的链方向
        }

        transposable_elements.append(te_record)

    return transposable_elements
