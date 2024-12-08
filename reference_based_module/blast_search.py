import os
from dataclasses import dataclass
from typing import List, Dict, Tuple
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline, NcbimakeblastdbCommandline
from io import StringIO
import tempfile
import os

from GFFEntry_class import GFFEntry, parse_gff_attributes




# run blastn search, and return a list of GFFEntry objects
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
    """
    Run a BLASTN search and return the results as a list of GFFEntry objects.

    Args:
        target_sequence_file (str): Path to the target sequence file.
        genome_fasta_path (str): Path to the genome FASTA file.
        evalue_threshold (float): E-value threshold for BLASTN.
        identity_threshold (float): Identity threshold for BLASTN.
        alignment_length_threshold (int): Alignment length threshold for BLASTN.
        max_target_seqs (int): Maximum number of target sequences for BLASTN.
        max_hsps (int): Maximum number of HSPs for BLASTN.
        threads (int): Number of threads to use.
    
    Returns:
        list: List of GFFEntry objects representing the BLASTN results.
    """

    # define the path for the BLAST database
    genome_blast_db_path = os.path.splitext(genome_fasta_path)[0]

    # build a BLAST database if it does not exist
    if not os.path.exists(genome_blast_db_path + ".nin"):
        makeblastdb_cmd = NcbimakeblastdbCommandline(
            dbtype='nucl',
            input_file=genome_fasta_path,
            out=genome_blast_db_path
        )
        makeblastdb_cmd()

    # run BLASTN search 
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

    # parse the BLAST output and filter the results
    gff_results = []
    with open(blast_output_path, 'r') as blast_file:
        for line in blast_file:
            fields = line.strip().split('\t')
            # define the required fields
            qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = fields
            if (float(evalue) <= evalue_threshold and
                float(pident) >= identity_threshold and
                int(length) >= alignment_length_threshold):
                
                # judge the strand based on the start and end positions
                sstart = int(sstart)
                send = int(send)
                # if the start position is less than the end position, the strand is '+'
                if sstart < send:
                    strand = '+'
                # if the start position is greater than the end position, the strand is '-'
                else:
                    strand = '-'
                
                # create a GFFEntry object
                gff_entry = GFFEntry(
                    seqid=sseqid,
                    source='BLASTN',
                    type='transposable_element',
                    start=min(sstart, send),  # make sure start is not greater than end
                    end=max(sstart, send),    # make sure end is not less than start
                    score=float(evalue),
                    strand=strand,  # set the strand
                    phase='.',
                    attributes={'Name': 'Transposable_Element'}
                )
                gff_results.append(gff_entry)

        return gff_results


def extract_transposable_elements(gff_entries: List[GFFEntry], genome_fasta_path: str, extension: int = 5000) -> List[SeqRecord]:
    genome_dict = SeqIO.to_dict(SeqIO.parse(genome_fasta_path, "fasta"))
    transposable_elements = []
    unique_sequences = set()  

    for entry in gff_entries:
        chrom = entry.seqid
        extended_start = max(0, entry.start - extension)  
        extended_end = min(len(genome_dict[chrom].seq), entry.end + extension)  

        te_sequence = genome_dict[chrom].seq[extended_start:extended_end]  

        if str(te_sequence) in unique_sequences:  
            continue
        unique_sequences.add(str(te_sequence))  

        te_id = f"TE_{chrom}_{extended_start}_{extended_end}"
        te_description = f"{chrom}:{extended_start}-{extended_end}"

        te_record = SeqRecord(
            te_sequence,
            id=te_id,
            description=te_description
        )

        te_record.annotations['blast_match'] = {
            'chromosome': chrom,
            'start': extended_start,
            'end': extended_end,
            'original_match_start': entry.start,
            'original_match_end': entry.end,
            'strand': entry.strand
        }

        transposable_elements.append(te_record)

    return transposable_elements



