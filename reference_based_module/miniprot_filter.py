import os
import subprocess
from dataclasses import dataclass
from typing import List, Dict, Tuple
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline, NcbimakeblastdbCommandline
import tempfile
import os


from GFFEntry_class import GFFEntry, parse_gff_attributes

import os
import subprocess
from dataclasses import dataclass
from typing import List, Dict, Tuple
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline, NcbimakeblastdbCommandline
import tempfile
import os


from GFFEntry_class import GFFEntry, parse_gff_attributes

# This section is the implementation of the miniprot analysis, to identify protein region (transposase)
def run_miniprot_analysis(
    miniprot_path: str,
    Tpase: str,
    transposable_elements: List[SeqRecord],
    threads: int = 16
) -> List[GFFEntry]:
    """
    Run the miniprot analysis to identify protein regions in the transposable elements.

    Args:
        miniprot_path (str): Path to the miniprot executable.
        Tpase (str): Path to the Tpase sequences file.
        transposable_elements (list): List of SeqRecord objects representing the transposable elements.
        threads (int): Number of threads to use.
    
    Returns:
        list: List of GFFEntry objects representing the protein regions.
    """
    try:
        miniprot_results = []

        # Create a temporary FASTA file to store all transposable elements
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fasta') as all_te_fasta:
            SeqIO.write(transposable_elements, all_te_fasta, "fasta")
            all_te_fasta_path = all_te_fasta.name

        # Create a temporary file to store the miniprot output in GFF format
        miniprot_output = tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.gff').name

        # run miniprot command, with default parameters
        miniprot_cmd = f"{miniprot_path} -ut {threads} -c 50000 -m 10 -p 0.2 -N 200 -O 3 -J 8 -F 8 -K 5M --outn=5000 --outs=0.5 --outc=0.03 --gff {all_te_fasta_path} {Tpase} > {miniprot_output}"
        result = subprocess.run(miniprot_cmd, shell=True, capture_output=True, text=True)

        if result.returncode != 0:
            print(f"miniprot run with error: {result.returncode}")
            return []

        # parse the miniprot output in GFF format
        with open(miniprot_output, 'r') as miniprot_file:
            for line in miniprot_file:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue
                seqid, source, feature, start, end, score, strand, phase, attributes = fields
                if feature == 'CDS':
                    # Search for the corresponding transposable element record
                    te_record = next((te for te in transposable_elements if te.id == seqid), None)

                    if te_record is None:
                        print(f"Warming: No TE record found for {seqid}")
                        continue

                    # Get the original BLAST match information
                    blast_strand = te_record.annotations['blast_match']['strand']

                    # If the BLAST strand is minus, reverse the strand
                    if blast_strand == '-':
                        if strand == '+':
                            strand = '-'
                        elif strand == '-':
                            strand = '+'

                    gff_entry = GFFEntry(
                        seqid=seqid,  # use the TE ID
                        source='miniprot',
                        type='protein_coding_region',
                        start=int(start),
                        end=int(end),
                        score=float(score) if score != '.' else 0.0,
                        strand=strand, 
                        phase=phase,
                        attributes=parse_gff_attributes(attributes)
                    )
                    miniprot_results.append(gff_entry)

        # clean up the temporary files
        os.remove(all_te_fasta_path)
        os.remove(miniprot_output)

        return miniprot_results

    except Exception as e:
        print(f"Miniprot runtime error: {str(e)}")
        return []
