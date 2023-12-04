#!/usr/bin/env python3

import magnumopus
import argparse
import os


def amplicon_aligned(assembly_file1: str, assembly_file2: str, primer_file: str, max_amplicon_size: int, match: int, mismatch: int, gap: int):

	""" Takes in two assembly files, a primer file, the max amplicon size, the match score, mismatch score, and gap score. It performs a BLAST on each of the assembly files to find the hits, finds the hit pairs, and extracts the amplicons for these. It then aligns the two resulting sequences based on the needleman wunsch algorithm and returns the aligned sequences with the corresponding score. 

	Args:
		assembly_file1 (str): First assembly file to be analyzed (FASTA).
		assembly_file2 (str): Second assembly file to be analyzed (FASTA).
		primer_file (str): Primer  file to use in BLAST against the assembly files provided.
		max_amplicon_size (int): Int representing the maximum allowed size for an amplicon.
		match (int): Int representing the score assigned when bases in the same position in the input sequences match. 
		mismatch (int): Int representing the score assigned when bases in the same position in the input sequences do not match.
		gap (int): Int representing the score assigned when an insertion or deletion (gap) is introduced. This aligns sequences of different lengths.

	Returns: 
		None.

	Example Usage:
		For help:
			./amplicon_align.py -h

				usage: amplicon_align.py [-h] -1 ASSEMBLY1 -2 ASSEMBLY2 -p PRIMERS -m
				MAX_AMPLICON_SIZE --match MATCH --mismatch MISMATCH --gap GAP
				Perform in-silico PCR on two assemblies and align the amplicons
				options:
				-h, --help show this help message and exit
				-1 ASSEMBLY1 Path to the first assembly file
				-2 ASSEMBLY2 Path to the second assembly file
				-p PRIMERS Path to the primer file
				-m MAX_AMPLICON_SIZE maximum amplicon size for isPCR
				--match MATCH match score to use in alignment
				--mismatch MISMATCH mismatch penalty to use in alignment
				--gap GAP gap penalty to use in alignment

		For results:
			/amplicon_align.py -1 data/Pseudomonas_aeruginosa_PAO1.fna -2 
			data/Pseudomonas_protegens_CHA0.fna -p data/rpoD.fna -m 2000 
			--match 1 --mismatch=-1 --gap=-1
	"""


    for i in range(1, 3):
        max_amp_size = max_amplicon_size
        if i == 1:
            assembly = assembly_file1
            amplicons_assembly1 = magnumopus.ispcr(primer_file, assembly, max_amp_size)
            lines = amplicons_assembly1.split('\n')
            seq_lines = [line for line in lines if not line.startswith('>')]
            amplicons_assembly1 = ''.join(seq_lines)


        elif i == 2:
            assembly = assembly_file2
            amplicons_assembly2 = magnumopus.ispcr(primer_file, assembly, max_amp_size)
            lines = amplicons_assembly2.split('\n')
            seq_lines = [line for line in lines if not line.startswith('>')]
            amplicons_assembly2 = ''.join(seq_lines)



    # seq1 = amplicons_assembly1
    # seq2 = amplicons_assembly2

    seq1 = amplicons_assembly1.replace('-', '')
    seq2 = amplicons_assembly2.replace('-', '')
 

    aln, score = magnumopus.needleman_wunsch(seq1, seq2, match, mismatch, gap)

     # reverse complement seq2 and calculate alignment score
    #rev_seq2 = reverse_complement(seq2)

    complements = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    seq = read_fasta(assembly_file2)
    rev_seq = seq[::-1]
    rev_comp = ''.join(complements.get(base, base) for base in rev_seq)

    rev_assembly_file2_path = 'data/rev_' + os.path.basename(assembly_file2)


    write_fasta(rev_assembly_file2_path, rev_comp, header="Reverse Assembly 2")



    seq2_rev = magnumopus.ispcr(primer_file, rev_assembly_file2_path, max_amplicon_size)
    lines = seq2_rev.split('\n')
    seq_lines = [line for line in lines if not line.startswith('>')]
    seq2_rev = ''.join(seq_lines)



    aln_rev, score_rev = magnumopus.needleman_wunsch(seq1, seq2_rev, match, mismatch, gap)


    if score >= score_rev:
        best_aln = aln
        best_score = score
    else:
        best_aln = aln_rev
        best_score = score_rev

    print(f"{best_aln[0]}\n{best_aln[1]}\n{best_score}")



def read_fasta(file_path: str):

	""" Takes in a file path of the desired FASTA to be read and returns the sequence without header. 
	Args:
		file_path (str): Path of FASTA file to be read.
		
	Returns: 
		sequence (str): String of the sequence within the FASTA file -- no headers 

	Example Usage:
		read_fasta(assembly_file2)
	"""


	sequence = ""
	with open(file_path, 'r') as file:
		for line in file:
			if not line.startswith('>'):
				sequence += line.strip()
	return sequence


def write_fasta(file_path: str, sequence: str, header: str):

	""" Takes in a file path of the desired output name for FASTA to be written, sequence of interest, and header and writes the header and sequence to the FASTA file.
	Args:
		file_path (str): Output path of FASTA file.
		sequence (str): Sequence to write to FASTA file.
		header (str): Header to write to FASTA file. 
		
	Returns: 
		None.

	Example Usage:
		write_fasta(rev_assembly_file2_path, rev_comp, header="Reverse Assembly 2")

	"""

	with open(file_path,'w') as file:
		if header:
			file.write(f">{header}\n")
		for i in range(0, len(sequence), 60):
			file.write(sequence[i:i+60]+'\n')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        usage="%(prog)s -1 ASSEMBLY1 -2 ASSEMBLY2 -p PRIMERS -m MAX_AMPLICON_SIZE --match MATCH --mismatch MISMATCH --gap GAP",
        description="Perform in-silico PCR on two assemblies and align the amplicons."
    )

    parser.add_argument("-1", "--ASSEMBLY1", help="Path to the first assembly file")
    parser.add_argument("-2", "--ASSEMBLY2", help="Path to the second assembly file")
    parser.add_argument("-p", "--PRIMERS", help="Path to the primer file")
    parser.add_argument("-m", "--MAX_AMPLICON_SIZE", help="maximum amplicon size for isPCR")
    parser.add_argument("--match", help="match score to use in alignment", type=int)
    parser.add_argument("--mismatch", help="mismatch penalty to use in alignment", type=int)
    parser.add_argument("--gap", help="gap penalty to use in alignment", type=int)

    args = parser.parse_args()

    amplicon_aligned(
        args.ASSEMBLY1,
        args.ASSEMBLY2,
        args.PRIMERS,
        int(args.MAX_AMPLICON_SIZE),
        int(args.match),
        int(args.mismatch),
        int(args.gap),
    )
