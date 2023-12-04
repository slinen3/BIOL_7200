#!/usr/bin/env python3
# ispcr/blast.py

import subprocess
from typing import List, Tuple
import tempfile
from subprocess import PIPE

def BLAST_function(primer_file: str, assembly_file: str, output_file: str):

	""" Takes in specified primer, assembly, and output files and runs BLAST on them, outputting to the specified output file.

	Args:
		primer_file (str): Name of primer file (FASTA) to be analyzed.
		assembly_file (str): Name of assembly file (FASTA) to be analyzed.
		output_file (str): Desired name of output file for BLAST results.


	Returns: 
		None.

	Example Usage:
		BLAST_function(primer_file, assembly_file, 'blast_output.txt')
	"""

	command = ['blastn', '-query', primer_file,
		'-subject', assembly_file,
		'-task','blastn-short',
		'-outfmt','6 std qlen'
		]

	with open(output_file, 'w') as output:
		subprocess.run(command, stdout=output)

def filter_blast(blast_output: str, filtered_output: str):

	""" Takes in output from BLAST function (or any BLAST outputfile) and filters it to yield only hits that are full length and with at least (>=) 80 percent identity.

	Args:
		blast_output (str): Name of BLAST file to be analyzed (this comes from step_one()).
		filtered_output (str): Desired name of file to output filtered results to. 

	Returns: 
		None.

	Example Usage:
		filter_blast('blast_output.txt', filtered_blast)

	"""

	awk_filter = [ 'awk','{if ($13==$4 && $3>=80) print}'] #alignment positions must match and pident >= 80

	with open(filtered_output,'w') as filtered:
		subprocess.run(awk_filter, stdin=open(blast_output,'r'),stdout=filtered)


def step_one(primer_file: str, assembly_file: str) -> List[List[str]]:

	""" Takes in primer file and assembly file, runs BLAST on them, and filters blast results using the functions BLAST_function and filter_blast functions. 
	Args:
		primer_file (str): Name of primer file (FASTA) to be analyzed.
		assembly_file (str): Name of assembly file (FASTA) to be analyzed.

	Returns: 
		result (List[List[str]]): List of unique, filtered BLAST hits sorted by ascending matched position (5').

	Example Usage:
		sorted_good_hits = step_one(
		primer_file="data/general_16S_515f_806r.fna",
		assembly_file="data/Vibrio_cholerae_N16961.fna"
		)

		for hit in sorted_good_hits:
			print(hit)

	"""

	BLAST_function(primer_file, assembly_file, 'blast_output.txt')
	filtered_blast = 'filtered_blast_output.txt'
	filter_blast('blast_output.txt', filtered_blast)

	unique_entries = set()
	result = []

	with open(filtered_blast, 'r') as filtered:
		for line in filtered:
			entry = line.strip().split('\t')
			entry_tuple = tuple(entry)

			if entry_tuple not in unique_entries:
				unique_entries.add(entry_tuple)
				result.append(entry)

	result.sort(key=lambda x: int(x[8]))
	return result


def step_two(sorted_hits: List[str], max_amplicon_size: int) -> List[Tuple[List[str]]]:

	""" Takes in sorted hits (output from step_one) and max amplicon size (), matches hits that are less than max amplicon size apart and pointing towards eachother.
	Args:
		sorted_hits (List[str]): List of sorted BLAST hits from step_one().
		max_amplicon_size (int): Maximum distance allowed between hits to form an amplicon.

	Returns: 
		result (List[Tuple[List[str]]]): List of tuples with pairs of hits forming amplicons based on the criteria specified.

	Example Usage:
		paired_hits = step_two(
		sorted_hits=sorted_good_hits,
		max_amplicon_size=1000
		)

	for hit in paired_hits:
		print(hit)

	"""
	# result = []
	# for i in range(len(sorted_hits) - 1):
	# 	current_hit = sorted_hits[i]
	# 	three_prime_current = int(current_hit[9])
	# 	five_prime_current = int(current_hit[8])

	# 	for j in range(i + 1, len(sorted_hits)):
	# 		following_hit = sorted_hits[j]
	# 		three_prime_next = int(following_hit[9])
	# 		five_prime_next = int(following_hit[8])

	# 		if (five_prime_current < five_prime_next and abs(three_prime_next - three_prime_current) < max_amplicon_size):
	# 			result.append((current_hit, following_hit))

	result = [(sorted_hits[i], sorted_hits[j]) for i in range(len(sorted_hits) - 1) for j in range(i + 1, len(sorted_hits)) if sorted_hits[i][8] < sorted_hits[j][8] and abs(sorted_hits[j][8]-sorted_hits[i][9]) < max_amplicon_size]

	return result


def step_three(hit_pairs: List[Tuple[List[str]]], assembly_file: str) -> str:

	""" Takes in hit pairs (output from step_two) and assembly file, converts hit pairs from a list of tuples to a BED file with the contig, start, and stop positions, and runs seqtk using the BED file generated and assembly file provided.

	Args:
		hit_pairs (List[Tuple[List[str]]]): List of tuples indicating pairs of hits forming amplicons.
		assembly_file (str): Name of assembly file (FASTA) to be analyzed.

	Returns: 
		None.

	Example Usage:
		amplicons = step_three(
		hit_pairs=hit_pairs,
		assembly_file="data/Vibrio_cholerae_N16961.fna"
		)

		print(amplicons)

	"""

	bed_file_path = 'bed_file.bed'
	count = 0


	with open(bed_file_path,'w') as bed_file:
		for i in hit_pairs:
			contig = i[0][1]
			start = i[0][9]
			end = int(i[1][9])-1
			count+=1

			bed_file.write(f"{contig}\t{start}\t{end}\n")

			# bed_file.write(f"{contig}\t{start}\t{end}")
			# if count =len(hit_pairs):
			# 	bed_file.write('\n')
			# count += 1

	command = ['seqtk', 'subseq', assembly_file, bed_file_path]

	result = subprocess.run(command,check=True,stdout=PIPE,text=True)

	return result.stdout.strip()



if __name__ == "__main__":

	result = step_one(primer_file,assembly_file)
	print(result)






