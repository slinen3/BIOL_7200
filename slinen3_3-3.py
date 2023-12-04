#!/usr/bin/env python3
import sys
import pandas as pd


def read_bedfile(bedfile: str):
	""" Reads data from the .bed file passed into the terminal and organizes it into a DataFrame using pandas. 

	Args:
		bedfile (str): Path to .bed file.

	Returns: 
		DataFrame containing parsed data with columns ['Chromosome', 'Start', 'Stop', 'Name','.','Direction']. 'Name' is the index of the DataFrame.

	Example Usage:
		>>> bed_file = Vibrio_cholerae_N16961.bed
		>>> bed_data = read_bedfile(bed_file)
		>>> print(bed_data)
				                   Chromosome    Start     Stop  . Direction
		Name                                                        
		dnaN            NZ_CP028827.1     1462     2559  .         +
		recF            NZ_CP028827.1     2580     3668  .         +
		gyrB            NZ_CP028827.1     3692     6106  .         +
		N16961_RS00025  NZ_CP028827.1     6534     6953  .         +
		N16961_RS00035  NZ_CP028827.1     7358     7792  .         +
		...                       ...      ...      ... ..       ...
		N16961_RS19570  NZ_CP028828.1  1065954  1066634  .         -
		N16961_RS19575  NZ_CP028828.1  1066680  1068200  .         -
		N16961_RS19580  NZ_CP028828.1  1068214  1068813  .         -
		N16961_RS19585  NZ_CP028828.1  1068990  1071020  .         -
		N16961_RS19590  NZ_CP028828.1  1071241  1072269  .         +

		[3647 rows x 5 columns]
	"""

	bed_data = pd.read_csv(bedfile,sep='\t')
	headers = ['Chromosome', 'Start', 'Stop', 'Name','.','Direction']
	bed_data.columns = headers[:len(bed_data.columns)]

	bed_data.set_index('Name', inplace=True)

	return bed_data


def read_fasta(assemblyfile: str):
	""" Reads data from the .fna file passed into the terminal and organizes the sequence for each chromosome into a single string. 

	Args:
		fasta (.fna) file (str): Path to fasta file.

	Returns:
		seq_str_chr1, seq_str_chr2 (str): String(s) of sequence correlated to each chromosome for the assembly file.

	Note:
		Designed for the Vibrio_cholerae_N16961.fna fasta file specifically (has 2 chromosomes). To modify you will need to update the specific chromosomes perfinent to your file.
	
	Example:
		>>> fasta_file = Vibrio_cholerae_N16961.fna
		>>> seq1, seq2 = read_fasta(fasta_file)
		>>> print(len(seq1),len(seq2))
	"""

	file = open(assemblyfile)
	seq_str_chr1 = ''
	seq_str_chr2 = ''
	for line in file:
		if ">" in line:
			chr_name = line[1:14] 
		elif chr_name == 'NZ_CP028827.1':
			seq_str_chr1 += line
		else:
			seq_str_chr2 += line

	seq_str_chr1 = seq_str_chr1.replace('\n','')
	seq_str_chr2 = seq_str_chr2.replace('\n','')

	return(seq_str_chr1, seq_str_chr2)



def find_homologs(blastfile: str,bedfile: str,assemblyfile: str,outputfile: str):

	""" Takes in blastfile, bedfile, fatsafile, and name of a designated output file and:
		 - processes the BLAST output to keep only hits with >30% identity and 90% length 
		 - identifies BED features (i.e., genes) that have a BLAST hit entirely within the boundaries of the feature start and end (same as last exercise)
 		 - extracts the sequence of identified homologous genes from the assembly sequence
 		 - if the gene is encoded on the - strand, reverse complements the sequence
 		 - writes the sequences of the homologous genes to the specified output file in FASTA format (gene name as header).

	Args:
		 - blastfile (str): Path to BLAST file.
		 - bedfile (str): Path to .bed file.
		 - assemblyfile (str): Path to .fna file.
		 - outputfile (str): Desired path for output (.fna) file.

	Returns: 
		none

	Note:
		works for the Vibrio_cholerae_N16961.fna fasta file specifically. To modify you will need to update the specific chromosomes perfinent to your file.

	Example Usage:
		>>> find_homologs("blast_output.txt", "genes.bed", "assembly.fna", "output.fna")
	"""


	blastfile = sys.argv[1]
	bedfile = sys.argv[2]
	bed_data = read_bedfile(bedfile)

	assemblyfile = sys.argv[3]
	seq_str_chr1, seq_str_chr2 = read_fasta(assemblyfile)

	outputfile = sys.argv[4]


	# read blast file
	hits = []
	with open(blastfile) as fin:
	    for line in fin: # .readlines() is the default iter method for the open file class
	        
	        # unpack and convert types of desired columns. This is ugly. We'll revisit later...
	        _, sid, pcnt, matchlen, _, _, _, _, sstart, send, _, _, qlen = line.split()
	        pcnt = float(pcnt)
	        matchlen = int(matchlen)
	        sstart = int(sstart)
	        send = int(send)
	        qlen = int(qlen)
	    
	        # Keep hits that could be homologs
	        if pcnt > 30 and matchlen > 0.9*qlen:
	            # We could store matches as a list or tuple.
	            # We won't want to modify the elements so a tuple is "safer" in that we then can't modify it by mistake
	            hits.append((sid, sstart, send))

	# Now read the bed file
	feats = []
	with open(bedfile) as fin:
	    for line in fin:
	        bed_sid, bed_start, bed_end, gene, *_ = line.split() # an asterisk before a variable name when unpacking makes that variable store remaining elements
	        bed_start = int(bed_start)
	        bed_end = int(bed_end)
	        
	        feats.append((bed_sid, bed_start, bed_end, gene))

	# Now we have our two datasets read in, we can loop over them to find matches
	homologs = []
	for blast_sid, blast_sstart, blast_send in hits: # unpack our blast data
	    for bed_sid, bed_start, bed_end, gene in feats:
	        # Don't bother checking the rest if the sid doens't match
	        if blast_sid != bed_sid:
	            continue
	        
	        # Once we are dealing with features at higher index locations than our hit, go to the next hit (break loop over feats)
	        if blast_sstart <= bed_start or blast_send <= bed_start:
	            break
	        
	        # Otherwise, check if the hit is inside the feature
	        if (blast_sstart > bed_start
	            and blast_sstart <= bed_end
	            and blast_send > bed_start
	            and blast_send <= bed_end
	        ):
	            homologs.append(gene)
	            break # Each BLAST hit will only be in one feature so move to next hit once you've found it

	# Get the unique homologs using a set()
	unique_homologs = set(homologs)


	#extract sequence of identified homologous genes from assembly sequence and reverse complement if gene is encoded on the '-' strand
	seq_homolog_genes = dict()
	output_contents = ''
	for homolog in unique_homologs:
		bed_row = bed_data.loc[homolog]
		direction = bed_row['Direction']
		chrom = bed_row['Chromosome']
		if chrom == 'NZ_CP028827.1': #checks which chromosome gene is in (there are 2 in this file)
			bed_row = bed_data.loc[homolog]
			seq_start = bed_row['Start']
			seq_end = bed_row['Stop']
			direction = bed_row['Direction']
	


			gene_seq = seq_str_chr1[seq_start:seq_end+1]
			gene_seq = gene_seq.replace('\n','')


			if direction == '-':
				gene_seq = seq_str_chr1[seq_start-1:seq_end]
				gene_seq = gene_seq.replace('\n','')
				rev_comp_seq = reverse_complement(gene_seq)
				seq_homolog_genes[homolog] = rev_comp_seq
			else:
				seq_homolog_genes[homolog] = gene_seq
		else:
			bed_row = bed_data.loc[homolog]
			seq_start = bed_row['Start']
			seq_end = bed_row['Stop']
			direction = bed_row['Direction']

			gene_seq = seq_str_chr2[seq_start:seq_end+1]
			gene_seq = gene_seq.replace('\n','')


			if direction == '-':
				gene_seq = seq_str_chr2[seq_start-1:seq_end]
				gene_seq = gene_seq.replace('\n','')
				rev_comp_seq = reverse_complement(gene_seq)
				seq_homolog_genes[homolog] = rev_comp_seq

			else:
				seq_homolog_genes[homolog] = gene_seq

	#Write output to designated file
	write_seqs(outputfile, seq_homolog_genes)



def reverse_complement(gene_seq: str):
	""" Creates the reverse complement of DNA sequence input.

	Args:
		gene_seq (str): DNA sequence input.

	Returns: 
		rev_comp_seq: string of the reverse complement of the DNA sequence input.

	Example Usage:
		>>> reverse_complement('ACCGTATGGACCTAAG')
		>>> 'CTTAGGTCCATACGGT'
	"""

	complements = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

	bases = [base for base in gene_seq if base in complements]

	reversed_seq = ''.join(bases)[::-1]
	complement = [complements[base] for base in reversed_seq]

	rev_comp_seq = ''.join(complement)
	return rev_comp_seq 


def write_seqs(outputfile: str, seq_homolog_genes: dict):
	""" Writes sequences of homologous genes to the output file specified using command line input in FASTA format.

	Args:
		outputfile (str): specified output file from command line input 
		seq_homolog_genes (dict): dictionary containing all of the genes and their corresponding sequences as keys and values respectively. 

	Returns: 
		none

	Example Usage:
		>>> write_seqs("output_Ecoli.fna", {"gene1":"CCTGA", "gene2":"AAGTCTG"})
		>>> cat output_Ecoli.fna
		>gene1
		CCTGA
		>gene2
		AAGTCTG
	"""
	with open(outputfile,"w") as f:
		for gene in seq_homolog_genes:
			f.write(">"+gene+'\n'+seq_homolog_genes[gene]+'\n')





if __name__ == '__main__':
	find_homologs(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])
		