#!/bin/bash
eval "$(conda shell.bash hook)"


#Run all genomes at once
genomes=("Escherichia_coli_K12" "Vibrio_cholerae_N16961" "Wolbachia")

bold=$(tput bold)
normal=$(tput sgr0)

if [[ "$#" -eq 1 ]]; then
	for genome in "${genomes[@]}"; do 
		command="./identify_homologs.sh $1 ${genome}.fna ${genome}.bed ${genomes}_output.txt"

		echo "Identifying which genes in BED file contain the identified homologous histidine kinase domains for ${bold}$genome${normal}. Output File Name: ${genome}_output.txt"

		current_genome_fna="${genome}.fna"
		current_genome_bed="${genome}.bed"
		current_genome_output="${genome}_output.txt"




		#find genes
		temp_hom="${genome}_temp_hom.txt"
		nonuniq_genes="${genome}_all_genes.fa"

		#usage: find_perfect_matches.sh <query file><subject file><output file>blastn_ex2.fna

		
		conda activate week2ex 

		(tblastn -query $1 -subject $current_genome_fna -task tblastn-fast -outfmt '6 std qlen') | awk '$3>30.000 && $4>0.9*$13' > $temp_hom

		#loop through BED file 
		while read -r seqID start stop gene dot strand; do

			has_domain=false
			#loop through temp file to see overlap between start and ends ... think i need to awk for something... we want the whole domain inside gene
			while read -r qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore; do 
	        #echo "Checking gene: $gene, start: $start, stop: $stop, sstart: $sstart, send: $send"

			#check if sequence is within gene 
			if [[ "$start" -le "$sstart" && "$stop" -ge "$send" ]]; then 
						#echo "start: $start sstart: $sstart end: $end send: $send"
				has_domain=true
				break 
			fi 

			done < "$temp_hom"

				#add genes to temp file 
			if [[ "$has_domain" == true ]]; then
					#echo "Appending gene $gene to $nonuniq_genes"
				echo "$gene" >> "$nonuniq_genes"
			fi

				#wc -l < "$4"
		done < "$current_genome_bed"

		#ensure only unique gene values in output
		sort "$nonuniq_genes" | uniq > "$current_genome_output"

		#echo number of genes found
		echo "	Total Genes Identified for ${bold}$genome${normal}: $(wc -l < "$current_genome_output")"
		echo


	done
elif [[ $# -lt 1 ]]; then
	echo 'Too few arguments. Please enter one argument for query file.'
else
	echo 'Only one file can be provided to each argument.'

fi


#col 1 of BED files = seqID corresponding to FASTA headers in associated assembly 

#col 2 & 3 = start (contains location of the boundary of the feature closest to the first base of the FASTA file... significant in col 6) and stop position in sequence where gene is encoded

#col 4 = name or label for feature 

#col 5 = score ... arbitrary "." character for all genes

#col 6 = strand (which of 2 DNA molecules contains the feauture) on which feature is present; implied directionality (+ -> transcription left to right)
	#for all genes with (-) transcription BEGINS at the stop end of the BED file and ENDS at the start column 


##Blast output columns
	 #  1.  qseqid      query or source (gene) sequence id
	 #  2.  sseqid      subject or target (reference genome) sequence id
	 #  3.  pident      percentage of identical positions
	 #  4.  length      alignment length (sequence overlap)
	 #  5.  mismatch    number of mismatches
	 #  6.  gapopen     number of gap openings
	 #  7.  qstart      start of alignment in query
	 #  8.  qend        end of alignment in query
	 #  9.  sstart      start of alignment in subject
	 #  10.  send        end of alignment in subject
	 #  11.  evalue      expect value
	 #  12.  bitscore    bit score