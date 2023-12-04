#!/usr/bin/env python

import sys 


blast_output=sys.argv[1]             #user input blast output file
bed_file=sys.argv[2]     
file_name = bed_file[:-4]   		 #user input bed file
output_file=sys.argv[3]              #output file with matches

blast_output = open(blast_output,'r')
blast_list = blast_output.readlines()

bed_file = open(bed_file,'r')
bed_list = bed_file.readlines()

blast_list_filtered = []

for line in blast_list:
    cols = line.split()

    if len(cols) >= 12:
        identity = float(cols[2])
        length = float(cols[3])
        # query_alignment = float(cols[7]) - float(cols[8])  
        query_alignment = float(cols[12])

        if identity > 30 and length >= 0.9*query_alignment:
            blast_list_filtered.append(line)

found_hom = []

for line in blast_list_filtered:
    columns = line.split()
    homolog_start = float(columns[8])
    homolog_stop = float(columns[9])
    #print(homolog_stop,homolog_start)

    for entry in bed_list:
        line_bed = entry.split()

        locus_start = float(line_bed[1])
        locus_end = float(line_bed[2])
        hom_name = line_bed[3]

        if homolog_start >= locus_start and homolog_stop <= locus_end:
            found_hom.append(hom_name)
            #count += 1
        # elif locus_start >= homolog_stop:
            # break 

found_hom_uniq = []
count = 0
for name in found_hom:
	if name not in found_hom_uniq:
		found_hom_uniq.append(name)
		count += 1

output_file = open(output_file, 'w')
for gene in found_hom_uniq:
	output_file.write(gene+'\n')


blast_output.close()
bed_file.close()
output_file.close()

print(len(found_hom_uniq), 'matches were identified for ', file_name+".")


