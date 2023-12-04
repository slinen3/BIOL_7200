#!/usr/bin/env python3

import sys 


if len(sys.argv) != 2: 
	print("Execute in this format: script_name.py <FASTA file>.")
	sys.exit()
input_file = sys.argv[1]	

input_file_opened = open(input_file,'r')
line_list = input_file_opened.readlines()

input_file_opened.close()


seqs = []
current_seq = ""

for line in line_list:
	line = line.strip()
	if line[0] == ">":
		if current_seq:
			seqs.append(current_seq)
		current_seq = ""
	else:
		current_seq += line

if current_seq:
	seqs.append(current_seq)

# if len(seqs < 2):

new_file = open("matched.txt","w+")
new_file.write(seqs[0]+'\n')
print(seqs[0])

match = ""

for i in range(len(seqs[0])):
	for j in range(1, len(seqs)):
		if seqs[j][i] == seqs[0][i]:
			match += "|"

		else:
			match += " "
	# for sequence in seqs:
	# 	print(separator, end="")

new_file.write(match+'\n')
print(match)
new_file.write(seqs[1])
print(seqs[1])



input_file_opened.close()
new_file.close()
	
