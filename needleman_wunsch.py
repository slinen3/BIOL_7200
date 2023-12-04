import subprocess
import os
import sys
import tempfile
import numpy as np
from collections import defaultdict
from typing import List, Tuple

def needleman_wunsch(seq_a: str, seq_b: str, match: int, mismatch: int, gap: int) -> Tuple[Tuple[str,str],int]:

	""" Takes in two sequences, the match score, mismatch score, and gap score. It aligns the two sequences based on the needleman wunsch algorithm and returns a tuple of aligned sequences with the corresponding score. 

	Args:
		seq_a (str): String of the first sequence to be analyzed.
		seq_a (str): String of the second sequence to be analyzed.
		match (int): Int representing the score assigned when bases in the same position in the input sequences match. 
		mismatch (int): Int representing the score assigned when bases in the same position in the input sequences do not match.
		gap (int): Int representing the score assigned when an insertion or deletion (gap) is introduced. This aligns sequences of different lengths.

	Returns: 
		((aligned_1, aligned_2), score) Tuple[Tuple[str,str],int]: Tuple containing tuple of the aligned sequences and the alignment score. 

	Example Usage:
		magnumopus.needleman_wunsch(seq1, seq2, match, mismatch, gap)
	"""

	row = len(seq_a) + 1
	col = len(seq_b) + 1

	matrix = np.zeros((row,col))

	for i in range(0, row):
		matrix[i][0] = i * gap
	for j in range(0, col):
		matrix[0][j] = j * gap

	for k in range(1, row):
		for l in range(1, col):
			sub_score = match if seq_a[k-1] == seq_b[l -1] else mismatch
			match_score = matrix[k-1][l-1] + sub_score
			remove = matrix[k-1][l] + gap
			insert = matrix[k][l-1] + gap
			matrix[k][l] = max(match_score, remove, insert)

	aligned_a, aligned_b = '', ''
	len_a = len(seq_a)
	len_b = len(seq_b)

	while len_a > 0 or len_b > 0:
		#diagonal
		if len_a > 0 and len_b > 0 and matrix[len_a][len_b] == matrix[len_a - 1][len_b - 1] + (1 if seq_a[len_a - 1] == seq_b[len_b - 1] else -1):
			aligned_a += seq_a[len_a - 1]
			aligned_b += seq_b[len_b - 1]
			len_a -= 1
			len_b -= 1
		#vertical
		elif len_a > 0 and matrix[len_a][len_b] == matrix[len_a - 1][len_b] + gap:
			aligned_a += seq_a[len_a - 1]
			aligned_b += '-'
			len_a -= 1
		#horizontal
		elif j > 0 and matrix[len_a][len_b] == matrix[len_a][len_b - 1] + gap:
			aligned_a += '-'
			aligned_b += seq_b[len_b - 1]
			len_b -= 1


	aligned_1 = aligned_a[::-1]
	aligned_2 = aligned_b[::-1]
	score = matrix[row -1][col - 1]



	return (aligned_1, aligned_2), score