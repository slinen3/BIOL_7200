#!/usr/bin/env python3
import sys

def generate_triangle(char: str, num_lines: int):

	""" Creates a triangle of specified length using a specific character.

	Args:
		char (str): Character used to create the triangle
		num_lines (int): Number of lines in triangle (must be > 0)
		

	Prints: 
		Triangle pattern with specified number of lines of specified char.

	Example Usage:
		>>> generate_triangle(X, 9)
		X
		XX
		XXX
		XXXX
		XXXXX
		XXXX
		XXX
		XX
		X

		>>> generate_triangle($, 12)
		$
		$$
		$$$
		$$$$
		$$$$$
		$$$$$$
		$$$$$$
		$$$$$
		$$$$
		$$$
		$$
		$

	"""
	char = sys.argv[1]
	num_lines = sys.argv[2]
	num_lines = int(num_lines)

	i = 1
	j = num_lines // 2

	for num in range(0,num_lines):
		if num < num_lines/2:
			
			print(char*i)
			i += 1
		else:
			print(char*j)
			j -= 1


if __name__ == '__main__':
	generate_triangle(sys.argv[1],sys.argv[2])