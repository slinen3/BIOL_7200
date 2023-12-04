#!/usr/bin/env python3
import sys

def check_parentheses(test_str: str):
	""" Takes in a string of parentheses and checks if they are properly paired (i.e. every opening parenthese '(' has a complementary closing parenthese ')' that comes after it at some point.).

	Args:
		test_str (str): Stirng containing parentheses to be checked. 

	Prints: 
		"PAIRED" if parentheses all have a proper pair
		"NOT PAIRED" if parentheses are not properly paired

	Example Usage:
		>>> check_parentheses("(()())))()))((((")
		NOT PAIRED 

		>>> check_parentheses(")")
		NOT PAIRED 

		>>> check_parentheses("((()))")
		PAIRED
	"""
	test_str = sys.argv[1]

	open_list = ["("]
	close_list = [")"]
	stack = []

	for char in test_str:
		if char in open_list:
			stack.append(char)
		elif char in close_list:
			posit = close_list.index(char)
			if ((len(stack) > 0) and (open_list[posit] == stack[len(stack)-1])):
				stack.pop()
			else:
				print("NOT PAIRED")
				sys.exit()
	if len(stack) == 0:
		print("PAIRED")
	else:
		print("NOT PAIRED")





if __name__ == '__main__':
	check_parentheses(sys.argv[1])