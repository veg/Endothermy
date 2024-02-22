########################################################
### SCRIPT TO LABEL ALL UNLABELED BRANCHES AS BACKGROUND
########################################################

import sys
import re

input_tree = sys.argv[1]
cleaned_tree = sys.argv[2]

with open(input_tree, 'r') as i_t:
    tree_string = i_t.read().rstrip()  # Combine reading and rstrip()

pattern = re.compile(r'([a-zA-Z0-9_.-]+):')

def replace_function(match):
	return match.group(1) + '{BACKGROUND}:'  # Append to captured group

tree_string = pattern.sub(replace_function, tree_string)

with open(cleaned_tree, 'w') as c_t:
    c_t.write(tree_string)
