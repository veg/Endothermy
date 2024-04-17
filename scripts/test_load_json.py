"""
Combine analysis results for filtered and unfiltered aBSREL, MEME, and PAML emulator

Author:
    Avery Selberg (avery@temple.edu)

Version:
    v0.0.1 (2023-11-27)


"""

import argparse
import random
import os
import json
import sys
import math


arguments = argparse.ArgumentParser(description='Process BUSTED-PH results')

arguments.add_argument('-i', '--input',  help = 'Directories ', required = True, type = str)
arguments.add_argument('-o', '--output', help = 'Output file full path ', required = True, type = str) 

settings = arguments.parse_args()

by_file = {}
class CustomDecoder(json.JSONDecoder):
    def decode(self, s, **kwargs):
        try:
            # Try decoding the JSON string using the default decoder
            return super().decode(s, **kwargs)
        except json.JSONDecodeError:
            # If decoding fails, handle non-standard values
            s = s.replace('inf', 'Infinity')
            s = s.replace('-inf', '-Infinity')
            return json.loads(s)

file_name, file_ext =  os.path.splitext (settings.input)
if file_ext == '.json':
	#denote key identifier (usually gene name)
	key = file_name.split ('.')[0]
	with open (settings.input, "r") as fh:
		print('settings.input: ', settings.input, '\n\n')
		try:
			res = json.load (fh, cls=CustomDecoder)  
			print(res)
			exit()
			if not key in by_file:
				by_file [key] = {}
			if 'BUSTED-PH.json' in each_file:
				by_file[key]['BUSTED-PH_' + filt] = {}
				by_file[key]['BUSTED-PH_' + filt]["fits"] = res["fits"]
		except json.JSONDecodeError as e:
			print(f"Error decoding JSON in {os.path.join(file_name + file_ext)}: {e}", file=sys.stderr)
		except Exception as e:
			print (e, file = sys.stderr)
		#	exit()
			pass

print(settings.output)
with open(settings.output, 'w') as f:
	json.dump (by_file, f)


        
