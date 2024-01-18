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

#import progressbar

#bar = progressbar.ProgressBar(max_value=progressbar.UnknownLength)


cutoff = 100.

arguments = argparse.ArgumentParser(description='Process BUSTED-PH results')

arguments.add_argument('-i', '--input',  help = 'Directories ', required = True, type = str, action = 'append')
arguments.add_argument('-o', '--output', help = 'Output file full path ', required = True, type = str) 

settings = arguments.parse_args()

by_file = {}

for basedir in settings.input:
	for root, dirs, files in os.walk(basedir):
		#print(files)
		for each_file in files:
			file_name, file_ext =  os.path.splitext (each_file)
			if file_ext == '.json':
				if os.path.getsize(os.path.join (root, each_file)) == 0:
					continue
				key = file_name.split ('.')[0]
				with open (os.path.join (root, each_file), "r") as fh:
					try:
						res = json.load (fh)  
						if not key in by_file:
							by_file [key] = {}
						if "filter" in each_file:
							filt = "FILTERED"
						else:
							filt = "UNFILTERED"
						if 'BUSTED.json' in each_file:
							by_file[key]['BUSTED'] = {}
							by_file[key]['BUSTED']['fits'] = res['fits']
						if 'BUSTED-E.json' in each_file:
							by_file[key]['BUSTED_E'] = {}
							by_file[key]['BUSTED_E']['fits'] = res['fits']
						if 'BUSTED-PH.json' in each_file:
							by_file[key]['BUSTED-PH_' + filt] = {}
							by_file[key]['BUSTED-PH_' + filt]["fits"] = res["fits"]
							by_file[key]['BUSTED-PH_' + filt]["tree"] = res["input"]["trees"]["0"]
							by_file[key]['BUSTED-PH_' + filt]["sequences"] = res["input"]["number of sequences"]
							by_file[key]['BUSTED-PH_' + filt]["sites"] = res["input"]["number of sites"]
							by_file[key]['BUSTED-PH_' + filt]["tests"] = {
							    "FG" : res["test results"],
							    "BG" : res["test results background"],
							    "DIFF" : res["test results shared distributions"]
							}
							by_file [key]['BUSTED-PH_' + filt]["tested"] = res["tested"]
							by_file [key]['BUSTED-PH_' + filt]["distributions"] = {
							    'unconstrained' : res['fits']['Unconstrained model'],
							    'shared' : res['fits']['Shared distribution model']
							}
							by_file [key]['BUSTED-PH_' + filt]["lengths"] = dict ([(b,v['unconstrained']) for b,v in res["branch attributes"]["0"].items()])
							if 'constrained' in res['Evidence Ratios']:
							    by_file [key]['BUSTED-PH_' + filt]["fg_selected"] = len([k for k in res['Evidence Ratios']['constrained'][0] if k >= cutoff])
							else:
							    by_file [key]['BUSTED-PH_' + filt]["fg_selected"] = 0
							if 'constrained background' in res['Evidence Ratios']:
								by_file [key]['BUSTED-PH_' + filt]["bg_selected"] = len([k for k in res['Evidence Ratios']['constrained background'][0] if k >= cutoff])
							else:
							    by_file [key]["bg_selected"]  = 0
					except Exception as e:
						print (e, file = sys.stderr)
						print(os.path.join(basedir, each_file))
					#	exit()
						pass

print(settings.output)
with open(settings.output, 'w') as f:
	json.dump (by_file, f)


        
