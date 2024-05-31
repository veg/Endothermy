import argparse
import random
import os
import json
import sys
import math
import statsmodels.api as sm
import scipy

cutoff = 100.

arguments = argparse.ArgumentParser(description='Process BUSTED-PH results')

arguments.add_argument('-i', '--input', help='Directories', required=True, type=str, action='append')
arguments.add_argument('-o', '--output', help='Output file full path', required=True, type=str) 

settings = arguments.parse_args()

by_file = {}

class CustomDecoder(json.JSONDecoder):
    def decode(self, s, **kwargs):
        try:
            return super().decode(s, **kwargs)
        except json.JSONDecodeError:
            s = s.replace('inf', '"Infinity"')
            s = s.replace('-inf', '"-Infinity"')
            return json.loads(s)

file_long = '.fasta.RD.SA.codons.cln.trim67.fa.'
busted_e_fdr = {}

for basedir in settings.input:
	#for basedir, dirs, all_files in os.walk(basedir):
		#files_endo = [x for x in all_files if x.endswith('BUSTED.json')]
		#files = [x for x in files_endo if 'Endothermy' not in x]
	files = [f for f in os.listdir(basedir) if os.path.isfile(os.path.join(basedir, f)) and f.endswith('BUSTED.json')]
	##### REMOVE LATER #####
	#files = ['E2148_NT_Cleaned_Aligned.fasta.RD.SA.codons.cln.trim67.fa.BUSTED.json']
	##### REMOVE LATER #####

	for each_file in files:
		file_name, file_ext = os.path.splitext(each_file)
		##### Get BUSTED file results ##### 
		if 'BUSTED.json' in each_file:
			if os.path.getsize(os.path.join(basedir, each_file)) == 0:
				continue
			key = file_name.split('.')[0]
			p_val_e = 0
			with open(os.path.join(basedir, each_file), "r") as fh:
				try:
					res = json.load(fh, cls=CustomDecoder)
					if key not in by_file:
						by_file[key] = {}
					by_file[key]['BUSTED'] = {}
					by_file[key]['BUSTED']['fits'] = res['fits']
					by_file[key]['BUSTED']["test results"] = res['test results']
				except Exception as e:
					print(e, file=sys.stderr)
					print(os.path.join(basedir, each_file))
			##### Get BUSTED-E file results #####
			busted_e_file_path = os.path.join(basedir, key + file_long + "BUSTED-E.json")
			#print('busted-e path: ', busted_e_file_path, '\n')
			if 'Endothermy' in busted_e_file_path:
				print('ENDOTHERMY busted-e', busted_e_file_path, '\n')
				continue
			if os.path.isfile(busted_e_file_path):
				with open(busted_e_file_path, "r") as fh:
					try:
						res = json.load(fh, cls=CustomDecoder)
						by_file[key]['BUSTED_E'] = {}
						by_file[key]['BUSTED_E']['fits'] = res['fits']
						by_file[key]['BUSTED_E']["test results"] = res['test results']
						b_e_logl = by_file[key]['BUSTED_E']['fits']['Unconstrained model']['Log Likelihood']
						b_logl = by_file[key]['BUSTED']['fits']['Unconstrained model']['Log Likelihood']
						##### Perform LRT on BUSTED-E vs BUSTED #####
						LR = -2 * (b_logl - b_e_logl)    # alex check
						p_val_e = scipy.stats.chi2.sf(LR, 2) / 2    # alex check
						busted_e_fdr[key] = p_val_e    # alex check
						by_file[key]['BUSTED_E']["test for error"] = {"LRT": LR, "p-value": p_val_e}    # alex check
					except Exception as e:
						print(e, file=sys.stderr)
						print('BUSTED-E error for: ', busted_e_file_path)
			
			##### Choose which BUSTED-PH file to get, based on BUSTED-E/BUSTED LRT #####
			if p_val_e < 0.05:    # alex check
				filt = '.filtered'    # alex check
			if p_val_e >= 0.05:    # alex check
				filt = ''    # alex check
			for i in range(1,5):
				scenario = str(i)
				busted_ph_file_path = os.path.join(basedir, key + file_long + "Scenario_" + scenario + filt + '.BUSTED-PH.json')    # alex check
				#print('busted-ph no fdr path: ', busted_ph_file_path, '\n')
				if 'Endothermy' in busted_ph_file_path:
					print('ENDOTHERMY busted-ph no fdr', busted_ph_file_path, '\n')
					continue
				##### Get BUSTED-PH results based on BUSTED-E/BUSTED LRT #####
				if os.path.isfile(busted_ph_file_path):
					with open(busted_ph_file_path, 'r') as fh:
						try:
							res = json.load(fh, cls=CustomDecoder)
							by_file[key]['BUSTED-PH_' + scenario] = {}
							by_file[key]['BUSTED-PH_' + scenario]['filtered?'] = "unfiltered"  if filt == "" else "filtered"    # alex check
							by_file[key]['BUSTED-PH_' + scenario]["fits"] = res["fits"]
							by_file[key]['BUSTED-PH_' + scenario]["tree"] = res["input"]["trees"]["0"]
							by_file[key]['BUSTED-PH_' + scenario]["sequences"] = res["input"]["number of sequences"]
							by_file[key]['BUSTED-PH_' + scenario]["sites"] = res["input"]["number of sites"]
							by_file[key]['BUSTED-PH_' + scenario]["tests"] = {
								"FG" : res["test results"],
								"BG" : res["test results background"],
								"DIFF" : res["test results shared distributions"]
							}
							by_file [key]['BUSTED-PH_' + scenario]["tested"] = res["tested"]
							by_file [key]['BUSTED-PH_' + scenario]["distributions"] = {
								'unconstrained' : res['fits']['Unconstrained model'],
								'shared' : res['fits']['Shared distribution model'] 
							}
							by_file [key]['BUSTED-PH_' + scenario]["lengths"] = dict ([(b,v['unconstrained']) for b,v in res["branch attributes"]["0"].items()])
						except Exception as e:
							print(e, file =sys.stderr)
							print('BUSTED-PH filtered error for: ', busted_ph_file_path)
	
##### Do FDR calculation on	BUSTED-E vs BUSTED LRT p-values #####	
sorted_befdr = dict(sorted(busted_e_fdr.items(), key=lambda item: item[1]))    # alex check
fdr_calc = {    # alex check
	key: {    # alex check
	        'uncorrected': value,    # alex check
			'uncorrected significance' : value < 0.05,    # alex check
			'corrected': (rank + 1) / len(sorted_befdr) * 0.05,    # alex check
			'fdr significance': value < (rank + 1) / len(sorted_befdr) * 0.05,    # alex check
			'new file needed' : (value < (rank + 1) / len(sorted_befdr) * 0.05) != (value < 0.05)    # alex check
	}    # alex check
	for rank, (key, value) in enumerate(sorted(sorted_befdr.items(), key=lambda x: x[1]))    # alex check
}    # alex check
#print('fdr_calc dict', fdr_calc, '\n')

##### Choose which BUSTED-PH results file (filtered/unfiltered) to get based on FDR #####
##### Avery's thought process: check if FDR leads to the same results or different ('new file needed')
##### (was significant for error but now is not) and get other busted-ph files if necessary
busted_ph_fdr = {'Foreground' : {}, 'Background' : {}, 'Difference' : {}}
scenarios = ['1', '2', '3', '4']
bph_3_tests = ['Foreground', 'Background', 'Difference']
busted_ph_fdr = {sce: {test: {} for test in bph_3_tests} for sce in scenarios}
for key, value in fdr_calc.items():
	for i in range(1, 5):
		scenario = str(i)
		# For cases where busted-ph scenario does not exist but BUSTED-E files does:
		if 'BUSTED-PH_' + scenario not in by_file[key].keys():
			continue

		if value['new file needed'] == False:
			by_file[key]['BUSTED-PH_' + scenario + '_fdr'] = by_file[key]['BUSTED-PH_' + scenario]
			busted_ph_fdr[scenario]['Foreground'][key] = by_file[key]['BUSTED-PH_' + scenario + '_fdr']["tests"]["FG"]['p-value']
			busted_ph_fdr[scenario]['Background'][key] = by_file[key]['BUSTED-PH_' + scenario + '_fdr']["tests"]["BG"]['p-value']
			busted_ph_fdr[scenario]['Difference'][key] = by_file[key]['BUSTED-PH_' + scenario + '_fdr']["tests"]["DIFF"]['p-value']
		elif value['new file needed'] == True:
			filt = 'AVERY'
			if by_file[key]['BUSTED-PH_' + scenario]['filtered?'] == 'unfiltered':
				filt = '.filtered'
			elif by_file[key]['BUSTED-PH_' + scenario]['filtered?'] == 'filtered':
				filt = ''
			##### Get BUSTED-PH results (FDR filtered version) #####
			busted_ph_file_path = os.path.join(basedir, key + file_long + "Scenario_" + scenario + filt + '.BUSTED-PH.json')

			if 'Endothermy' in busted_ph_file_path:
				print('ENDOTHERMY BUSTED-PH FDR', busted_ph_file_path, '\n')
				continue
			if os.path.isfile(busted_ph_file_path):
				with open(busted_ph_file_path, 'r') as fh:
					try:
						res = json.load(fh, cls=CustomDecoder)
						by_file[key]['BUSTED-PH_' + scenario + '_fdr'] = {}
						by_file[key]['BUSTED-PH_' + scenario + '_fdr']['filtered?'] = "unfiltered" if filt == "" else "filtered"
						by_file[key]['BUSTED-PH_' + scenario + '_fdr']["fits"] = res["fits"]
						by_file[key]['BUSTED-PH_' + scenario + '_fdr']["tree"] = res["input"]["trees"]["0"]
						by_file[key]['BUSTED-PH_' + scenario + '_fdr']["sequences"] = res["input"]["number of sequences"]
						by_file[key]['BUSTED-PH_' + scenario + '_fdr']["sites"] = res["input"]["number of sites"]
						by_file[key]['BUSTED-PH_' + scenario + '_fdr']["tests"] = {
							"FG": res["test results"],
							"BG": res["test results background"],
							"DIFF": res["test results shared distributions"]
						}
						by_file[key]['BUSTED-PH_' + scenario + '_fdr']["tested"] = res["tested"]
						by_file[key]['BUSTED-PH_' + scenario + '_fdr']["distributions"] = {
							'unconstrained': res['fits']['Unconstrained model'],
							'shared': res['fits']['Shared distribution model']
						}
						by_file[key]['BUSTED-PH_' + scenario + '_fdr']["lengths"] = dict((b, v['unconstrained']) for b, v in res["branch attributes"]["0"].items())
						busted_ph_fdr[scenario]['Foreground'][key] = by_file[key]['BUSTED-PH_' + scenario + '_fdr']["tests"]["FG"]['p-value']
						busted_ph_fdr[scenario]['Background'][key] = by_file[key]['BUSTED-PH_' + scenario + '_fdr']["tests"]["BG"]['p-value']
						busted_ph_fdr[scenario]['Difference'][key] = by_file[key]['BUSTED-PH_' + scenario + '_fdr']["tests"]["DIFF"]['p-value']
					except Exception as e:
						print(e, file=sys.stderr)
						print('BUSTED-PH filtered error for: ', busted_ph_file_path)


fdr_bph ={sce: {test: {} for test in bph_3_tests} for sce in scenarios}
#print(fdr_bph)
for scenario in scenarios:
	for test in bph_3_tests:
		for k, v in busted_ph_fdr[scenario].items():
			sorted_bphfdr = dict(sorted(busted_ph_fdr[scenario][test].items(), key=lambda item: item[1]))    # alex check
			fdr_bph[scenario][test] = { 
						key: {    # alex check
						        'uncorrected': value,    # alex check
								'uncorrected significance' : value < 0.05,    # alex check
								'corrected': (rank + 1) / len(sorted_bphfdr) * 0.05,    # alex check
								'fdr significance': value < (rank + 1) / len(sorted_bphfdr) * 0.05,    # alex check
						}    # alex check
						for rank, (key, value) in enumerate(sorted(sorted_bphfdr.items(), key=lambda x: x[1]))    
					}


# Add BUSTED-PH FDR correction to full dictionary/json
for key, value in by_file.items():
	for scenario in scenarios:
		busted_ph_key = 'BUSTED-PH_' + scenario + '_fdr'
		if busted_ph_key in by_file[key]:
			by_file[key]['BUSTED-PH_' + scenario + '_fdr']['fdr significance'] = {}
			by_file[key]['BUSTED-PH_' + scenario + '_fdr']['fdr significance'] = {
					'FG': fdr_bph[scenario]['Foreground'][key]['fdr significance'],
					'BG': fdr_bph[scenario]['Background'][key]['fdr significance'],
					'DIFF':fdr_bph[scenario]['Difference'][key]['fdr significance'],
				}

with open(settings.output, 'w') as f:
    json.dump(by_file, f)
