"""
Combine analysis results

Author:
    Sergei L Kosakovsky Pond (spond@temple.edu)

Version:
    v0.0.1 (2021-01-17)


"""
### TO USE: 
### python result_summary.py -i /path/to/results/ > results.csv

import argparse
import csv
import random
import os
import json
import sys
import re
import math
import numpy

import progressbar

bar = progressbar.ProgressBar(max_value=progressbar.UnknownLength)

from collections import Counter

random.seed ()

arguments = argparse.ArgumentParser(description='Combine alignments into a single file, adding a reference sequence as well')

arguments.add_argument('-i', '--input',  help = 'Directories ', required = True, type = str, nargs = '*')

settings = arguments.parse_args()

by_file = {}

timer = 0
count = 0
tags = {}

def get_omega3 (fit):
    omegas = fit["fits"]["Unconstrained model"]["Rate Distributions"]["Test"]
    return omegas[str (len (omegas) - 1)]

def get_omega2 (fit):
    omegas = fit["fits"]["Unconstrained model"]["Rate Distributions"]["Test"]
    return omegas[str (len (omegas) - 2)]

def get_omega1 (fit):
    omegas = fit["fits"]["Unconstrained model"]["Rate Distributions"]["Test"]
    return omegas[str (len (omegas) - 3)]
    
def get_omega_e (fit):
    omegas = fit["fits"]["Unconstrained model"]["Rate Distributions"]["Test"]
    return omegas["0"]

def get_srv (fit):
    try:
        srv = fit["fits"]["Unconstrained model"]["Rate Distributions"]["Synonymous site-to-site rates"]
        s  = 0
        s2 = 0
        for k, v in srv.items():
            s += v["rate"] * v["proportion"] 
            s2 += v["rate"]*v["rate"] * v["proportion"] 
        return math.sqrt ((s2-s*s)/s)
    except:
        return 0

 
i = 0

    
for basedir in settings.input:
     for root, dirs, files in os.walk(basedir):
        for each_file in sorted(files):
            #if i == 10: break
            file_name, file_ext =  os.path.splitext (each_file)
            #print (file_ext, file = sys.stderr)
            if file_ext == '.json':
                parts = file_name.split ('.')
                #print (parts)
                with open (os.path.join (root, each_file), "r") as fh:
                    try:
                        i += 1
                        bar.update(i)

                        results = json.load (fh)
                        pv = results['test results']['p-value']
                        lrt = results['test results']['LRT']
                        logl = results["fits"]["Unconstrained model"]["Log Likelihood"]
                        aic = results["fits"]["Unconstrained model"]["AIC-c"]
                        omega3 = get_omega3
                        file_key = parts[0]

                        if not file_key in by_file:
                            by_file[file_key] = {}
                        
                        model_key = parts[-1].split('-')
                        
                        
                        if len(model_key) > 1:
                            model_key = "BUSTED-E"
                        else:
                            model_key = "BUSTED"
                            
                       
                        by_file[file_key][model_key] = {'p' : pv, 
                                                        'LR' : lrt, 
                                                        'logL' : logl, 
                                                        'AIC' : aic, 
                                                        'omega'  : get_omega3 (results),
                                                        'omega1' : get_omega1 (results),
                                                        'omega2' : get_omega2 (results),
                                                        'SRV' : get_srv (results),
                                                        'runtime' : results["timers"]["Overall"]["timer"]}
                                                        
                        if model_key == "BUSTED-E":
                            by_file[file_key][model_key]['N'] = results['input']['number of sequences']
                            by_file[file_key][model_key]['S'] = results['input']['number of sites']
                            by_file[file_key][model_key]['T'] = sum([v['unconstrained'] if 'unconstrained' in v else 0 for k,v in results['branch attributes']['0'].items()])
                            
                            by_file[file_key][model_key]['omega_e'] = get_omega_e (results)
                                
                            
                        
                        prior = by_file[file_key][model_key]['omega']
                        if prior ['proportion'] == 1: 
                                prior = 0
                        else:
                            prior = prior ['proportion'] / (1 -  prior ['proportion'])
                        EB = 0
                        if prior > 0:   
                            for b, data in results['branch attributes']['0'].items():
                                if "Posterior prob omega class by site" in data:
                                    for s in data["Posterior prob omega class by site"][-1]:
                                        if s > 0 and s < 1 and s/(1-s) / prior >= 100:
                                            EB+=1
                                            
                        by_file[file_key][model_key]['EB'] = EB
                        by_file[file_key][model_key]['LRF'] = None
                        if "Site Log Likelihood" in results:
                            if "unconstrained" in results["Site Log Likelihood"] and "optimized null" in results["Site Log Likelihood"]:
                                lldiff = [ results["Site Log Likelihood"]["unconstrained"][0][k] - i for (k,i) in enumerate (results["Site Log Likelihood"]["optimized null"][0])]
                                S = sum(lldiff)
                                if S > 2:
                                    SD = S * 0.8
                                    for idx, v in enumerate (sorted (lldiff, key = lambda x: -x)):
                                        SD -= v
                                        if SD <= 0:
                                            break
                                    by_file[file_key][model_key]['LRF'] = idx + 1

                            
                    except Exception as e:
                        print (e, each_file, file = sys.stderr)
                        pass
                    
                    


output_writer = csv.writer (sys.stdout)
headers = ['File','N','S','T','BUSTED_AIC','BUSTED_E_AIC','p','p_e']
for model in ['BUSTED','BUSTED-E']:
    headers.append ("omega1_%s" % model)
    headers.append ("p1_%s" % model)
    headers.append ("omega2_%s" % model)
    headers.append ("p2_%s" % model)
    headers.append ("omega3_%s" % model)
    headers.append ("p3_%s" % model)
    
headers.append ("omega_e")
headers.append ("prop_e")

for model in ['BUSTED','BUSTED-E']:
    headers.append ("Runtime_%s" % model)

for model in ['BUSTED','BUSTED-E']:
    headers.append ("EB_%s" % model)

for model in ['BUSTED','BUSTED-E']:
    headers.append ("LRF_%s" % model)

for model in ['BUSTED','BUSTED-E']:
     headers.append ("SRV_%s" % model)

pv = Counter ()

output_writer.writerow (headers)

for file, result in sorted (by_file.items()):
    #print (len (result))
    if len (result) == 2:  
        try:
            row = [file, str (result['BUSTED-E']['N']), str (result['BUSTED-E']['S']), "%.2f" % result['BUSTED-E']['T']]
            row.append ("%.5f" % result['BUSTED']["AIC"])
            row.append ("%.5f" % (result['BUSTED-E']["AIC"] - result['BUSTED']["AIC"]))
            row.append ("%.5f" % result['BUSTED']["p"])
            row.append ("%.5f" % result['BUSTED-E']["p"])
            for model in ['BUSTED','BUSTED-E']:
                row.append ("%.5f" % (result[model]["omega1"]["omega"]))
                row.append ("%.5f" % (result[model]["omega1"]["proportion"]))
                row.append ("%.5f" % (result[model]["omega2"]["omega"]))
                row.append ("%.5f" % (result[model]["omega2"]["proportion"]))
                row.append ("%.5f" % (result[model]["omega"]["omega"]))
                row.append ("%.5f" % (result[model]["omega"]["proportion"]))

            row.append ("%.5f" % (result['BUSTED-E']["omega_e"]["omega"]))
            row.append ("%.5f" % (result['BUSTED-E']["omega_e"]["proportion"]))


            for model in ['BUSTED','BUSTED-E']:
                row.append ("%.5f" % (result[model]["runtime"]))

            for model in ['BUSTED','BUSTED-E']:
                row.append ("%d" % (result[model]["EB"]))

            for model in ['BUSTED','BUSTED-E']:
                if result[model]["LRF"] == None:
                   row.append ("N/A")
                else:
                    row.append ("%d" % (result[model]["LRF"]))

            for model in ['BUSTED','BUSTED-E']:
               row.append ("%.5f" % (result[model]["SRV"]))
            
            for model in ['BUSTED','BUSTED-E']:
               pv[model] += 1 if result [model]['p'] <= 0.05 else 0
            pv['N'] += 1
            #print (result['MH']['p'], file = sys.stderr)
            output_writer.writerow (row)
        except Exception as e:
            print (e, file, file = sys.stderr)
            pass
                    
    
print (pv, file = sys.stderr)

                        
