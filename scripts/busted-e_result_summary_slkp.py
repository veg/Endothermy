"""
Combine analysis results

Author:
    Sergei L Kosakovsky Pond (spond@temple.edu)

Version:
    v0.0.1 (2021-01-17)


"""

import argparse
import csv
import random
import os
import json
import sys
import re
import math
import numpy
import gzip
from scipy.stats import chi2
from scipy.stats import mstats


import progressbar

bar = progressbar.ProgressBar(max_value=progressbar.UnknownLength)

from collections import Counter

random.seed ()

arguments = argparse.ArgumentParser(description='Combine alignments into a single file, adding a reference sequence as well')

arguments.add_argument('-i', '--input',     help = 'Directories ', required = True, type = str, nargs = '*')
arguments.add_argument('-s', '--sim',       help = 'Simulation mode ', action = "store_true")
arguments.add_argument('-z', '--gzip',      help = '.gzipped input', action = "store_true")
arguments.add_argument('-j', '--json',      help = 'JSON output', required = True, type = str)
arguments.add_argument('-p', '--prefix',    help = 'Summary CSV prefix', required = False, type = str, default = "File,Omega,Weight,Filtered")
arguments.add_argument('-o', '--output',    help = 'Summary CSV output', required = True, type = str)
arguments.add_argument('-F', '--file',      help = 'file parts to include in names', required = False, type = int, default = 1)

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
    try:
        return omegas[str (len (omegas) - 3)]
    except:
        return omegas[str (len (omegas) - 2)]
    
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


def newick_parser(nwk_str, bootstrap_values, track_tags):
    clade_stack = []
    automaton_state = 0
    current_node_name = ""
    current_node_attribute = ""
    current_node_annotation = ""
    quote_delimiter = None
    name_quotes = {
      "'": 1,
      '"': 1
    }
    
    def add_new_tree_level():
      new_level = {
        "name": None
      };
      the_parent = clade_stack[len(clade_stack) - 1]
      if (not "children" in the_parent):
        the_parent["children"] = [];
      
      clade_stack.append (new_level);
      the_parent["children"].append(clade_stack[len(clade_stack) - 1]);
      clade_stack[len(clade_stack)-1]["original_child_order"] = len(the_parent["children"])
    

    def finish_node_definition():
      nonlocal current_node_name
      nonlocal current_node_annotation
      nonlocal current_node_attribute
      
      this_node = clade_stack.pop()
      if (bootstrap_values and "children" in this_node):
        this_node["bootstrap_values"] = current_node_name
      else:
        this_node["name"] = current_node_name
      
      this_node["attribute"] = current_node_attribute
      this_node["annotation"] = current_node_annotation
      
      try:
      
          if not 'children' in this_node:
            node_tag = "background"
            for k, v in tags.items():
                if this_node["name"].find (k) >= 0:
                    node_tag = v
                    break
          else:
            '''
            counts = {}
            node_tag = ""
            for n in this_node['children']:
                counts[n["tag"]] = 1 + (counts[n["tag"]] if n["tag"] in counts  else 0)
            if len (counts) == 1:
                node_tag = list (counts.keys())[0]
            '''
            node_tag = "test"
        
          this_node["tag"] = node_tag
      except Exception as e:
        print ("Exception ", e)
        
      if track_tags is not None:
        track_tags[this_node["name"]] = [this_node["tag"], 'children' in this_node]
       
      current_node_name = ""
      current_node_attribute = ""
      current_node_annotation = ""
    

    def generate_error(location):
      return {
        'json': None,
        'error':
          "Unexpected '" +
          nwk_str[location] +
          "' in '" +
          nwk_str[location - 20 : location + 1] +
          "[ERROR HERE]" +
          nwk_str[location + 1 : location + 20] +
          "'"
      }


    tree_json = {
      "name" : "root"
    }
    
    clade_stack.append(tree_json);

    space = re.compile("\s")

    for char_index in range (len(nwk_str)):
      try:
        current_char = nwk_str[char_index]
        if automaton_state == 0:
           #look for the first opening parenthesis
           if (current_char == "("):
              add_new_tree_level()
              automaton_state = 1
        elif automaton_state == 1 or automaton_state == 3:
            #case 1: // name
            #case 3: { // branch length
            #reading name
            if (current_char == ":"):
              automaton_state = 3;
            elif current_char == "," or current_char == ")":
              try:
                finish_node_definition()
                automaton_state = 1
                if (current_char == ","):
                  add_new_tree_level()
              except Exception as e:
                return generate_error(char_index)
              
            elif (current_char == "("):
              if len(current_node_name) > 0:
                return generate_error(char_index);
              else:
                add_new_tree_level()
              
            elif (current_char in name_quotes):
              if automaton_state == 1 and len(current_node_name) == 0 and len (current_node_attribute) == 0 and len (current_node_annotation) == 0:
                automaton_state = 2
                quote_delimiter = current_char
                continue
              return generate_error(char_index)
            else:
              if (current_char == "{"):
                if len (current_node_annotation):
                  return generate_error(char_index)
                else:
                  automaton_state = 4
              else:
                if (automaton_state == 3):
                  current_node_attribute += current_char;
                else:
                  if (space.search(current_char)):
                    continue;
                  if (current_char == ";"):
                    char_index = len(nwk_str)
                    break
                  current_node_name += current_char;
        elif automaton_state == 2: 
            # inside a quoted expression
            if (current_char == quote_delimiter):
              if (char_index < len (nwk_str - 1)):
                if (nwk_str[char_index + 1] == quote_delimiter):
                  char_index+=1
                  current_node_name += quote_delimiter;
                  continue;

              quote_delimiter = 0
              automaton_state = 1
              continue
            else:
              current_node_name += current_char;
        elif automaton_state == 4:
           ##inside a comment / attribute
            if (current_char == "}"):
              automaton_state = 3
            else:
              if (current_char == "{"):
                return generate_error(char_index);
              current_node_annotation += current_char;
      except Exception as e:
        return generate_error(char_index);

    if (len (clade_stack) != 1):
      return generate_error(len (nwk_str) - 1);

    if (len (current_node_name)):
        tree_json['name'] = current_node_name;

    return {
      'json': tree_json,
      'error': None
    }
    

class result_reader (object):
    def __init__(self, file_name):
        self.file_name = file_name
     
    def __enter__(self):
        if settings.gzip:
            self.file = gzip.open(self.file_name)
        else:
            self.file = open(self.file_name, 'r')
        return self.file
 
    def __exit__(self, *args):
        self.file.close()
        
    
for basedir in settings.input:
     for root, dirs, files in os.walk(basedir):
        for each_file in sorted(files):
            #if i == 10: break
            file_name, file_ext =  os.path.splitext (each_file)
            #print (file_ext, file = sys.stderr)
            if settings.gzip:
                file_name,file_ext =  os.path.splitext (file_name)
            if file_ext == '.json':
                parts = file_name.split ('.')
                #print (parts)
                with result_reader (os.path.join (root, each_file)) as fh:
                    try:
                        i += 1
                        bar.update(i)

                        results = json.load (fh)
                        pv = results['test results']['p-value']
                        lrt = results['test results']['LRT']
                        logl = results["fits"]["Unconstrained model"]["Log Likelihood"]
                        aic = results["fits"]["Unconstrained model"]["AIC-c"]
                        omega3 = get_omega3
                        if settings.sim:
                            file_key = "|".join (parts[0:-1])
                        else:
                            file_key = ".".join (parts[0:settings.file])
                            
                        #print (file_name, file = sys.stderr)

                        if not file_key in by_file:
                            by_file[file_key] = {}
                        
                        model_key = parts[-1].split('-')
                        
                        
                        if len(model_key) > 1:
                            model_key = "BUSTED-E"
                        else:
                            model_key = "BUSTED"
                            
                        tree = newick_parser (results["input"]["trees"]["0"],{},{})["json"]
                        
                        descendants_by_name = {}
                        
                        def traverse_children (n):
                            
                            if "children" in n:
                                descendants_by_name [n["name"]] = set ()
                                for c in n["children"]:
                                    descendants_by_name [n["name"]].update (traverse_children (c))
                            else:
                                descendants_by_name [n["name"]] = set ([n["name"]])
                                
                            return descendants_by_name [n["name"]]
                        
                        
                       
                        traverse_children (tree)
                        
                        gaps_by_seq  = Counter ()
                        gaps_by_site = [] 
                            
                        for s,subs in results["substitutions"]["0"].items():
                            gaps = set ()
                            for n, state in subs.items():
                                if state == '---':
                                    gaps.update (descendants_by_name [n])
                            gaps_by_site.append (len (gaps))
                            for s in gaps:
                                gaps_by_seq [s] += 1


                        #print (gaps_by_seq)
                        #print (gaps_by_site)
                                                        
                        #sys.exit (0)
                       
                        by_file[file_key][model_key] = {'p' : pv, 
                                                        'LR' : lrt, 
                                                        'logL' : logl, 
                                                        'AIC' : aic, 
                                                        'omega'  : get_omega3 (results),
                                                        'omega1' : get_omega1 (results),
                                                        'omega2' : get_omega2 (results),
                                                        'SRV' : get_srv (results),
                                                        'runtime' : results["timers"]["Overall"]["timer"]
                                                        
                                                        }
                                                        
                        if model_key == "BUSTED-E":
                            by_file[file_key][model_key]['N'] = results['input']['number of sequences']
                            by_file[file_key][model_key]['S'] = results['input']['number of sites']
                            by_file[file_key][model_key]['T'] = sum([v['MG94xREV with separate rates for branch sets'] if 'MG94xREV with separate rates for branch sets' in v else 0 for k,v in results['branch attributes']['0'].items()])
                            
                            by_file[file_key][model_key]['omega_e'] = get_omega_e (results)
                            
                            by_file[file_key]['gaps'] = {'sequence' : gaps_by_seq,'sites' : gaps_by_site}
                                
                            
                        
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
                    
                    

settings.output_file = open (settings.output, "w")

output_writer = csv.writer (settings.output_file)

headers = ['File','N','S','T','BUSTED_AIC','BUSTED_E_AIC','p','p_e','logL','logL_e']
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
means = [0,0,0,0,0]
medians = [[],[],[],[]]
ma_medians = [[], []]

output_writer.writerow (headers)

for file, result in sorted (by_file.items()):
    #print (len (result))
    if len (result) == 3:  
        try:
            row = [file, str (result['BUSTED-E']['N']), str (result['BUSTED-E']['S']), "%.2f" % result['BUSTED-E']['T']]
            row.append ("%.5f" % result['BUSTED']["AIC"])
            row.append ("%.5f" % (result['BUSTED-E']["AIC"] - result['BUSTED']["AIC"]))
            row.append ("%.5f" % result['BUSTED']["p"])
            row.append ("%.5f" % result['BUSTED-E']["p"])
            
            for model in ['BUSTED','BUSTED-E']:
                row.append ("%.5f" % (result[model]["logL"]))
                
            lrt = 2*(result['BUSTED-E']["logL"]-result['BUSTED']["logL"])
            if chi2.cdf (lrt,2) >= 0.95:
                pv['-E is better'] += 1
            else:
                pv['-E is better'] += 0
            
        
            k = 0

            for model in ['BUSTED','BUSTED-E']:
                row.append ("%.5f" % (result[model]["omega1"]["omega"]))
                row.append ("%.5f" % (result[model]["omega1"]["proportion"]))
                row.append ("%.5f" % (result[model]["omega2"]["omega"]))
                row.append ("%.5f" % (result[model]["omega2"]["proportion"]))
                row.append ("%.5f" % (result[model]["omega"]["omega"]))
                row.append ("%.5f" % (result[model]["omega"]["proportion"]))
                means[k] += result[model]["omega"]["omega"]
                means[k+1] += result[model]["omega"]["proportion"]
                medians[k].append (result[model]["omega"]["omega"])
                medians[k+1].append (100.*result[model]["omega"]["proportion"])
                k += 2
                
            means[4] += result['BUSTED-E']["omega_e"]["proportion"]

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
               
            min_AIC = min (result['BUSTED']["AIC"], result['BUSTED-E']["AIC"])
            
            aw = {
                'BUSTED' : math.exp ((min_AIC-result['BUSTED']["AIC"])*0.5),
                'BUSTED-E' : math.exp ((min_AIC-result['BUSTED-E']["AIC"])*0.5)
            }
            
            
            
            aws = sum (aw.values())
            pv_MA = 0
            omega_p = 0
            weight_p = 0
            
            for i,model in enumerate (['BUSTED','BUSTED-E']):
               pv_MA += result [model]['p'] * aw[model] / aws
               omega_p  += medians[2*i][-1] * aw[model] / aws
               weight_p += medians[2*i+1][-1] * aw[model] / aws
               
            ma_medians[0].append (omega_p)
            ma_medians[1].append (weight_p)
            
               
            pv['MA'] += 1 if pv_MA <= 0.05 else 0
               
            pv['N'] += 1
            #print (result['MH']['p'], file = sys.stderr)
            output_writer.writerow (row)
        except Exception as e:
            print ("Processing", e, file, file = sys.stderr)
            pass
            
            
with open (settings.json, "w") as fh:
    json.dump (by_file, fh)                   
    
for k,n in pv.items():
    if k != 'N':
        print ("%s = %.2g" % (k, n / pv['N']), file = sys.stderr)
    else:
        print ("N = %d" % pv['N'], file = sys.stderr)

mean_labels = ['BUSTED omega','BUSTED fraction','BUSTED-E omega', 'BUSTED-E fraction', 'error fraction']
for i,k in enumerate(means):
    print (mean_labels[i], file = sys.stderr)
    if i % 2 == 0 and i < 4:
        print ("\t%.2g" % (k / pv['N']), file = sys.stderr)
        if i < 4:
            print ("\tmedian = %.2g" % tuple(mstats.mquantiles(medians[i],[0.5])), file = sys.stderr)
    else:
        print ("\t%.2g" % (k / pv['N']*100.), file = sys.stderr)
        if i < 4:
            print ("\tmedian = %.2g" % tuple(mstats.mquantiles(medians[i],[0.5])), file = sys.stderr)

busted_medians = mstats.mquantiles(medians[0],[0.25,0.75])
busted_frac_medians = mstats.mquantiles(medians[1],[0.25,0.75])
bustede_medians = mstats.mquantiles(medians[2],[0.25,0.75])
bustede_frac_medians = mstats.mquantiles(medians[3],[0.25,0.75])
bustedma_medians = mstats.mquantiles(ma_medians[0],[0.25,0.75])
bustedma_frac_medians = mstats.mquantiles(ma_medians[1],[0.25,0.75])
    
prfx = "\t".join (settings.prefix.split (","))
    
print (
    "%s\tBUSTED\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g" % (prfx, pv["BUSTED"]/pv["N"], 0.,              numpy.mean (medians[0]), busted_medians[0], busted_medians[1], numpy.mean (medians[1]), busted_frac_medians[0], busted_frac_medians[1]),
    "\n%s\tBUSTED-E\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g" % (prfx, pv["BUSTED-E"]/pv["N"], means[4], numpy.mean (medians[2]), bustede_medians[0], bustede_medians[1], numpy.mean (medians[3]), bustede_frac_medians[0], bustede_frac_medians[1]),
    "\n%s\tMA\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g" % (prfx, pv["MA"]/pv["N"], 0.,                    numpy.mean (ma_medians[0]), bustedma_medians[0], bustedma_medians[1], numpy.mean (ma_medians[1]), bustedma_frac_medians[0], bustedma_frac_medians[1]),
    file = sys.stdout
)

                        
