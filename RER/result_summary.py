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
import gzip

import progressbar

bar = progressbar.ProgressBar(max_value=progressbar.UnknownLength)

from collections import Counter

random.seed ()

arguments = argparse.ArgumentParser(description='Combine alignments into a single file, adding a reference sequence as well')

arguments.add_argument('-o', '--output',     help = 'Output ', required = True, type = str)
arguments.add_argument('-i', '--input',     help = 'Directories ', required = True, type = str, nargs = '*')
arguments.add_argument('-z', '--gzip',      help = '.gzipped input', action = "store_true")

settings = arguments.parse_args()

by_file = {}

timer = 0
count = 0
tags = {}

i = 0

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
   
output_file = csv.writer (open (settings.output, "w"))

output_file.writerow (['gene','fg','bg','ratio','p','branch','alternative','null']) 
    
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
                
                     
                
                with result_reader (os.path.join (root, each_file)) as fh:
                    try:
                        i += 1
                        bar.update(i)

                        results = json.load (fh)
                        rates = results['fits']['Proportional Partitioned']['Rate Distributions']
                        rfg = rates['branch length scaler for test']
                        rbg = rates['branch length scaler for background']
                        pv = results['test results']['Proportional:Proportional Partitioned']['Uncorrected P-value']
                        branch =','.join(map(str, [d for d in results['branch level analysis']]))
                        alternative = ','.join(map(str, [results['branch level analysis'][d]['alternative']['p-value']for d in results['branch level analysis']]))
                        null = ','.join(map(str, [results['branch level analysis'][d]['null']['p-value'] for d in results['branch level analysis']]))
                        row = [str(k) for k in [parts[0],rfg,rbg,rfg/rbg,pv,branch, alternative, null]]
                        output_file.writerow (row)
                        

                            
                    except Exception as e:
                        print (e, each_file, file = sys.stderr)
                        pass
                    
                    


                        
