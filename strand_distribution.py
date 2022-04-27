#!usr/bin/env python3

#this script is to be ran from the results folder of Ribosemap run to output unique read percentage
#Activate the conda enviroment for Ribosemap before running this script

import subprocess
import pandas as pd


import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-p","--files",type=argparse.FileType('r'),nargs='+', help="File(s) to be analyzed for strand info")
parser.add_argument("-o","--outprefix", help="Prefix of output file")
parser.add_argument("-s","--orgsep", help="separate the nucleus and mitochondria",action='store_true')
args= parser.parse_args()
if args.orgsep:
    df = pd.DataFrame(columns=['File_name', 'nucl_pos','nucl_neg','chrM_pos','chrM_neg'])
else:
    df = pd.DataFrame(columns=['File_name','pos','neg'])

#subprocess.call("conda activate conda_env", shell=True, universal_newlines=True)
for file in args.files:
    print("Analyzing "+file.name)
    if args.orgsep:
        nucl_pos=subprocess.check_output("grep -v chrM "+file.name+"| grep '\s+' | wc -l", shell=True, universal_newlines=True)
        nucl_neg=subprocess.check_output("grep -v chrM "+file.name+"| grep '\s-' | wc -l", shell=True, universal_newlines=True)
        chrM_pos=subprocess.check_output("grep chrM "+file.name+"| grep '\s+' | wc -l", shell=True, universal_newlines=True)
        chrM_neg=subprocess.check_output("grep chrM "+file.name+"| grep '\s-' | wc -l", shell=True, universal_newlines=True)
        df_row={'File_name':file.name, 'nucl_pos':nucl_pos.strip('\n'),'nucl_neg':nucl_neg.strip('\n'),'chrM_pos':chrM_pos.strip('\n'),'chrM_neg':chrM_neg.strip('\n')}
    else:
        pos=subprocess.check_output("cat "+file.name+"| grep '\s+' | wc -l", shell=True, universal_newlines=True)
        neg=subprocess.check_output("cat "+file.name+"| grep '\s-' | wc -l", shell=True, universal_newlines=True)
        df_row={'File_name':file.name, 'pos':pos.strip('\n'),'neg':neg.strip('\n')}

    
    new_row=pd.Series(df_row)
    df=df.append(new_row, ignore_index=True)

df.to_csv('{}_strandcounts.tsv'.format(args.outprefix), sep='\t', index=False)
print("Stats saved in {}_strandcounts.tsv".format(args.outprefix))