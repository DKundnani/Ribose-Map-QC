#!usr/bin/env python3
#Author: Deepali L. Kundnani
#Institution: Georgia Institute of Technology (Georgia Tech)
#This script is to be ran from the results folder of Ribosemap run to output unique read percentage
#Activate the conda enviroment for Ribosemap before running this script(or any that provides samtools and bedtools)

import subprocess
import pandas as pd
import os
import glob
import argparse
parser = argparse.ArgumentParser(description='Quality statitics of ribosemap runs for post Alignment, Coordinate, Composition modules')

parser.add_argument("-f", "--file", required=True, default='/storage/home/hcoda1/5/dkundnani3/p-fstorici3-0/rich_project_bio-storici/AGS/files', help="tab separated file with first column with library id mentioned in Ribosemap, and additional sample names if to be attached int the final file.")
parser.add_argument("-r", "--referencefasta" ,required=True, default='/storage/home/hcoda1/5/dkundnani3/p-fstorici3-0/rich_project_bio-storici/reference/sacCer2/sacCer2.fa', help="Specify the fasta file for reference of the sequenced Libraries. Make sure you have Ribosemap environment activated or have bedtools available to be used") 
args= parser.parse_args()

sample_names=pd.read_csv(args.file,sep='\t', header=None)

print("Sample information provided")
print(sample_names)

df = pd.DataFrame(columns=['Library','Counts_All','Counts_Nuclear','Counts_Mito','Comp_nucl_rA','Comp_nucl_rC','Comp_nucl_rG','Comp_nucl_rU','Freq_nucl_rA','Freq_nucl_rC','Freq_nucl_rG','Freq_nucl_rU','Comp_mito_rA','Comp_mito_rC','Comp_mito_rG','Comp_mito_rU','Freq_mito_rA','Freq_mito_rC','Freq_mito_rG','Freq_mito_rU','pos>1_nucleus/totpos', 'pos>1_mito/totpos'], dtype=object)

#subprocess.call("conda activate conda_env", shell=True, universal_newlines=True)
for lib in sample_names[0]:
    file=glob.glob("*{}*".format(lib))[0]
    print("Analyzing "+file)
    #Parameter SetI:Processing Statitics: Raw reads from fastq files, Aignment rate, percent reads barcoded, unique read percent
    
    #Paramter Set II: Counts, Composition and Frequency in Nucleus and Mitochondria
    '''
    chrX_reads=subprocess.check_output("grep chrX "+file+"/coordinate30/"+file+".bed | wc -l", shell=True, universal_newlines=True)
    #parameter 9: Reads in Chr Y
    chrY_reads=subprocess.check_output("grep chrY "+file+"/coordinate30/"+file+".bed | wc -l", shell=True, universal_newlines=True)
    
    
    #Paramter set 10:  Perc of rN(A,C,G,U) in mitochondria and nucleus.
    subprocess.call("bedtools getfasta -s -bedOut -fi "+args.referencefasta+" -bed "+file+"/coordinate30/"+file+"-chrM.bed > "+file+"/coordinate30/"+file+"-chrM_seq.bed",shell=True, universal_newlines=True)
    subprocess.call("bedtools getfasta -s -bedOut -fi "+args.referencefasta+" -bed "+file+"/coordinate30/"+file+"-chromosomes.bed > "+file+"/coordinate30/"+file+"-chromosomes_seq.bed",shell=True, universal_newlines=True)
    
    for ribo in ['A','C','G','T','N']:
         prefNTP='comp_nucl'+ribo
         temp=subprocess.check_output("grep -v chrM "+file+"/coordinate30/"+file+"-chrM_seq.bed | grep -i '\s"+ribo+"$'| wc -l", shell=True, universal_newlines=True)
         vars()[prefNTP]=int(temp.strip('\n'))
    for ribo in ['A','C','G','T','N']:
         prefNTP='comp_mito'+ribo
         temp=subprocess.check_output("grep chrM "+file+"/coordinate30/"+file+"-chromosomes_seq.bed | grep -i '\s"+ribo+"$'| wc -l", shell=True, universal_newlines=True)
         vars()[prefNTP]=int(temp.strip('\n'))
    '''
    mito_total=int(subprocess.check_output("cat "+file+"/coordinate30/"+file+"-chrM.coords.bed | wc -l", shell=True, universal_newlines=True))
    nucl_total=int(subprocess.check_output("cat "+file+"/coordinate30/"+file+"-chromosomes.coords.bed | wc -l", shell=True, universal_newlines=True))
    total_ribos=mito_total+nucl_total

    comp_mito=subprocess.check_output("cut -f2 "+file+"/composition30/"+file+"-chrM.counts.txt", shell=True, universal_newlines=True)
    comp_mito=comp_mito.split('\n')
    comp_mito=[int(i) for i in comp_mito[:-1]]
    comp_nucl=subprocess.check_output("cut -f2 "+file+"/composition30/"+file+"-chromosomes.counts.txt", shell=True, universal_newlines=True)
    comp_nucl=comp_nucl.split('\n')
    comp_nucl=[int(i) for i in comp_nucl[:-1]]


    #Paramter Set III:  Frequency of rN(A,C,G,U) in mitochondria and nucleus
    freq_mito=subprocess.check_output("cut -f2 "+file+"/composition30/"+file+"-chrM.frequencies.txt", shell=True, universal_newlines=True)
    freq_mito=freq_mito.split('\n')
    freq_nucl=subprocess.check_output("cut -f2 "+file+"/composition30/"+file+"-chromosomes.frequencies.txt", shell=True, universal_newlines=True)
    freq_nucl=freq_nucl.split('\n')
    
    #Parameter Set IV: Percentage of counts greater than one in nucleus and mitochondria
    positions_nucl=subprocess.check_output("grep -v 'chrM' "+file+"/coordinate30/"+file+"-chromosomes.counts.tab | wc -l", shell=True, universal_newlines=True)
    positions_nucl=int(positions_nucl.strip('\n'))
    positions1_nucl = subprocess.check_output("grep -v 'chrM' "+file+"/coordinate30/"+file+"-chromosomes.counts.tab | awk '$7>1' | wc -l", shell=True, universal_newlines=True)
    positions1_nucl=int(positions1_nucl.strip('\n'))
    
    positions_mito=subprocess.check_output("grep 'chrM' "+file+"/coordinate30/"+file+"-chrM.counts.tab | wc -l", shell=True, universal_newlines=True)
    positions_mito=int(positions_mito.strip('\n'))
    positions1_mito= subprocess.check_output("grep 'chrM' "+file+"/coordinate30/"+file+"-chrM.counts.tab | awk '$7>1' | wc -l", shell=True, universal_newlines=True)
    positions1_mito=int(positions1_mito.strip('\n'))
    
    pos_mito=str(positions1_mito)+"/"+str(positions_mito)
    pos_nucl=str(positions1_nucl)+"/"+str(positions_nucl)
    
    df_row={'Library':file,'Counts_All':total_ribos,'Counts_Nuclear':nucl_total,'Counts_Mito':mito_total,'Comp_nucl_rA':comp_nucl[1]/nucl_total*100,'Comp_nucl_rC':comp_nucl[2]/nucl_total*100,'Comp_nucl_rG':comp_nucl[2]/nucl_total*100,'Comp_nucl_rU':comp_nucl[3]/nucl_total*100,'Freq_nucl_rA':freq_nucl[0],'Freq_nucl_rC':freq_nucl[1],'Freq_nucl_rG':freq_nucl[2],'Freq_nucl_rU':freq_nucl[3],'Comp_mito_rA':comp_mito[0]/mito_total*100,'Comp_mito_rC':comp_mito[1]/mito_total*100,'Comp_mito_rG':comp_mito[2]/mito_total*100,'Comp_mito_rU':comp_mito[3]/mito_total*100,'Freq_mito_rA':freq_mito[0],'Freq_mito_rC':freq_mito[1],'Freq_mito_rG':freq_mito[2],'Freq_mito_rU':freq_mito[3], 'pos>1_nucleus/totpos':pos_nucl, 'pos>1_mito/totpos':pos_mito}
    new_row=pd.Series(df_row)
    df=df.append(new_row, ignore_index=True)

df.to_csv('QC.tsv', sep='\t', index=False)
print("Stats saved in QC.tsv")





