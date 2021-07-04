#!usr/bin/env python3
#Author: Deepali L. Kundnani
#Institution: Georgia Institute of Technology (Georgia Tech)
#This script is to be ran from the results folder of Ribosemap run to output unique read percentage
#Activate the conda enviroment for Ribosemap before running this script(or any that provides samtools and bedtools)

import subprocess
import pandas as pd

import glob
sample_names=glob.glob('FS*')
sample_names.sort()

import argparse
parser = argparse.ArgumentParser(description='Quality statitics of ribosemap runs for post Alignment, Coordinate, Composition modules')
parser.add_argument("-r", "--referencefasta", required=True, help="Specify the fasta file for reference of the sequenced Libraries. Make sure you have Ribosemap environment activated or have bedtools available to be used")
args= parser.parse_args()

df = pd.DataFrame(columns=['Library', 'raw_reads', 'alignment_rate','barcoded_percent','unique_percent','Counts_All','Counts_Nuclear','Counts_Mito','Comp_nucl_rA','Comp_nucl_rC','Comp_nucl_rG','Comp_nucl_rU','Freq_nucl_rA','Freq_nucl_rC','Freq_nucl_rG','Freq_nucl_rU','Comp_mito_rA','Comp_mito_rC','Comp_mito_rG','Comp_mito_rU','Freq_mito_rA','Freq_mito_rC','Freq_mito_rG','Freq_mito_rU','pos>1_nucleus/totpos', 'pos>1_mito/totpos'])

#subprocess.call("conda activate conda_env", shell=True, universal_newlines=True)
for file in sample_names:
    print("Analyzing "+file)
    #Parameter SetI:Processing Statitics: Raw reads from fastq files, Aignment rate, percent reads barcoded, unique read percent
    raw_reads = subprocess.check_output("grep ^@ ../readyfiles/cut_"+file+".f*q | wc -l", shell=True, universal_newlines=True)
    alignment_rate = subprocess.check_output("tail -1 "+file+"/alignment/alignment.log", shell=True, universal_newlines=True)
    alignment_rate=alignment_rate[0:6]
    barcode_reads=subprocess.check_output("grep ^@ "+file+"/alignment/demultiplexed1.fq | wc -l", shell=True, universal_newlines=True)
    barcode_perc=int(barcode_reads)/int(raw_reads)*100
    barcode_perc=str(barcode_perc)+"%"
    numerator = subprocess.check_output("samtools view -c "+file+"/alignment/"+file+".bam", shell=True, universal_newlines=True)
    denominator = subprocess.check_output("samtools view -c -F 4 "+file+"/alignment/sorted.bam", shell=True, universal_newlines=True)
    unique_read_ratio=int(numerator)/int(denominator)*100
    unique_read_percent=str(unique_read_ratio)+"%"
    
    #Paramter Set II: Counts, Composition and Frequency in Nucleus and Mitochondria
    '''
    chrX_reads=subprocess.check_output("grep chrX "+file+"/coordinate30/"+file+".bed | wc -l", shell=True, universal_newlines=True)
    #parameter 9: Reads in Chr Y
    chrY_reads=subprocess.check_output("grep chrY "+file+"/coordinate30/"+file+".bed | wc -l", shell=True, universal_newlines=True)
    
    '''
    #Paramter set 10:  Perc of rN(A,C,G,U) in mitochondria and nucleus.
    subprocess.call("bedtools getfasta -s -bedOut -fi "+args.referencefasta+" -bed "+file+"/coordinate30/"+file+".bed > "+file+"/coordinate30/"+file+"_seq.bed",shell=True, universal_newlines=True)
    
    for ribo in ['A','C','G','T','N']:
         prefNTP='comp_nucl'+ribo
         temp=subprocess.check_output("grep -v chrM "+file+"/coordinate30/"+file+"_seq.bed | grep -i '\s"+ribo+"$'| wc -l", shell=True, universal_newlines=True)
         vars()[prefNTP]=int(temp.strip('\n'))
    for ribo in ['A','C','G','T','N']:
         prefNTP='comp_mito'+ribo
         temp=subprocess.check_output("grep chrM "+file+"/coordinate30/"+file+"_seq.bed | grep -i '\s"+ribo+"$'| wc -l", shell=True, universal_newlines=True)
         vars()[prefNTP]=int(temp.strip('\n'))
    
    nucl_total=comp_nuclA+comp_nuclC+comp_nuclG+comp_nuclT+comp_nuclN
    mito_total=comp_mitoA+comp_mitoC+comp_mitoG+comp_mitoT+comp_mitoN
    total_ribos=mito_total+nucl_total
    
    #Paramter Set III:  Frequency of rN(A,C,G,U) in mitochondria and nucleus
    freq_mito=subprocess.check_output("cut -f2 "+file+"/composition30/"+file+"-chrM.frequencies.txt", shell=True, universal_newlines=True)
    freq_mito=freq_mito.split('\n')
    freq_nucl=subprocess.check_output("cut -f2 "+file+"/composition30/"+file+"-chromosomes.frequencies.txt", shell=True, universal_newlines=True)
    freq_nucl=freq_nucl.split('\n')
    
    #Parameter Set IV: Percentage of counts greater than one in nucleus and mitochondria
    positions_nucl=subprocess.check_output("grep -v 'chrM' "+file+"/coordinate30/"+file+".counts.tab | wc -l", shell=True, universal_newlines=True)
    positions_nucl=int(positions_nucl.strip('\n'))
    positions1_nucl = subprocess.check_output("grep -v 'chrM' "+file+"/coordinate30/"+file+".counts.tab | awk '$7>1' | wc -l", shell=True, universal_newlines=True)
    positions1_nucl=int(positions1_nucl.strip('\n'))
    
    positions_mito=subprocess.check_output("grep 'chrM' "+file+"/coordinate30/"+file+".counts.tab | wc -l", shell=True, universal_newlines=True)
    positions_mito=int(positions_mito.strip('\n'))
    positions1_mito= subprocess.check_output("grep 'chrM' "+file+"/coordinate30/"+file+".counts.tab | awk '$7>1' | wc -l", shell=True, universal_newlines=True)
    positions1_mito=int(positions1_mito.strip('\n'))
    
    pos_mito=str(positions1_mito)+"/"+str(positions_mito)
    pos_nucl=str(positions1_nucl)+"/"+str(positions_nucl)
    
    df_row={'Library':file, 'raw_reads':raw_reads.strip('\n'), 'alignment_rate':alignment_rate,'barcoded_percent':barcode_perc,'unique_percent':unique_read_percent,'Counts_All':total_ribos,'Counts_Nuclear':nucl_total,'Counts_Mito':mito_total,'Comp_nucl_rA':comp_nuclA/nucl_total*100,'Comp_nucl_rC':comp_nuclC/nucl_total*100,'Comp_nucl_rG':comp_nuclG/nucl_total*100,'Comp_nucl_rU':comp_nuclT/nucl_total*100,'Freq_nucl_rA':freq_nucl[0],'Freq_nucl_rC':freq_nucl[1],'Freq_nucl_rG':freq_nucl[2],'Freq_nucl_rU':freq_nucl[3],'Comp_mito_rA':comp_mitoA/mito_total*100,'Comp_mito_rC':comp_mitoC/mito_total*100,'Comp_mito_rG':comp_mitoG/mito_total*100,'Comp_mito_rU':comp_mitoT/mito_total*100,'Freq_mito_rA':freq_mito[0],'Freq_mito_rC':freq_mito[1],'Freq_mito_rG':freq_mito[2],'Freq_mito_rU':freq_mito[3], 'pos>1_nucleus/totpos':pos_nucl, 'pos>1_mito/totpos':pos_mito}
    new_row=pd.Series(df_row)
    df=df.append(new_row, ignore_index=True)

df.to_csv('QC.tsv', sep='\t', index=False)
print("Stats saved in QC.tsv")





