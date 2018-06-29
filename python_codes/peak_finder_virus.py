# -*- coding: utf-8 -*-
"""
Created on Tue Jun 12 10:40:53 2018
@author: axel KournaK 
To find peacks of contacts of virus (adenovirus) and extract the DNA sequence.
We choose a threshold so that: contats signal of virus / Hi-C coverage > threshold 
"""
import numpy as np
import matplotlib
from pylab import *
import pandas as pd
import matplotlib.gridspec as gridspec
import random
from Bio import SeqIO
from Bio.SeqUtils import GC
from Bio import Restriction
from Bio.Restriction import *
import os

# Adeno : J7                  
list_HiC_files = ['/media/axel/d0a28364-6c64-4f8e-9efc-f332d9a0f1a9/pcrfree/output_files_alignment_iterative_adenovirus/output_alignment_idpt_BC176_CGAT_WT3Adeno_fused5runs.dat.indices.filtered.pcr5.positions',
                  '/media/axel/RSG53/Next_seq_15_november2017/output/output_alignment_idpt_BC108_TGGT_PHH345_J7_Adeno.dat.indices.filtered.pcr5.positions'] 

list_HiC_files = ['/media/axel/d0a28364-6c64-4f8e-9efc-f332d9a0f1a9/pcrfree/output_files_alignment_iterative/output_alignment_DADE.BC164.dat.indices.filtered.pcr5',
                  '/media/axel/d0a28364-6c64-4f8e-9efc-f332d9a0f1a9/pcrfree/output_files_alignment_iterative/output_alignment_DADE.BC176_CGAT_WT2_fused20runs.dat.indices.filtered.pcr5',
                  '/media/axel/RSG53/Next_seq_15_november2017/output/output_alignment_idpt_BC176_CGAT_PHH399_WT_fused_4captures.dat.indices.filtered.pcr5.pcr5',
                  '/media/axel/RSG53/Next_seq_15_november2017/output/output_alignment_idpt.BC164_GTGT_fused3captures.dat.indices.filtered.pcr5.pcr5'] 
#------------------------------------------------------------------------------
virus_chosen = 'adenovirus'
#virus_chosen = 'HBVayw'

list_all_chrms= ("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14",
"chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX")

d_virus = pd.DataFrame()
d_human = pd.DataFrame()

for fi in  list_HiC_files:
    #file_input = sys.argv[fi]
    file_input = fi
    print("sparsing file: ")
    print(file_input)
    df=pd.read_table(file_input,header=None, delimiter=" ")
    shape(df)
    df.head()
    
    # Virus contacts on the human genome :   
    da = df[ (df[0] == virus_chosen) & ( df[4].isin(list(list_all_chrms) )   )]
    db = df[ (df[4] == virus_chosen) & ( df[0].isin(list(list_all_chrms) )   )]
    
    daa = pd.DataFrame( transpose([da[4], da[10], da[11]  ]) )
    shape(daa)
    dbb = pd.DataFrame( transpose([db[0], db[8],  db[9]]) )
    shape(dbb)    
    
    d_virus = d_virus.append( daa.append(dbb) )
    shape(d_virus)
    
    # Data HiC to have the HiC coverage for the Null model:     
    da = df[ ( df[0].isin(list(list_all_chrms) ) ) & ( df[4].isin(list(list_all_chrms) )   ) ]
    daa = pd.DataFrame( transpose([da[0], da[8],  da[9] ]) )
    dbb = pd.DataFrame( transpose([da[4], da[10], da[11]  ]) )
    
    d_human = d_human.append( daa.append(dbb) )
    shape(d_human)
    
# Histogram by chromosome:

ka={}
kb={}
for c in  list_all_chrms:
    print(c)
    ka[c] = d_virus.loc[(d_virus[0] == c)]
    kb[c] = d_human.loc[(d_human[0] == c)]

# Histogram of the signals 
BIN = 2000   #  bin for 1D histogram
hist1={}
hist2={}
bins1={}
bins2={}
list_nb_peacks = []
list_len_chr = []

bins_peaks={}
threshold = 5.0
fout = open("positions_peaks_all_chrs.txt","w")

for c in list_all_chrms:
    print(c)
    va= np.array(ka[c][1])
    va = [ int(x) for x in va ]
    vb= np.array(kb[c][1])
    vb = [ int(x) for x in vb ]
    hist1[c], bins1[c] = np.histogram(va,bins= range(0,int(max(vb)),BIN),density="True")
    hist2[c], bins2[c] = np.histogram(vb,bins= range(0,int(max(vb)),BIN),density="True")
    
    # Filtering of bins enriched with virus 
    h1 = hist1[c] / hist2[c]
    h1[np.isnan(h1)] = 1.0
    b1 = bins1[c] 
    bins_peaks[c]  = b1[ h1 > threshold ]
    len(bins_peaks[c]  )
    list_nb_peacks.append(len(bins_peaks[c]  ))
    list_len_chr.append(int(max(vb)))  
    for p in bins_peaks[c] :
        print(p) 
        fout.write(str(c)+"\t"+str(p)+"\t"+str(p+BIN)+"\t"+"+"+"\n")
    
    record = SeqIO.read(open("/media/axel/RSG53/data_celine/human_genome_Virus/"+c+".fa"), "fasta")
    fh1 = open("Sequences_peaks_"+c+".txt","w")
    for p in bins_peaks[c] :
        print(p)
        seq = record.seq[ p : p+BIN ]
        fh1.write(">sequence_"+ c +"_" + str(p) + "\n")
        fh1.write(str(seq) )
        fh1.write("\n")
        
fout.close()

    #np.savetxt('pos_peaks'+c+".txt", bins_peaks[c], fmt='%d')    
    
# various plots:
plot(b1[:-1],hist1[c] / hist2[c])  
plot(bins_peaks[c],bins_peaks[c]*0.,'o')
    
plot(list_nb_peacks,'o')  
plot(list_len_chr,'o') 

plot(list_len_chr,list_nb_peacks,'o') 


# Extraction of sequences: 
