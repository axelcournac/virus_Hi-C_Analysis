# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 15:34:26 2018
@author: Axel KournaK
To compute and plot Histone Marks or Transcription enrichment
"""
import numpy as np
import matplotlib
from pylab import *
import pandas as pd
import matplotlib.gridspec as gridspec
import operator
import json
import pybedtools
from pybedtools import BedTool
import matplotlib.ticker as ticker
import os
import time

os.chdir('/media/axel/RSG5/python_codes_virus/')
from histo_r import *

# 3C data:
nb_file_input = int(sys.argv[1])
virus_chosen = sys.argv[2]

list_all_chrms = ("chr1","chr10","chr11","chr11_gl000202_random","chr12","chr13","chr14","chr15","chr16",
"chr17","chr17_ctg5_hap1","chr17_gl000203_random","chr17_gl000204_random","chr17_gl000205_random","chr17_gl000206_random",
"chr18","chr18_gl000207_random","chr19","chr19_gl000208_random","chr19_gl000209_random","chr1_gl000191_random",
"chr1_gl000192_random","chr2","chr20","chr21","chr21_gl000210_random","chr22","chr3","chr4","chr4_ctg9_hap1",
"chr4_gl000193_random","chr4_gl000194_random","chr5","chr6","chr6_apd_hap1","chr6_cox_hap2","chr6_dbb_hap3","chr6_mann_hap4",
"chr6_mcf_hap5","chr6_qbl_hap6","chr6_ssto_hap7","chr7","chr7_gl000195_random","chr8","chr8_gl000196_random","chr8_gl000197_random",
"chr9","chr9_gl000198_random","chr9_gl000199_random","chr9_gl000200_random","chr9_gl000201_random","chrM","chrUn_gl000211",
"chrUn_gl000212","chrUn_gl000213","chrUn_gl000214","chrUn_gl000215","chrUn_gl000216","chrUn_gl000217","chrUn_gl000218",
"chrUn_gl000219","chrUn_gl000220","chrUn_gl000221","chrUn_gl000222","chrUn_gl000223","chrUn_gl000224","chrUn_gl000225",
"chrUn_gl000226","chrUn_gl000227","chrUn_gl000228","chrUn_gl000229","chrUn_gl000230","chrUn_gl000231","chrUn_gl000232",
"chrUn_gl000233","chrUn_gl000234","chrUn_gl000235","chrUn_gl000236","chrUn_gl000237","chrUn_gl000238","chrUn_gl000239",
"chrUn_gl000240","chrUn_gl000241","chrUn_gl000242","chrUn_gl000243","chrUn_gl000244","chrUn_gl000245","chrUn_gl000246",
"chrUn_gl000247","chrUn_gl000248","chrUn_gl000249","chrX","chrY")

def f_neg(x):
    if x >= 0:
        x = x
    else :
        x = 0
    return x 

#   PARAMETERS : 
virus_chosen = 'HBVayw'
virus_chosen = 'adenovirus'

print("Virus chosen")
print(virus_chosen)

# HBV 1rst patient:    
list_HiC_files1 = ['/media/axel/d0a28364-6c64-4f8e-9efc-f332d9a0f1a9/pcrfree/output_files_alignment_iterative/output_alignment_DADE.BC164.dat.indices.filtered.pcr5',
                  '/media/axel/d0a28364-6c64-4f8e-9efc-f332d9a0f1a9/pcrfree/output_files_alignment_iterative/output_alignment_DADE.BC176_CGAT_WT2_fused20runs.dat.indices.filtered.pcr5',
                  '/media/axel/RSG5/Next_seq_15_november2017/output/output_alignment_idpt_BC176_CGAT_PHH399_WT_fused_4captures.dat.indices.filtered.pcr5.pcr5',
                  '/media/axel/RSG5/Next_seq_15_november2017/output/output_alignment_idpt.BC164_GTGT_fused3captures.dat.indices.filtered.pcr5.pcr5'] 

# HBV 2d patient:  
list_HiC_files2 = ['/media/axel/RSG5/Next_seq_27_november_2017/pcrfree/indices1/output_alignment_idpt_BC182_AGAT_PHH342_WT-1.dat.indices1.pcr5',
                  '/media/axel/RSG5/Next_seq_27_november_2017/pcrfree/indices1/output_alignment_idpt_BC186_ACGT_PHH342_WT-2.dat.indices1.pcr5',
                  '/media/axel/RSG5/Next_seq_2_december_2017/axel_fastq/pcrfree/output_alignment_idpt_BC182_AGAT_PHH342_WT-1.dat.indices.filtered.pcr5',
                  '/media/axel/RSG5/Next_seq_2_december_2017/axel_fastq/pcrfree/output_alignment_idpt_BC186_ACGT_PHH342_WT-2.dat.indices.filtered.pcr5', 
                  '/media/axel/RSG5/Next_seq_1_december_2017/Pierrick/output_alignment_idpt_BC182_AGAT_PHH342_WT_1.dat.indices1.filtered.pcr52',
                  '/media/axel/RSG5/Next_seq_1_december_2017/Pierrick/output_alignment_idpt_BC186_ACGT_PHH342_WT_2.dat.indices1.filtered.pcr5']
           
list_HiC_files =  list_HiC_files1 +  list_HiC_files2       
           
# Adeno J7:               
list_HiC_files = ['/media/axel/d0a28364-6c64-4f8e-9efc-f332d9a0f1a9/pcrfree/output_files_alignment_iterative_adenovirus/output_alignment_idpt_BC176_CGAT_WT3Adeno_fused5runs.dat.indices.filtered.pcr5']                   

# Adeno J4:                  
list_HiC_files = ['/media/axel/d0a28364-6c64-4f8e-9efc-f332d9a0f1a9/pcrfree/output_files_alignment_iterative_adenovirus/output_alignment_idpt_BC172_CGGT_NI3Adeno_fused5runs.dat.indices.filtered.pcr5'] 
 
#------------------------------------------------------------------------------   
area_contact = 3500
d_virus = pd.DataFrame()
d_human = pd.DataFrame()

for fi in  list_HiC_files:
    #file_input = sys.argv[fi]
    file_input = fi
    print("sparsing file: ")
    print(file_input)
    df=pd.read_table(file_input,header=None, delimiter=" ")
    shape(df)
    
    # Virus contacts on the human genome :   
    da = df[ (df[0] == virus_chosen) & ( df[7].isin(list(list_all_chrms) )   )]
    db = df[ (df[7] == virus_chosen) & ( df[0].isin(list(list_all_chrms) )   )]
    
    daa = pd.DataFrame( transpose([da[7], da[8] ]) )
    shape(daa)
    dbb = pd.DataFrame( transpose([db[0], db[1] ]) )
    shape(dbb)    
    
    d_virus = d_virus.append( daa.append(dbb) )
    shape(d_virus)
    
    # Data HiC to have the HiC coverage for the Null model:     
    da = df[ ( df[0].isin(list(list_all_chrms) ) ) & ( df[7].isin(list(list_all_chrms) )   ) ]
    daa = pd.DataFrame( transpose([da[0], da[1] ]) )
    dbb = pd.DataFrame( transpose([da[7], da[8] ]) )
    
    d_human = d_human.append( daa.append(dbb) )
    shape(d_human)

# Conversion into a bedtool: 
contacts_virus_bed = BedTool("\n".join(("\t".join((i,str( f_neg(int(j)-area_contact)), str(j+area_contact)))) for i,j in zip(d_virus[0],d_virus[1]) ), from_string=True)

# control positif: 
#contacts_virus_bed = pybedtools.BedTool('/media/axel/RSG5/R_codes_virus/data_UCSC_hg19/refGene_hg19_TSS.bed')  # TSS sites
#contacts_virus_bed = (contacts_virus_bed.sort()).merge(d=0)

len(contacts_virus_bed)

# Chip-seq signals analysis : 

# Input data:
file_input = '/media/axel/RSG5/PIERRICK_DATA/Bernstein_HepG2_Control/SRR227552.fastq.sam.MQ30' 
file_input = '/media/axel/RSG5/PIERRICK_DATA/ChIP_HepG2_INPUT/p1.sam.MQ30'
df2=pd.read_table(file_input,header=None, delimiter=" ")
shape(df2)
input_signal_bed = BedTool("\n".join(("\t".join(( str(i),str(j), str(j) ))) for i,j in zip(df2[0],df2[1]) ), from_string=True)

# Computation of the sum of Input signal via the intersect function
intersection_virus_input = contacts_virus_bed.intersect(input_signal_bed)
sum_virus_input = intersection_virus_input.count()

# Chip-seq signals:
LIST_CHIP_seq = ['/media/axel/RSG5/PIERRICK_DATA/ChIP_Bernstein_HepG2_H3K4me3/SRR227563.fastq.sam.MQ30',
                 '/media/axel/RSG5/PIERRICK_DATA/Bernstein_HepG2_H3K4me2/SRR227466.fastq.sam.MQ30',
                 '/media/axel/RSG5/PIERRICK_DATA/Bernstein_HepG2_H4K20me1/SRR227468.fastq.sam.MQ30',
                 '/media/axel/RSG5/PIERRICK_DATA/Bernstein_HepG2_H3K36me3/SRR227447.fastq.sam.MQ30',
                 '/media/axel/RSG5/PIERRICK_DATA/Bernstein_HepG2_H3K27ac/SRR227575.fastq.sam.MQ30',
                 '/media/axel/RSG5/PIERRICK_DATA/Bernstein_HepG2_CTCF/SRR227363.fastq.sam.MQ30',
                 '/media/axel/RSG5/PIERRICK_DATA/Bernstein_HepG2_H3K27me3/SRR227598.fastq.sam.MQ30',
                 '/media/axel/RSG5/PIERRICK_DATA/Bernstein_HepG2_H3K79me2/SRR227354.fastq.sam.MQ30',
                 '/media/axel/RSG5/PIERRICK_DATA/ChIP_seq_Broad_ChipSeq_HepG2_H3K9me3/SRR568329.fastq.sam.MQ30']
                 
LIST_CHIP_seq = ['/media/axel/RSG5/PIERRICK_DATA/ChIP-seq_HepG2_PolII_REP2/p1.sam.MQ30',
                 '/media/axel/RSG5/PIERRICK_DATA/ChIP-seq_HepG2_cMYC_REP1/p1.sam.MQ30',
                 '/media/axel/RSG5/PIERRICK_DATA/ChIP_seq_Broad_ChipSeq_HepG2_H3K9me3/SRR568329.fastq.sam.MQ30']   
                 
LIST_CHIP_seq = ['/media/axel/RSG5/PIERRICK_DATA/ChIP_Bernstein_HepG2_H3K4me3/SRR227563.fastq.sam.MQ30',
                 '/media/axel/RSG5/PIERRICK_DATA/Bernstein_HepG2_H3K4me2/SRR227466.fastq.sam.MQ30',
                 '/media/axel/RSG5/PIERRICK_DATA/Bernstein_HepG2_H3K27ac/SRR227575.fastq.sam.MQ30',
                 '/media/axel/RSG5/PIERRICK_DATA/Bernstein_HepG2_H3K27me3/SRR227598.fastq.sam.MQ30',
                 '/media/axel/RSG5/PIERRICK_DATA/ChIP_seq_Broad_ChipSeq_HepG2_H3K9me3/SRR568329.fastq.sam.MQ30']              
                 
                 
len(LIST_CHIP_seq)
t0 =time.time()

S = list()
for chip_file in LIST_CHIP_seq :
    file_input = chip_file
    print(file_input)
    df1=pd.read_table(file_input,header=None, delimiter=" ")
    shape(df1)
    # Conversion into a bedtool: 
    HM_signal_bed = BedTool("\n".join(("\t".join(( str(i),str(j), str(j) ))) for i,j in zip(df1[0],df1[1]) ), from_string=True)
    
    # Computation of the sum of Histone Marks via the intersect function
    intersection_virus_HM = contacts_virus_bed.intersect(HM_signal_bed)
    sum_virus_HM = intersection_virus_HM.count()
    
    ratio_virus = float(sum_virus_HM) / float(sum_virus_input)
    
    nber_realisations = 100
    list_random_enrichment = list()
    # Random Null model :
    for r in range(nber_realisations) :
        print(r)
        contacts_random = d_human.sample(n= len(d_virus), replace=True, weights=None, random_state=None, axis=None)
        contacts_random_bed =  BedTool("\n".join(("\t".join(( str(i),str(j), str(j) ))) for i,j in zip(contacts_random[0],contacts_random[1]) ), from_string=True)
        # Computation of the sum of Histone Marks via the intersect function:
        intersection_random_HM = contacts_random_bed.intersect(HM_signal_bed)
        sum_random_HM = intersection_random_HM.count()
    
        # Computation of the sum of Input signal via the intersect function:
        intersection_random_input = contacts_random_bed.intersect(input_signal_bed)
        sum_random_input = intersection_random_input.count()
    
        list_random_enrichment.append( float(sum_random_HM) / float(sum_random_input) )
    
    ratio_random = mean(list_random_enrichment)
    
    # Computation of the log2 Ratio for this HM: 
    E = log(ratio_virus / ratio_random)
    
    S.append(E)
    
print time.time() - t0, "seconds wall time"
    
S2 = np.reshape(S, (1,len(S)))


# Plots : 
plt.figure()
ax = plt.gca();
imshow(S2,interpolation="none", cmap ="seismic",vmin=-1.,vmax=1.0)

# Minor ticks
ax.set_xticks(np.arange(-.5, shape(S2)[1], 1), minor=True);
ax.grid(which='minor', color='black', linestyle='-', linewidth=2)

ax.set_xticks( range(len(S)) )

#names_h = ['H3K9me3','H3K4me3']
#names_h = ['H3K4me3','H3K4me2','H4K20me1','H3K36me3','H3K27ac','CTCF','H3K27me3','H3K79me2','H3K9me3']
#names_h = ['PolII','cMYC','H3K9me3']
names_h = ['H3K4me3','H3K4me2','H3K27ac','H3K27me3','H3K9me3']
ax.set_xticklabels(names_h, rotation=45)

plt.yticks([], [])
ax.tick_params(axis='both', width = 0.)

cbar = plt.colorbar(shrink = 0.25,orientation="horizontal", ticks=[-1,0,1]);
cbar.set_label('Log2 enrichment')

title("Histone Modifications Enrichment (nb of contacts of "+str(virus_chosen)+" = "+str(len(contacts_virus_bed)) + ")" +
"\n"+"Total number of contacts = " + str(len(d_human)) + " Nber of contacts of virus in the control: "+str(len(intersection_virus_input))+" Nber input "+str(len(input_signal_bed))+"\n")
#"Files used = " + fi)

