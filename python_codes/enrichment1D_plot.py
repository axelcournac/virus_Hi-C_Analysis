# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 15:35:52 2017
@author: Axel KournaK
To plot and study signals from PHH data to see if the HBV is located at specific places.
5: on average autour des positions de contacts du virus 
Janvier 2018: on réadapte un peu le script pour prendre directement les fichiers d'output d'alignement
histo_signals_PHH622_janvier2018.py
"""
import numpy as np
import matplotlib
from pylab import *
import pandas as pd
import matplotlib.gridspec as gridspec

#-----------------
# Adeno : J7                  
list_HiC_files = ['/media/axel/d0a28364-6c64-4f8e-9efc-f332d9a0f1a9/pcrfree/output_files_alignment_iterative_adenovirus/output_alignment_idpt_BC176_CGAT_WT3Adeno_fused5runs.dat.indices.filtered.pcr5',
                  '/media/axel/RSG53/Next_seq_15_november2017/output/output_alignment_idpt_BC108_TGGT_PHH345_J7_Adeno.dat.indices.filtered.pcr5']                 

list_HiC_files = ['/media/axel/d0a28364-6c64-4f8e-9efc-f332d9a0f1a9/pcrfree/output_files_alignment_iterative/output_alignment_DADE.BC164.dat.indices.filtered.pcr5',
                  '/media/axel/d0a28364-6c64-4f8e-9efc-f332d9a0f1a9/pcrfree/output_files_alignment_iterative/output_alignment_DADE.BC176_CGAT_WT2_fused20runs.dat.indices.filtered.pcr5',
                  '/media/axel/RSG5/Next_seq_15_november2017/output/output_alignment_idpt_BC176_CGAT_PHH399_WT_fused_4captures.dat.indices.filtered.pcr5.pcr5',
                  '/media/axel/RSG5/Next_seq_15_november2017/output/output_alignment_idpt.BC164_GTGT_fused3captures.dat.indices.filtered.pcr5.pcr5'] 


virus_chosen = 'HBVayw'
virus_chosen = 'adenovirus'


#------------------------------------------------------------------------------
list_all_chrms= ("chr1","chr10","chr11","chr11_gl000202_random","chr12","chr13","chr14","chr15","chr16",
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
    da = df[ (df[0] == virus_chosen) & ( df[7].isin(list(list_all_chrms) )   )]
    db = df[ (df[7] == virus_chosen) & ( df[0].isin(list(list_all_chrms) )   )]
    
    daa = pd.DataFrame( transpose([da[7], da[8] ] ) )
    shape(daa)
    dbb = pd.DataFrame( transpose([db[0], db[1] ] ) )
    shape(dbb)    
    
    d_virus = d_virus.append( daa.append(dbb) )
    shape(d_virus)
    
    # Data HiC to have the HiC coverage for the Null model:     
    da = df[ ( df[0].isin(list(list_all_chrms) ) ) & ( df[7].isin(list(list_all_chrms) )   ) ]
    daa = pd.DataFrame( transpose([db[0], db[1] ]) )
    dbb = pd.DataFrame( transpose([da[7], da[8]  ]) )
    
    d_human = d_human.append( daa.append(dbb) )
    shape(d_human)


# ---------------------------------------------------------------------------------------------------------------------- 
# !! Positif controls (ChIPseq data:)
d_virus = pd.read_table('/media/axel/RSG5/PIERRICK_DATA/ChIP-seq_HepG2_PolII_REP2/p1.sam.MQ30',header=None, delimiter=" ") # PolII
d_human = pd.read_table('/media/axel/RSG5/PIERRICK_DATA/ChIP_HepG2_INPUT/p1.sam.MQ30',header=None, delimiter=" ")  # INPUT 

shape(d_virus)
shape(d_human)

d_virus = d_virus.iloc[range(0,1500000)]

# ---------------------------------------------------------------------------------------------------------------------- 
# Selection for each chromosome and histogram 
# we take the positions of the signal for each chromosome 
ka={}
kb={}
for c in  list_all_chrms:
    print(c)
    ka[c] = d_virus.loc[(d_virus[0] == c)]
    kb[c] = d_human.loc[(d_human[0] == c)]


# Histogram of the signals 
BIN = 50000   #  bin for 1D histograms of signals
hist1={}
hist2={}
bins1={}
bins2={}

for c in list_all_chrms:
    print(c)
    va= np.array(ka[c][1])
    va = [ int(x) for x in va ]
    vb= np.array(kb[c][1])
    vb = [ int(x) for x in vb ]
    # alternative computation: 
    if len(va) > 0 : 
        hist1[c], bins1[c] = np.histogram(va,bins= range(0,int(max(vb)),BIN),range=None, normed=False, weights=None, density=None)
    else : 
        hist1[c] = 0
    if len(vb) > 0 :    
        hist2[c], bins2[c] = np.histogram(vb,bins= range(0,int(max(vb)),BIN),range=None, normed=False, weights=None, density=None)
    else : 
        hist2[c] = 0    
    
#    hist1[c], bins1[c] = np.histogram(va,bins= range(0,int(max(va)),BIN))
#    hist2[c], bins2[c] = np.histogram(vb,bins= range(0,int(max(vb)),BIN))    
    
    # Filtering of outliers (put them at the mediam):   FACULTATIVE !
#    h1 = hist1[c]
#    h1[ hist1[c] > (median(hist1[c]) + 4 * std(hist1[c]) ) ] = median(hist1[c])
#    hist1[c] = h1
#    h2 = hist2[c]
#    h2[ hist2[c] > (median(hist2[c]) + 4 * std(hist2[c]) ) ] = 0
#    h2[ hist2[c] < (median(hist2[c]) - 4 * std(hist2[c]) ) ] = 0
#    hist2[c] = h2
    
#-------------------------------------------------------------------------------------------------
#  INPUT: set of positions where we compute the average signal: 
pos_virus = pd.read_table('/media/axel/RSG53/R_codes_virus/data_UCSC_hg19/cpgIslandExt.txt3',header=None, delimiter=" ")      #  CpG 
pos_virus = pd.read_table('/media/axel/RSG53/R_codes_virus/data_UCSC_hg19/refGene_hg19_TSS.bed2',header=None, delimiter=" ")  # TSS sites
pos_virus = pd.read_table('/media/axel/RSG53/R_codes_virus/data_enhancers/hepatocytes-DS32057A.peaks.fdr0.01.hg19.bed3',header=None, delimiter=" ") # enhancers
pos_virus = pd.read_table('/media/axel/RSG53/R_codes_virus/data_UCSC_hg19/wa.HepG2.rep-1.J7.hg19.pks.bed',header=None, delimiter="\t")   # origins
pos_virus = pd.read_table('/media/axel/RSG53/CTCF_peacks/SL3582_Peaks.bed.broadPeak',header=None, delimiter="\t")  # CTCF

pos_virus = pd.read_table('/media/axel/RSG53/R_codes_virus/GRList_exons.bed',header=None, delimiter="\t")  # exons
pos_virus = pd.read_table('/media/axel/RSG53/R_codes_virus/GRList_introns.bed',header=None, delimiter="\t")  # introns

pos_virus = pd.read_table('/media/axel/RSG53/R_codes_virus/data_UCSC_hg19/laminB1Lads.txt',header=0, delimiter="\t")  # lamins

#pos_virus = pd.read_table('/media/axel/RSG5/CONCATENATED_PIERRICK_DATA_ITERATIV/output_alignment_idpt_BC164_BC176_cap_BC176_cap_BC164_cap_BC176.dat_capBC164_capBC176.indices.filtered.pcr5.HBVpos',header=None, delimiter=" ")
#pos_virus = pd.read_table('/media/axel/RSG5/R_codes_virus/data_UCSC_hg19/wa.HepG2.rep-1.J7.hg19.pks.bed',header=None, delimiter="\t")   # origins                   

#pos_virus = pd.read_table('/home/axel/Bureau/comparison_contacts_VIRUS_vs_transcriptome/VO_enriched_HBV_cpG_sorted_thr1.txt.th1',header=None, delimiter=" ")
#pos_virus = pd.read_table('/media/axel/RSG5/R_codes_virus/data_UCSC_hg19/tRNA_genes.txt3',header=None, delimiter=" ") 
#pos_virus = pd.read_table('/media/axel/RSG5/R_codes_virus/data_UCSC_hg19/enhancers_promoters_free.bed',header=None, delimiter="\t")    # enhancers promoters free 
#pos_virus = pd.read_table('/media/axel/RSG5/PIERRICK_DATA/ChIP-seq_HepG2_CTCF_REP1/GSM748538_HepG2_UTAChIPSeqCTCF_align_rep1.bed',header=None, delimiter="\t")    # CTCF

len(pos_virus)
pos_virus_uniq = pos_virus.drop_duplicates()  # to remove duplicate entries
len(pos_virus_uniq)
pos_virus = pos_virus_uniq 

# ------------------------------------------------------------------------------

#    b = np.array( (b[1]+b[2])/2 )   # with centering for CpG/CTCF notably   b = np.array( b[1] ) 
    #b = np.array( (b[5]+b[8])/2 )   # with centering for exons/introns notably 

area = 50 #  number of bins to compute the signal around the positions
signal_averaged1 = np.zeros(area*2 +1) 
signal_averaged2 = np.zeros(area*2 +1)
signal_occ       = np.zeros(area*2 +1)

positions_object=0
for c in list_all_chrms :
    print(c)
    b = pos_virus.loc[(pos_virus[0] == c)]
    b = np.array( b[1] )   #   TSS
    for i in range(0, len(b) ) :
        #pos_i = b[i] * BIN_MAT 
        pos_binned =  int(b[i]/ BIN)
        pi=0
        if size(hist1[c]) > 1 and size(hist2[c]) > 1: 
            for p in range(pos_binned-area, pos_binned+area+1):
                if p >=0 and p < len(hist1[c]) and p < len(hist2[c]) :
                    if hist2[c][p] > 0 :
                        signal_averaged1[pi] = signal_averaged1[pi] + hist1[c][p]
                        signal_averaged2[pi] = signal_averaged2[pi] + hist2[c][p]
                        signal_occ[pi] = signal_occ[pi] +1
                pi=pi+1
            positions_object += 1
print(positions_object)

#------------------------------------------------------------------------------
area = 50 #  number of bins to compute the signal around the positions
signal_averaged1 = np.zeros(area*2 +1) 
signal_averaged2 = np.zeros(area*2 +1)
signal_occ       = np.zeros(area*2 +1)

positions_object=0
for c in list_all_chrms :
    print(c)
    b = pos_virus.loc[(pos_virus['chrom'] == c)]
    b = np.array( (b['chromStart']+b['chromEnd'])/2 )   # 
    for i in range(0, len(b) ) :
        #pos_i = b[i] * BIN_MAT 
        pos_binned =  int(b[i]/ BIN)
        pi=0
        if size(hist1[c]) > 1 and size(hist2[c]) > 1: 
            for p in range(pos_binned-area, pos_binned+area+1):
                if p >=0 and p < len(hist1[c]) and p < len(hist2[c]) :
                    if hist2[c][p] > 0 :
                        signal_averaged1[pi] = signal_averaged1[pi] + hist1[c][p]
                        signal_averaged2[pi] = signal_averaged2[pi] + hist2[c][p]
                        signal_occ[pi] = signal_occ[pi] +1
                pi=pi+1
            positions_object += 1
print(positions_object)


#------------------------------------------------------------------------------
#  Plot of averaged signals :
#s1= (signal_averaged1/signal_averaged2) / signal_occ   #  meth 1
s1= (signal_averaged1/signal_averaged2)       #  Meth 2   [now the one we use and maybe the correct one] 
#s1= signal_averaged1   #  meth 
#s1= signal_averaged1/ signal_occ   #  meth 

# PLOT for binnage of 100 bp:
#plot(signal_averaged1, linewidth = 3.0, label="Ori",color="blue")

axvline(area, color='k', linestyle='--')
tick_locs=(0,area,area*2+1)
tick_lbls=('-2500 kb','','+2500 kb')
xticks(tick_locs, tick_lbls,fontsize=15,rotation=45)
ylabel("Averaged Virus Contact signal")
legend(loc=4)


axvline(100, color='k', linestyle='--')
tick_locs=(0,100,200)
tick_lbls=('-500 kb','TSS','+500 kb')
xticks(tick_locs, tick_lbls,fontsize=15,rotation=45)
ylabel("Averaged Virus Contact signal")

xlim(4950,5050)
ylim(0,0.00035)

#legend()

plot(signal_averaged2)
plot(signal_averaged1)

plot(signal_averaged1/signal_averaged2)





