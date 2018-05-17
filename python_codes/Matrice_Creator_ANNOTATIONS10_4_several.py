# -*- coding: utf-8 -*-
"""
To convert an alignment output file into a matrice object
and normalise - correlation + DI + Eigen vector and virus concentration etc. 
10_2: without adding signal from capture output files.  
10_4: adding of the normalisation for the virus density signal: we divide by the general Hi-C coverage. 
November 2017: computation of sum of reads in A and B and sum of virus contact in A and B 
author: Axel KournaK
"""
import numpy as np
import matplotlib
from pylab import *
import os
import sys

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import scipy as sc
from scipy import signal
from scipy.linalg import eig
import scipy.stats as stats
import scipy.io as sio
from scipy.stats.stats import pearsonr

os.chdir("/home/axel/Bureau/z_python_scripts")
from histo_r import *
import scn_human
import distance_law_human
import directional_indice


nb_file_input = int(sys.argv[1])
BIN= int(sys.argv[2])    # Size (in bp) of the bin to build the contacts map 
virus_chosen = sys.argv[3]
name_bank = sys.argv[4]
indice= int(sys.argv[5])

#mat={};   # dictionary object to keep in memory all the contacts
mat = {}
maxi={};  # dictionary object to keep in memory the maximum reached by the chromosomes (in bins)
nbins_chr={};  # dictionary object to keep in memory the lenghts of the chromosomes (in bins)
i=0;


for fi in range(6,6+nb_file_input) :
    file_input = sys.argv[fi]
    print("sparsing file: ")
    print(file_input)

    # Creation of a dictionnary object to put into memory the cxontacts data 
    with open(file_input) as f: # open the file for reading output alignment file
        for line in f: # iterate over each line
            i=i+1;
            if i % 1000000 == 0:
                print(str(i)+" lines parsed.")
            chr1, locus1, sens1,indice1, chr2, locus2, sens2,indice2 = line.split()  # split it by whitespace
            locus1=int(locus1);sens1=int(sens1);
            locus2=int(locus2);sens2=int(sens2); 
            bin1 = int(locus1 /  BIN);
            bin2 = int(locus2 /  BIN);
            key1=(chr1, bin1, chr2, bin2)
            key2=(chr2, bin2, chr1, bin1)
            
            if key1 in mat:
                mat[key1] += 1;
            else:
                mat[key1] = 1  
            if key2 in mat:
                mat[key2] += 1
            else:
                mat[key2] = 1           
            if chr1 not in maxi:
                maxi[chr1] =0
            if chr2 not in maxi:
                maxi[chr2] =0     
            if bin1 > maxi[chr1]:
                maxi[chr1] =  bin1
            if bin2 > maxi[chr2]:
                maxi[chr2] =  bin2
            
print("All Hi-C contacts in memory in a dictionary object.")
#-----------------------------------------------------------------------------------------------

#C2 = sio.loadmat('pokemon_map.mat')
#C = C2['mycmap']
#cm = mpl.colors.ListedColormap(C)

#virus_chosen="HBVayw"
#virus_chosen="adenovirus"

# Computation of the matrice with the chromosome chosen and the virus:
list_all_chrms=('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10',
'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX') 

indice=0        # indice of the eigen vector to plot (can be 1 sometimes)

colorA="crimson"
colorB="royalblue"

colorVA="orange"
colorVB="darkturquoise"

indice_A_genome=0
indice_B_genome=0

V_in_A_genome=0
V_in_B_genome=0

H_in_A_genome = 0 
H_in_B_genome = 0

BIN_MIN_CHR={}  # dictionary object to keep in memory the min bin of each chromosome
BIN_MAX_CHR={}  # dictionary object to keep in memory the max bin of each chromosome
N_BINS=0;
for chr in list_all_chrms:
    nbins_chr[chr] = maxi[chr]+1
    BIN_MIN_CHR[chr] = N_BINS
    N_BINS = N_BINS + nbins_chr[chr]
    BIN_MAX_CHR[chr] = N_BINS-1
    print(chr,nbins_chr[chr],BIN_MIN_CHR[chr],BIN_MAX_CHR[chr])

#-------------------------------------------------------------------------------
# Building of the matrice with the whole genome:
'''
bin_mat1=0
bin_mat2=0
MATRICE=np.zeros( (N_BINS,N_BINS) )

for chr1 in list_all_chrms :
    for bin1 in range(0,maxi[chr1]+1) :
        bin_mat2=0;
        for chr2 in list_all_chrms  :  
            for bin2 in range(0,maxi[chr2]+1) :
                key1=(chr1, bin1, chr2, bin2);
                if key1 in mat:
                    #mat[key1] = mat[key1]; 
                    MATRICE[bin_mat1,bin_mat2]= mat[key1];
                else:
                    #mat[key1] = 0;
                    MATRICE[bin_mat1,bin_mat2]= 0;
                bin_mat2 +=1
        bin_mat1 +=1
#
MN=scn_human.scn_func(MATRICE,10)     # Normalisation of the whole genome 
#
##imshow(MN**0.2,interpolation="none",cmap=cm)
##colorbar()
##plt.savefig("Genome_CM_"+name_bank+"_"+virus_chosen+"_"+str(int(indice))+".eps", dpi=600, format='eps')
#
#
MATRICE_log = log(MN)
min_mn =  log(min(MN[MN>0]))
MATRICE_log[np.isinf(MATRICE_log)] = min_mn
fig, ax = plt.subplots()
imshow(MATRICE_log, interpolation="none",cmap="Reds")
colorbar()
savefig("Genome_CM_"+name_bank+"_"+virus_chosen+"_"+str(int(indice))+".eps", dpi=600, format='eps')
close('all')
np.savetxt("Genome_CM_"+name_bank+"_"+virus_chosen+"_"+str(int(indice))+"_SCN.txt",MN)

'''
#----------------------------------------------------------------------------------------------------------

V0= np.zeros(N_BINS)
Virus_genome= np.zeros(N_BINS)

#  For each chrm:
for chr_chosen in list_all_chrms:    # for each chromosome 
    print(chr_chosen)
    map_rows = range(BIN_MIN_CHR[chr_chosen],BIN_MAX_CHR[chr_chosen]+1)
    map_columns = range(BIN_MIN_CHR[chr_chosen],BIN_MAX_CHR[chr_chosen]+1)
    
    # Histo data:
    df=loadtxt("/media/axel/RSG53/PIERRICK_DATA/ChIP-seq_Miseq_H3K4me3_PNAS/p1.SRR2002335.sam.30."+chr_chosen)
    #df=loadtxt("/media/axel/RSG53/RNAseq_Novogene/data_release/rawdata/our-alignement/PHH_WT_1_1.fq.sam.MQ30."+chr_chosen)
    #df=loadtxt("/media/axel/RSG53/RNAseq_Novogene/data_release/rawdata/our-alignement/PHH_NI_1_1.fq.sam.MQ30."+chr_chosen)
    #df=loadtxt("/media/axel/RSG53/PIERRICK_DATA/ChIP_seq_Broad_ChipSeq_HepG2_H3K9me3/SRR568329.fastq.sam.MQ30."+chr_chosen)
    #df=loadtxt("/media/axel/RSG53/PIERRICK_DATA/rnaseq_SRX2463768/SRR5145370_1.sam.MQ30."+chr_chosen)
    #df=loadtxt("/media/axel/RSG53/PIERRICK_DATA/rnaseq_SRX2463768/SRR5145382_1.sam.MQ30."+chr_chosen)
    #df=loadtxt("/media/axel/RSG53/PIERRICK_DATA/ChIP-seq_CTCF/chrm_file/p1.sam.MQ30."+chr_chosen)
    #hist, bins = np.histogram(df,bins= int(max(df)/BIN) +1 )
    #hist, bins = np.histogram(df,bins= int(max(df)/BIN) +1 )     #  WRONG implementation
    hist, bins = np.histogram(df,bins=range(0,int(max(df)),BIN))   #  new implementation 20 / 12 / 2017

    # Check for the maximum of bins:  
    list_chr = [chr_chosen]
    N_BINS=0;
    for chr in list_chr:
        nbins_chr[chr] = maxi[chr]+1
        N_BINS=N_BINS + nbins_chr[chr]
        print(chr,nbins_chr[chr])
    print(N_BINS)
    
    #-------------------------------------------------------------------------------
    # Conversion of the Matrix in an numpy array object:
    bin_mat1=0;
    bin_mat2=0;
    MATRICE=np.zeros( (N_BINS,N_BINS) )
    
    for chr1 in list_chr :
        for bin1 in range(0,maxi[chr1]) :
            bin_mat2=0;
            for chr2 in list_chr  :  
                for bin2 in range(0,maxi[chr2]) :
                    key1=(chr1, bin1, chr2, bin2);
                    if key1 in mat:
                        #mat[key1] = mat[key1]; 
                        MATRICE[bin_mat1,bin_mat2]= mat[key1];
                    else:
                        #mat[key1] = 0;
                        MATRICE[bin_mat1,bin_mat2]= 0;
                    bin_mat2 +=1
            bin_mat1 +=1
    
    # Plot of the raw matrice and savings:
    #imshow(MATRICE**0.2,interpolation="none")
    coverage = MATRICE.sum(axis=0)
    
    nb_reads_chr_chosen = sum(MATRICE)
    
    # histogram of number of reads:
    #histo_r(MATRICE.sum(axis=0),100)
    
    # Normalisation of the contacts map:

    mn=scn_human.scn_func(MATRICE,10)
    #imshow(mn**0.2,interpolation="none")
    ## Computation of correlation matrice
    d=distance_law_human.dist_law(mn)
    n1=mn.shape[0]
    ## Computation of genomic distance law matrice:
    MAT_DIST =  np.zeros((n1, n1));
    for i in range(0,n1) :
        for j in range(0,n1) :
            MAT_DIST[i,j] =  d[abs(j-i)]  
    ##imshow(MAT_DIST**0.2)
    ## detrendage:
    MAT2=mn/MAT_DIST
    MAT2[np.isnan(MAT2)] = 0.0
    
    mat_positif = MAT2[MAT2.sum(axis=0)>0,:]
    mat_positif = mat_positif[:,MAT2.sum(axis=0)>0]
    
    mc=np.corrcoef(mat_positif)
    mc[np.isnan(mc)] = 0.0
  
    # Computation of eigen vector: 
    (V,D) = sc.linalg.eig(mc)
    v1_p=D[:,indice]
    
    v1 = np.zeros(n1)
    v1[np.where(MAT2.sum(axis=0)>0)] = v1_p
        
    # in order to have the histogram and the eigen vector with the same length:    
    hist1= np.zeros(n1)
    hist1[np.where(hist>-1)] = hist
    
    mc2= np.zeros( (n1,n1) ) 
    map_rows2    = list(np.where(MAT2.sum(axis=0)>0)[0])
    map_columns2 = list(np.where(MAT2.sum(axis=0)>0)[0])
    mc2[np.ix_(map_rows2, map_columns2)] = mc
    
    PC=pearsonr(v1,hist1)[0]
    if PC <0 :
        v1=-v1
    
    V0[map_rows] = v1

    list_chr = [virus_chosen,chr_chosen]
    N_BINS=0;
    for chr in list_chr:
        nbins_chr[chr] = maxi[chr]+1
        N_BINS=N_BINS + nbins_chr[chr]
        print(chr,nbins_chr[chr]) ;
    print(N_BINS)
    bin_mat1=0;
    bin_mat2=0;
    MATRICE=np.zeros( (N_BINS,N_BINS) )
    
    for chr1 in list_chr :
        print chr1
        for bin1 in range(0,nbins_chr[chr1]) :
            bin_mat2=0
            for chr2 in list_chr :  
                for bin2 in range(0,nbins_chr[chr2]) :
                    key1=(chr1, bin1, chr2, bin2);
                    if key1 in mat:
                        mat[key1] = mat[key1]; 
                        MATRICE[bin_mat1,bin_mat2]= mat[key1];
                    else:
                        mat[key1] = 0;
                        MATRICE[bin_mat1,bin_mat2]= 0;
                    bin_mat2 +=1
            bin_mat1 +=1
    #  Plot of density of the virus:
    virus_density = MATRICE[0,range(1,len(MATRICE[0,:]) )]

    # downloading of virus data:
    #df=loadtxt("/media/axel/RSG5/Next_seq_14decembre/Pierrick_data_december/seqs2/output_alignment_idpt_"+name_bank+".dat.indices.filtered.pcr4.inter.HBV."+chr_chosen)
    #hist2, bins = np.histogram(df,bins= int(max(df)/BIN) +1 )

    # in order to have the histogram and the eigen vector with the same length:    
    #hist21= np.zeros(n1)
    #hist21[np.where(hist>-1)] = hist2

    virus_density = virus_density

    Virus_genome[map_rows] = virus_density 
    BIN1000=BIN/1000
    nb_reads_virus = sum(virus_density)
    proportion_virus_chr_chosen =  nb_reads_virus / (nb_reads_chr_chosen + nb_reads_virus ) * 100 
        
    # Computation of proportions of virus reads in compartment A or B:
    indices_A = np.where( v1 > 0)  # number of bins attributed to compartment A 
    indices_B = np.where( v1 < 0)  # number of bins attributed to compartment B
    indice_A_genome= indice_A_genome + len(indices_A[0]) 
    indice_B_genome= indice_B_genome + len(indices_B[0])    
    
    V_in_A =  sum(virus_density[indices_A[0] ])  # number of reads from the virus in compartment A
    V_in_B =  sum(virus_density[indices_B[0] ])  # number of reads from the virus in compartment B
    H_in_A =  sum(coverage[indices_A[0]])
    H_in_B =  sum(coverage[indices_B[0]])  
    
    # Computation for the whole genome: 
    V_in_A_genome = V_in_A_genome + V_in_A 
    V_in_B_genome = V_in_B_genome + V_in_B

    H_in_A_genome = H_in_A_genome + H_in_A
    H_in_B_genome = H_in_B_genome + H_in_B
    
    # Binomial test: Perform a test that the probability of success is p.
    proportionA=float(len(indices_A[0]) ) / float((len(indices_A[0])+len(indices_B[0]) ) ) * 100
    proportionB=float(len(indices_B[0]) ) / float((len(indices_A[0])+len(indices_B[0]) ) ) * 100 
    
    pA = sc.stats.binom_test(V_in_A, V_in_A+V_in_B, p= float(len(indices_A[0]) ) / float((len(indices_A[0])+len(indices_B[0]) ) ) )
    pB = sc.stats.binom_test(V_in_B, V_in_A+V_in_B, p= float(len(indices_B[0]) ) / float((len(indices_A[0])+len(indices_B[0]) ) ) )
    
    pvA=V_in_A/(V_in_A+V_in_B)*100
    pvB=V_in_B/(V_in_A+V_in_B)*100

    # Plot of contacts Map and Eigen vector and virus density:
    matplotlib.rcParams.update({'font.size': 8});
    subplots_adjust(hspace= 0.24,top=0.97,bottom=0.09);
    plt.figure(num=None, figsize=(7.05, 9.92), dpi=80, facecolor='w', edgecolor='k')
    gs = gridspec.GridSpec(4, 1, height_ratios=[5.,0.6,0.6,0.6] )    
    ax0 = plt.subplot(gs[0])
    #ax0.imshow(mn**0.2,interpolation="none",cmap=cm)
    #ax0.imshow(mc2,interpolation="none")
    
    mn_log = log(mn)
    min_mn =  log(min(mn[mn>0]))
    mn_log[np.isinf(mn_log)] = min_mn    
    
    #ax0.imshow(mn**0.2,interpolation="none",cmap="afmhot_r")
    ax0.imshow(mn_log,interpolation="none",cmap="Reds")
    ax0.set_title("Contacts Map of "+ chr_chosen+" Nb chr: "+str(int(nb_reads_chr_chosen))+" Nb virus: "+str(int(nb_reads_virus))+" Pc: "+str(round(proportion_virus_chr_chosen,2) ) +"%" )
    ax0.set_xlim([0,shape(mn)[0] ])
    
    ax1 = plt.subplot(gs[1], sharex=ax0 )
    ax1.set_title("Transcriptome")
    width = 1.
    ax1.bar(range(0,len(hist)), hist, align='center', width=width, color='chocolate')
    
    ax2 = plt.subplot(gs[2], sharex=ax0 )    
    b1=0;
    b2=len(v1);
    ax2.set_xlim([b1, b2])
    ax2.set_title("Compartment attribution")
    ax2.bar(indices_A[0],v1[indices_A[0]], align='center', width=width, color=colorA)
    ax2.bar(indices_B[0],v1[indices_B[0]], align='center', width=width, color=colorB)
    #ax2.set_xlabel("Position along the genome - bins of " + str(BIN1000) + " kb" )
    ax2.set_ylabel("Eigen vector")
    
    ax3 = plt.subplot(gs[3], sharex=ax0 )
    b1=0;
    b2=len(virus_density)
    ax3.set_xlim([b1, b2])
    
    #ax3.set_title("Virus density")
    ax3.set_ylabel("Contact Score"+"\n"+"of "+ virus_chosen)
    ax3.set_xlabel("Position along the chromosome - bins of " + str(BIN1000) + " kb" )
    ax3.bar(indices_A[0],virus_density[indices_A[0] ]  , align='center', width=width, color=colorVA)
    ax3.bar(indices_B[0],virus_density[indices_B[0] ]  , align='center', width=width, color=colorVB)
    
    # Adding text of important informations :
    figtext(0.125,0.1735,  "Number of bins in compartment A: "+str(int(len(indices_A[0])))+", proportion: "+str(int(proportionA))+"%",color=colorA,fontweight="black")
    figtext(0.125,0.1585,"Number of bins in compartment B: "+str(int(len(indices_B[0])))+", proportion: "+str(int(proportionB))+"%",color=colorB,fontweight="black")
    
    figtext(0.125,0.035, "Number of "+virus_chosen+" reads in compartment A: "+str(int(V_in_A))+", proportion: "+str(pvA)+"%",color=colorVA,fontweight="extra bold")
    figtext(0.125,0.020, "Number of "+virus_chosen+" reads in compartment B: "+str(int(V_in_B))+", proportion: "+str(pvB)+"%",color=colorVB,fontweight="extra bold") 
    figtext(0.125,0.005,"Pvalue (from Binomial test): "+str(round(pA,4)),color='black',fontweight="extra bold")
    
    # Savings of pictures files:
    matplotlib.rcParams.update({'font.size': 8});
    subplots_adjust(hspace= 0.24,top=0.97,bottom=0.09);
    plt.savefig(chr_chosen+"_CM_"+name_bank+"_"+virus_chosen+"_"+str(int(indice))+".pdf", dpi=600, format='pdf')
    np.savetxt(chr_chosen+"_CM_"+name_bank+"_"+virus_chosen+"_"+str(int(indice))+".txt",mn)
    np.savetxt(chr_chosen+"_indicesA_"+name_bank+"_"+virus_chosen+"_"+str(int(indice))+".txt",indices_A)
    np.savetxt(chr_chosen+"_indicesB_"+name_bank+"_"+virus_chosen+"_"+str(int(indice))+".txt",indices_B)
    np.savetxt(chr_chosen+"_virus_pos_"+name_bank+"_"+virus_chosen+"_"+str(int(indice))+".txt",virus_density)
    print(chr_chosen+" was saved!")
    close('all')
    
    
# Genomic Results:
print(indice_A_genome)
print(indice_B_genome)

print(V_in_A_genome)
print(V_in_B_genome)

# Plot of general eigen vector:
fig2, ax2 = plt.subplots()
indices_A = np.where( V0 > 0)  # number of bins attributed to compartment A 
indices_B = np.where( V0 < 0)  # number of bins attributed to compartment B

b1=0;
b2=len(V0);
ax2.set_xlim([b1, b2])
ax2.set_title("Compartment attribution")
ax2.bar(indices_A[0],V0[indices_A[0]], align='center', width=width, color=colorA)
ax2.bar(indices_B[0],V0[indices_B[0]], align='center', width=width, color=colorB)
ax2.set_ylabel("Eigen vector")
ax2.set_xlabel("Position along the genome")
savefig("EigenVector_"+name_bank+"_"+virus_chosen+"_"+str(int(indice))+".eps", dpi=600, format='eps')
close('all')


# Plot of genomic virus:
fig3, ax3 = plt.subplots()
ax3.set_title("Virus density")
indices_A = np.where( V0 > 0)  # number of bins attributed to compartment A 
indices_B = np.where( V0 < 0)  # number of bins attributed to compartment B
ax3.bar(indices_A[0],Virus_genome[indices_A], align='center', width=width, color=colorVA)
ax3.bar(indices_B[0],Virus_genome[indices_B], align='center', width=width, color=colorVB)
ax3.set_ylabel("No. of reads")
ax3.set_xlabel("Position along the genome")
savefig("Virus_"+name_bank+"_"+virus_chosen+"_"+str(int(indice))+".eps", dpi=600, format='eps')
close('all')


# Plot of the results for the whole genome: 
print('Number of reads of the virus in A ')
print(V_in_A_genome)
print('Number of reads of the virus in B ')
print(V_in_B_genome)

print("nb_reads_human_in_A")
print(H_in_A_genome)
print("nb_reads_human_in_B")
print(H_in_B_genome)

# Fisher test on the whole genome:
oddsratio, pvalue = stats.fisher_exact([[V_in_A_genome, V_in_B_genome ], [H_in_A_genome, H_in_B_genome]])
print("pvalue Fisher test:")
print(pvalue)

# New barplot: 
data2 = [H_in_A_genome*1./ (H_in_A_genome+H_in_B_genome) , H_in_B_genome*1./ (H_in_A_genome+ H_in_B_genome) ]
data1 = [V_in_A_genome*1./ (V_in_A_genome+V_in_B_genome) , V_in_B_genome*1./ (V_in_A_genome+ V_in_B_genome) ]

propA_H = H_in_A_genome*1./ (H_in_A_genome + H_in_B_genome)
propA_V = V_in_A_genome*1./ (V_in_A_genome + V_in_B_genome)

fig, ax = plt.subplots()
width = 0.35 
ax.bar(0.5, data2[0],  width, color=colorVA)
ax.bar(0.5, data2[1], width, color=colorVB,bottom=data2[0])
ax.bar(1.5, data1[0],  width, color=colorVA)
ax.bar(1.5, data1[1], width, color=colorVB,bottom=data1[0])
ax.set_xticks( [0.5 ,1.5])
ax.set_xticklabels(('All \n contacts','HBV-human \n contacts'), rotation=45)
ax.set_title("Proportion of contacts in A and B"+"\n"+"Pvalue (Fisher test):"+str(pvalue,) )
text(0.5+width/2,0.2,str(round(propA_H,2)))
text(1.5+width/2,0.2,str(round(propA_V,2)))
savefig("Barplot2_"+name_bank+"_"+virus_chosen+"_"+str(int(indice))+".eps", dpi=600, format='eps')


# Bar plot of compartment: 

data1 = [indice_A_genome*1./ (indice_A_genome+indice_B_genome) , indice_B_genome*1./ (indice_A_genome+indice_B_genome)]
propA = indice_A_genome*1./ (indice_A_genome+indice_B_genome) 
fig, ax = plt.subplots()
width = 0.35 
ax.bar(0.5, data1[0], width, color=colorA)
ax.bar(0.5, data1[1], width, color=colorB, bottom=data1[0])
print(indice_A_genome)
print(indice_B_genome)
print(propA)
text(0.5,0.2,str(propA))
ax.set_title("Proportion of compartments A and B")
plt.gca().xaxis.set_major_locator(plt.NullLocator())
savefig("Barplot_compartments_"+name_bank+"_"+virus_chosen+"_"+str(int(indice))+".eps", dpi=600, format='eps')

close('all')




