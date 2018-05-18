#!/usr/bin/env Rscript

# Script to use the package hiCseq to identify borders in human contact maps - 24 / 02 / 2017

args = commandArgs(trailingOnly=TRUE)

library(HiCseg)
library(fields)

# m=read.table('/media/axel/RSG5/PLOTS_VIRUS/plots_concatenated/BC108/chr1_CM_BC108_concat_HBVayw_0.txt')
m=read.table(args[1])
# image.plot(m**0.2)

d=dim(m)
d=d[1]
n = d               #  size of the matrice 
Kmax=  round(n/5)   #  maximum number of borders

res=HiCseg_linkC_R(n,Kmax,"P",as.matrix(m), "Dplus")
#print(res)

v = res$t_hat
v = v[v>0]
v = head(v, -1)  # remove last element because it is the last element of the matrice 

write.table(file = paste("borders_POISSIAN_Dplus",args[2],".txt",sep="_"),  x = v, quote = FALSE,row.names = F,col.names = F)