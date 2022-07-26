#!/bin/env Rscript
# author: ph-u
# script: medNameExtraction.r
# desc: get standardized medication and active ingredient names
# in: Rscript medNameExtraction.r
# out: ../data/drug.csv
# arg: 0
# date: 20220726

##### env #####
p="../data/";a = read.csv(paste0(p,"medName-qcList.csv"), header=T)[,-c(2,3,6)];

##### drug/active ingredient identification #####
nAm = read.table(text=unique(tolower(c(a$QC.standardization,a$func.component))), sep=";", fill=T); # ln2318
n0 = c();for(i in 1:ncol(nAm)){n0 = c(n0,nAm[,i])};n0 = unique(n0); # ln1475

##### export #####
write.table(c("id",n0), paste0(p,"drug.csv"), sep=",", quote=F, row.names=F, col.names=F);
