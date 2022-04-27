#!/bin/env Rscript
# author: ph-u
# script: selCheck.r
# desc: select 30 manual check data rows in cf425
# in: Rscript selCheck.r 1> ../data/selCheckRec.txt
# out: NA
# arg: 0
# date: 20220419

##### env #####
load("../data/cf425MedMic.rda")
#load("../data/cf425FULL.rda")
nSel = 30

##### sample across time series 2008-20 #####
repeat{sAm = sample(1:nrow(mIcro),nSel,replace=F);if(length(unique(mIcro$year[sAm]))==length(unique(mIcro$year))){break}} # random sample, confirm each data structure have at least one sample

##### manual check #####
for(i in 1:nSel){
	cat("#######\ni =",i,"\n")
#	print(mEdic[sAm[i],which(mEdic[sAm[i],]!=0)])
	print(mIcro[sAm[i],which(mIcro[sAm[i],]!=0)])
}#i = min(i+1,nSel)
