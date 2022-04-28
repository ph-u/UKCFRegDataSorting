#!/bin/env Rscript
# author: ph-u
# script: micro2Genus.r
# desc: sort mIcro data in genus time-series
# in: Rscript micro2Genus.r
# out: genusCF425.csv
# arg: 0
# date: 20220428

##### f: import csv with colnames converted in standard format #####
iPt = function(x){
	x0 = read.csv(paste0("../data/",x,".csv"), header=T, stringsAsFactors=F)
	colnames(x0)[grep("[.][.]",colnames(x0))] = gsub("[.]$",")",gsub("[.][.]"," (",colnames(x0)[grep("[.][.]",colnames(x0))]))
	colnames(x0) = gsub("[.]"," ",colnames(x0))
	return(x0)
}

##### env #####
library(stringr)
mIcro = iPt("cf425Micro")#;mEdic = iPt("cf425Medic")
#load("../data/cf425MedMic.rda")

##### resort according to the list of genus #####
gEn = word(colnames(mIcro)[-c(1,2)],1)
uqGen = unique(gEn)
mcGen = as.data.frame(matrix(0,nr=nrow(mIcro),nc=2+length(uqGen)))
mcGen[,1:2] = mIcro[,1:2];colnames(mcGen) = c(colnames(mIcro)[1:2],uqGen)
for(i in 1:length(gEn)){
	tMp = grep(gEn[i],colnames(mcGen))
	mcGen[,tMp] = mcGen[,tMp] + mIcro[,2+i]
};for(i in 3:ncol(mcGen)){mcGen[,i] = ifelse(mcGen[,i]>0,1,0)}
## info check
#tMp = sample(1:nrow(mIcro),30,replace=F);i=1
#mIcro[tMp[i],which(mIcro[tMp[i],]!=0)];mcGen[tMp[i],which(mcGen[tMp[i],]!=0)];i = min(i+1,length(tMp))

##### time-series conversion #####
yR = unique(mIcro$year);yR = yR[order(yR)]
mcGenTS = as.data.frame(matrix(0,nr=length(yR),nc=ncol(mcGen)-1))
colnames(mcGenTS) = colnames(mcGen)[-1]
mcGenTS$year = yR
for(i0 in 1:nrow(mcGenTS)){
	i1 = mcGen[which(mcGen$year==mcGenTS$year[i0]),-c(1:2)]
	mcGenTS[i0,-1] = apply(i1,2,sum)/nrow(i1)*100
};rm(i0,i1)

##### export #####
write.csv(mcGenTS,"../data/genusTimeSeries.csv", row.names=F)
