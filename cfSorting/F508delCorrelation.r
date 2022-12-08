#!/bin/env Rscript
# author: ph-u
# script: F508delCorrelation.r
# desc: Correlation between CFTRm yearly diff with gLV model ability
# in: Rscript F508delCorrelation.r
# out: F508_Correlation.csv
# arg: 0
# date: 20221208

##### import #####
load("../data/cf425MedMic.rda");
rEf = read.csv("../data/cftrm_interaction.csv", header=T,stringsAsFactors=F);
anF = read.csv("../data/antiInfectives.csv", header=T, stringsAsFactors=F);
muT = read.csv("../data/mutPWCF.csv", header=T, stringsAsFactors=F);

##### env #####
yR = unique(mIcro$year);yR = yR[order(yR)]; # timeline

## CFTRm popularity by year
cftrM = c();for(i in strsplit(unique(rEf$cftrm),";")){cftrM = c(cftrM,i)};rm(i);cftrM = tolower(unique(cftrM)); # caftors list
mEd = ifelse(rowSums(mEdic[,which(colnames(mEdic) %in% cftrM)])>0,1,0); # 0 = no CFTRm

## plot settings
sEq = c(2,9,3,4,5,6,7,1,8)
cBp = palette.colors(palette = "Okabe-Ito", alpha=1, recycle = T)[sEq]

## F508del
F508del = muT$regid_anon[which(muT[,3]==1)]
mIcro = mIcro[which(mIcro$regid_anon %in% F508del),]
mEdic = mEdic[which(mEdic$regid_anon %in% F508del),]

## F508del simulation result
nRep = 7
simO = 500
ptIN = "../F508del_data/"
ptOT = "../graph/"
f = list.files(ptIN,"-eco")
f0 = strsplit(f,"_")
f1 = c();for(i in 1:length(f0)){f1 = c(f1,f0[[i]][2])};rm(i);f1 = unique(f1) # get medication group

##### correlation data #####
cNam = c("Year","sampleSize","CFTRmCount","ratioDiff",paste0(paste0("G",1:length(f1),"-"),rep(c("mutualism","overall"),each=8)))
yCor = as.data.frame(matrix(nr=length(yR),nc=length(cNam)))
colnames(yCor) = cNam;yCor[,1] = yR;rm(cNam)
for(i in 1:nrow(yCor)){
	yCor[i,2:3] = c(sum(mEdic$year==yCor$Year[i]), sum(mEd[which(mEdic$year==yCor$Year[i])]))
	if(i>2){
		for(i0 in 5:(4+length(f1))){
			fNam = paste0(ptIN,f[grep(paste0("_",f1[i0-4],"_*",substr(yCor$Year[i-2],3,4),substr(yCor$Year[i],3,4)),f)])
			if(length(grep("[.]csv$",fNam))>0){
				a = read.csv(fNam, header=T, stringsAsFactors=F)
				if(i==3 & i0==5){sP = length(unique(paste0(a$category1,a$category2)))}
				a0 = a[grep("mutu",a$c1_is),]
				yCor[i,c(i0,i0+length(f1))] = c(sum(a0$count),sum(a$count))/(sP*nRep*simO)
			}
		};rm(i0,a,a0)
	}};rm(i)
yCor$ratioDiff[-1] = yCor$CFTRmCount[-1]/yCor$sampleSize[-1] - yCor$CFTRmCount[-nrow(yCor)]/yCor$sampleSize[-nrow(yCor)]
write.csv(yCor,"../data/F508_Correlation.csv", quote=F, row.names=F)

##### correlation test #####
cor.test(rep(yCor$ratioDiff,8),c(yCor[,5],yCor[,6],yCor[,7],yCor[,8],yCor[,9],yCor[,10],yCor[,11],yCor[,12]), method="spearman") # non-sig
cor.test(rep(yCor$ratioDiff,8),c(yCor[,13],yCor[,14],yCor[,15],yCor[,16],yCor[,17],yCor[,18],yCor[,19],yCor[,20]), method="spearman")
matplot(yCor$ratioDiff,yCor[,13:20])
