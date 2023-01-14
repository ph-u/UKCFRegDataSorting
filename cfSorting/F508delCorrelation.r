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
samSize = read.csv("../graph/cftrm_samFreq.csv", header=T)

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

##### treatment distribution along time #####
x = c("Year","total","CFTRm","antimicrobials","chemicals","drug_interaction")
tReat = as.data.frame(matrix(nr=nrow(samSize), nc=length(x)))
colnames(tReat) = x;rm(x)
tReat$Year = samSize$year
tReat$total = apply(samSize[,-1],1,sum)
for(i in 3:ncol(tReat)){
        tReat[,i] = apply(samSize[,which(substr(colnames(samSize),i-1,i-1)==1)],1,sum)
};rm(i)

##### linear regression #####
glmm0 = data.frame("gLV_fit" = c(yCor[,13],yCor[,14],yCor[,15],yCor[,16],yCor[,17],yCor[,18],yCor[,19],yCor[,20]), "year" = samSize$year, "group" = rep(f1,each=nrow(yCor)))
for(i in 3:ncol(tReat)){
	glmm0[,ncol(glmm0)+1] = c(0,0,tReat[-c(1:2),i]/tReat$total[-c(1:2)]-tReat[1:(nrow(tReat)-2),i]/tReat$total[1:(nrow(tReat)-2)])
	colnames(glmm0)[ncol(glmm0)] = paste0(colnames(tReat)[i],"_rDiff")
};rm(i)
glmm0 = glmm0[!is.na(glmm0[,1]),]
yr2 = lm(glmm0[,1]~glmm0[,4]+glmm0[,5]+glmm0[,6]+glmm0[,7])
yr2.1 = lm(glmm0[,1]~glmm0[,6])

#glmm1 = data.frame("gLV_fit" = c(yCor[,13],yCor[,14],yCor[,15],yCor[,16],yCor[,17],yCor[,18],yCor[,19],yCor[,20]), "year" = samSize$year, "group" = rep(f1,each=nrow(yCor)))
#for(i in 3:ncol(tReat)){
#	glmm1[,ncol(glmm1)+1] = c(0,tReat[-1,i]/tReat$total[-1]-tReat[1:(nrow(tReat)-1),i]/tReat$total[1:(nrow(tReat)-1)])
#	colnames(glmm1)[ncol(glmm1)] = paste0(colnames(tReat)[i],"_rDiff")
#};rm(i)
#glmm1 = glmm1[!is.na(glmm1[,1]),]
#yr1 = lm(glmm1[,1]~glmm1[,4]+glmm1[,5]+glmm1[,6]+glmm1[,7])

##### correlation test #####
c1 = cor.test(glmm0[,4],glmm0[,1], method="spearman") # non-sig; gLV proportion vs CFTRm
c2 = cor.test(glmm0[,5],glmm0[,1], method="spearman") # non-sig; gLV proportion vs antimicrobials
c3 = cor.test(glmm0[,6],glmm0[,1], method="spearman") # alm-sig; gLV proportion vs chemicals
c4 = cor.test(glmm0[,7],glmm0[,1], method="spearman") # non-sig; gLV proportion vs interactions

cAll = data.frame(fac=colnames(tReat)[-c(1:2)],rho=c(c1$estimate, c2$estimate, c3$estimate, c4$estimate), S=c(c1$statistic, c2$statistic, c3$statistic, c4$statistic), adj.p=c(c1$p.value, c2$p.value, c3$p.value, c4$p.value)*4)
