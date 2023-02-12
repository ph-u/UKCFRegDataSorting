#!/bin/env Rscript
# author: ph-u
# script: F508delCorrelation.r
# desc: Correlation between CFTRm yearly diff with gLV model ability
# in: Rscript F508delCorrelation.r
# out: F508_Correlation.csv
# arg: 0
# date: 20221208

##### import #####
ptIN = "../p_raw_20230207/data/" # "../F508del_data/"
ptOT = "../p_raw_20230207/res/" # "../graph/"
samSize = read.csv(paste0(ptIN,"../F508dd_samFreq.csv"), header=T)
drugYrs = read.csv(paste0(ptOT,"medication.csv"), header=T)

## F508del simulation result
nRep = 7
simO = 500
f = list.files(ptIN,"-eco"); f = f[-grep("1101",f)] # group 1101 has 1 small sample size, unmodellable time-series

##### correlation data #####
cNam = c("group","start","end","gLV_fit",paste0(colnames(drugYrs)[-c(1:2)],"_rDiff"))
yCor = as.data.frame(matrix(0,nr=length(f),nc=length(cNam)))
colnames(yCor) = cNam
for(i in 1:length(f)){
	a = read.csv(paste0(ptIN,f[i]), header=T)
	if(i==1){
		nSp = unique(c(a$category1,a$category2)) # species list
		mPr = choose(length(nSp),2)+length(nSp) # pair of species
		eCo = unique(a$c1_is) # ecological interaction list
		}
	if(nrow(a)==0){yCor[i,4] = 0}else{yCor[i,4] = sum(a$count)}
## info record
	a0 = unlist(strsplit(f[i],"_"))
	yCor[i,1] = a0[2]
	yCor[i,2:3] = as.numeric(c(substr(a0[3],1,2),substr(a0[3],3,4)))+2000
	for(i0 in 5:ncol(yCor)){yCor[i,i0] = drugYrs[which(drugYrs[,1]==yCor[i,3]),i0-2]/drugYrs$total[which(drugYrs[,1]==yCor[i,3])]-drugYrs[which(drugYrs[,1]==yCor[i,2]),i0-2]/drugYrs$total[which(drugYrs[,1]==yCor[i,2])]} # endMedGroup/endTotal - startMedGroup/startTotal
};rm(i,i0,a,a0)
yCor[,4] = yCor[,4]/(nRep*simO*mPr*length(eCo)) # ratio of fit

##### correlation #####
cNam = c("factor","rho","S","p.val")
cAll = as.data.frame(matrix(0,nr=ncol(drugYrs)-2,nc=length(cNam)))
colnames(cAll) = cNam
cAll$factor = colnames(drugYrs)[-c(1:2)]
for(i in 1:nrow(cAll)){
	a = cor.test(yCor[,4], yCor[,i+4], method="spearman")
	cAll[i,-1] = c(a$estimate, a$statistic, a$p.value)
	};rm(i,a)
cAll$adj.p = p.adjust(cAll$p.val, "BH") # Benjamini-Hochberg false discovery rate

##### linear regression #####
summary(lm(yCor[,4]~yCor[,5]+yCor[,6]+yCor[,7]+yCor[,8]))

#cNam = c("Year","sampleSize","CFTRmCount","ratioDiff",paste0(paste0("G",1:length(f1),"-"),rep(c("mutualism","overall"),each=8)))
#yCor = as.data.frame(matrix(nr=length(yR),nc=length(cNam)))
#colnames(yCor) = cNam;yCor[,1] = yR;rm(cNam)
#for(i in 1:nrow(yCor)){
#	yCor[i,2:3] = c(sum(mEdic$year==yCor$Year[i]), sum(mEd[which(mEdic$year==yCor$Year[i])]))
#	if(i>2){
#		for(i0 in 5:(4+length(f1))){
#			fNam = paste0(ptIN,f[grep(paste0("_",f1[i0-4],"_*",substr(yCor$Year[i-2],3,4),substr(yCor$Year[i],3,4)),f)])
#			if(length(grep("[.]csv$",fNam))>0){
#				a = read.csv(fNam, header=T, stringsAsFactors=F)
#				if(i==3 & i0==5){sP = length(unique(paste0(a$category1,a$category2)))}
#				a0 = a[grep("mutu",a$c1_is),]
#				yCor[i,c(i0,i0+length(f1))] = c(sum(a0$count),sum(a$count))/(sP*nRep*simO)
#			}
#		};rm(i0,a,a0)
#	}};rm(i)
#yCor$ratioDiff[-1] = yCor$CFTRmCount[-1]/yCor$sampleSize[-1] - yCor$CFTRmCount[-nrow(yCor)]/yCor$sampleSize[-nrow(yCor)]
#write.csv(yCor,paste0(ptOT,"F508dd_Correlation.csv"), quote=F, row.names=F)

##### linear regression #####
#glmm0 = data.frame("gLV_fit" = c(yCor[,13],yCor[,14],yCor[,15],yCor[,16],yCor[,17],yCor[,18],yCor[,19],yCor[,20]), "year" = samSize$year, "group" = rep(f1,each=nrow(yCor)))
#for(i in 3:ncol(tReat)){
#	glmm0[,ncol(glmm0)+1] = c(0,0,tReat[-c(1:2),i]/tReat$total[-c(1:2)]-tReat[1:(nrow(tReat)-2),i]/tReat$total[1:(nrow(tReat)-2)])
#	colnames(glmm0)[ncol(glmm0)] = paste0(colnames(tReat)[i],"_rDiff")
#};rm(i)
#glmm0 = glmm0[!is.na(glmm0[,1]),]
#yr2 = lm(glmm0[,1]~glmm0[,4]+glmm0[,5]+glmm0[,6]+glmm0[,7])
#yr2.1 = lm(glmm0[,1]~glmm0[,6])

##### correlation test #####
#c1 = cor.test(glmm0[,4],glmm0[,1], method="spearman") # non-sig; gLV proportion vs CFTRm
#c2 = cor.test(glmm0[,5],glmm0[,1], method="spearman") # non-sig; gLV proportion vs antimicrobials
#c3 = cor.test(glmm0[,6],glmm0[,1], method="spearman") # alm-sig; gLV proportion vs chemicals
#c4 = cor.test(glmm0[,7],glmm0[,1], method="spearman") # non-sig; gLV proportion vs interactions

#cAll = data.frame(fac=colnames(tReat)[-c(1:2)],rho=c(c1$estimate, c2$estimate, c3$estimate, c4$estimate), S=c(c1$statistic, c2$statistic, c3$statistic, c4$statistic), adj.p=c(c1$p.value, c2$p.value, c3$p.value, c4$p.value)*4)
