#!/bin/env Rscript
# author: ph-u
# script: mutSort.r
# desc: sort CF mutation of pwCF
# in: Rscript mutSort.r {0,1}
# out: ../data/mutSort.{csv,txt} OR ../data/mutPWCF.csv
# arg: 0
# date: 20221126

##### import #####
argv=(commandArgs(T))
library(xlsx2dfs)
cat("import start:",date(),"\n")
mT = xlsx2dfs("../p_cf425A/425A_DataRequest_final.xlsx", rowNames=F, detectDates = T)
mtLib = xlsx2dfs("../p_cf425A/mutation_lookup_table_April2022.xlsx", rowNames=F, detectDates = T)
if(argv[1]!=0){
	mtRef = read.csv("../data/mutLib.csv", header=T, stringsAsFactors=F)
	mtRef[,2] = read.csv("../data/mutLib_r.csv", header=T, stringsAsFactors=F)[,2]
}

##### tabulize mutation data #####
cat("fusion start:",date(),"\n")
mT0 = mT[[2]];mT0$year = names(mT)[2]
for(i in 3:length(mT)){
	a0 = mT[[i]];a0$year = names(mT)[i]
	mT0 = merge(mT0, a0, all=T)

};rm(i,a0)
mT0$s02genotypeqty = NULL # this column irrelevant (recording quantity of mutation only)

##### capture latest year mutation record #####
cat("capture start:",date(),"\n")
pwCF = unique(mT0$regid_anon)
mT1 = as.data.frame(matrix(nr=length(pwCF),nc=ncol(mT0)))
colnames(mT1) = colnames(mT0)
mT1$regid_anon = pwCF
for(i in 1:nrow(mT1)){
	a0 = mT0[which(mT0$regid_anon==mT1$regid_anon[i]),]
	mT1[i,] = a0[which(a0$year==max(a0$year)),]
};rm(pwCF,i,a0)
for(i in 3:ncol(mT1)){mT1[,i] = gsub(",",";",mT1[,i])};rm(i)

##### match mutation notations #####
cat("matching start:",date(),"-")
if(argv[1]==0){
	mL0 = mtLib[[1]][,2:3]
	mSg = c();for(i in 3:ncol(mT1)){ cat(i,",");for(i0 in which(!is.na(mT1[,i]))){ if(
		is.na(suppressWarnings(as.numeric(mT1[i0,i])))==F & # record is a pure number
		suppressWarnings(as.numeric(mT1[i0,i]))<=max(as.numeric(mL0[,1])) & # number in lib range
		suppressWarnings(as.numeric(mT1[i0,i]))>=min(as.numeric(mL0[,1]))){ # number in lib range
			mSg = c(mSg,paste(mT1[i0,i],":",mL0[which(mL0[,1]==round(as.numeric(mT1[i0,i]),1)),2]))
			mT1[i0,i] = mL0[which(mL0[,1]==round(as.numeric(mT1[i0,i]),1)),2]
	}}};rm(i);cat("\n")

	i0=c();for(i in 3:ncol(mT1)){i0 = c(i0,mT1[,i])}
	write.table(unique(mSg),"../data/mutLib.txt", row.names=F, col.names=F,sep="\t", quote=F)
	write.csv(data.frame(mutation=gsub(",",";",unique(i0)[order(unique(i0))]),standard=""), "../data/mutLib.csv", row.names=F, quote=F)
}else{
## standardize numeric expressions for mutation types
	mT2 = as.data.frame(matrix(nr=nrow(mT1),nc=ncol(mT1)))
	colnames(mT2) = colnames(mT1)
	mT2[,1:2] = mT1[,1:2]
	cat("mapping ")
	for(i in 3:ncol(mT2)){ cat(i,",");for(i0 in which(!is.na(mT1[,i]))){ if(
		length(which(mtRef$mutation==mT1[i0,i])) == 1
	){
		mT2[i0,i] = mtRef$standard[which(mtRef$mutation==mT1[i0,i])]
	}else{
		mT2[i0,i] = mT1[i0,i]
	}}};rm(i,i0);cat("\n")

## convert mutation numeric expressions to binary presence/absence reference table
	cat("identifying pwCF mutation types:",date(),"\n")
	muTable = as.data.frame(matrix(0,nr=nrow(mT1),nc=nrow(mtLib[[1]])+3))
	colnames(muTable) = c(colnames(mT1)[1:2],gsub("==","=",gsub("NA","",gsub(",",";",paste(mtLib[[1]]$legacyname,mtLib[[1]]$proteinname,mtLib[[1]]$cdnaname, sep="=")))),"other")
	muTable[,1:2] = mT1[,1:2]
	for(i in 1:nrow(muTable)){
		i0 = strsplit(as.character(mT2[i,-c(1:2)]), ";")
		i2 = c();for(i1 in 1:length(i0)){if(!is.na(i0[[i1]][1])){i2 = c(i2,i0[[i1]])}}
		i0 = which(mtLib[[1]]$mutationcode %in% unique(as.numeric(gsub("other",999,tolower(i2))))) # get matching mutation identifier codes
		if(length(grep("other",tolower(i2)))>0){i0 = c(i0,ncol(muTable)-2)} # detect "other"
		muTable[i, 2+i0] = 1 # presence of a mutation category
	};rm(i,i0,i1,i2)

	cat("accounting for multi-types mutation categories:",date(),"\n")
	a = strsplit(colnames(muTable), "=")
	for(i in 3:ncol(muTable)){
		a0 = c();for(i0 in 1:length(a[[i]])){a0 = c(a0, grep(paste0("\\Q",a[[i]][i0]),colnames(muTable)))}
		if(length(unique(a0))>1){muTable[,i] = as.numeric(apply(muTable[,unique(a0)],1,sum) > 0)} # single mutation now contains multi-mutation type data
	};rm(a,i,i0,a0)

	write.csv(muTable,"../data/mutPWCF.csv", row.names=F, quote=F)
}
