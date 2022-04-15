#!/bin/env Rscript
# author: ph-u
# script: reArrange.r
# desc: rearrange columns to medical,microbe dataframes
# in: Rscript reArrange.r 1> ../data/reArrangeRec.txt
# out: ../data/cf425MedMic.rda, ../data/{otherSp-qcREF,cf425Medic,cf425Micro}.csv
# arg: 0
# date: 20220315,20220404

##### env #####
cat("read reference data:",date(),"\n")
library(stringr);
load("../data/cf425FULL.rda");
sP = read.csv("../data/otherSp-qcList.csv", header=T, stringsAsFactors=F);
gC = read.table("../data/genCols-qcList.tsv", sep="\t", header=T, stringsAsFactors=F);

##### construct medical,microbe dataframes #####
cat("construct result dataframes:",date(),"\n")
namMedi = unique(as.character(read.table(text=unique(gC$QC.standardization[which(gC$med.microbe.other=="medical")]),sep=";")[,1]))
namMedi = namMedi[order(namMedi)]
mEdic = as.data.frame(matrix(0, nr=nrow(rEc), nc=2+length(namMedi)))
colnames(mEdic) = c(colnames(rEc)[1:2], namMedi)

namMicro = unique(c(as.character(read.table(text=unique(gC$QC.standardization[which(gC$med.microbe.other=="microbe")]),sep=";")[,1]), as.character(read.table(text=unique(sP$QC.standardization),sep=";")[,1])))
namMicro = namMicro[order(namMicro)]
mIcro = as.data.frame(matrix(0, nr=nrow(rEc), nc=2+length(namMicro)))
colnames(mIcro) = c(colnames(rEc)[1:2], namMicro)

mEdic[,1:2] = mIcro[,1:2] = rEc[,1:2]

##### microbe reference dataframe ##### 20220323,20220327
cat("microbe reference dataframe:",date(),"\n")
qcREF = as.data.frame(matrix(0,nr=nrow(sP),nc=length(namMicro)+1))
colnames(qcREF) = c(colnames(sP)[1],namMicro);qcREF[,1] = sP[,1]
nEg = which(sP$note==-1)
for(i in 1:nrow(qcREF)){
	a = as.character(read.table(text=sP$QC.standardization[i],sep=";")[1,])
	qcREF[i,a] = qcREF[i,a] + ifelse(i %in% nEg,-1,1)
};rm(i,a)
qcTag = gsub("[$]","@",gsub("[(]","@",gsub("[)]","@",gsub("[+]","@",qcREF[,1]))))

## set text substrings to -ve 20220404
eLim=c();for(i in 1:nrow(qcREF)){if((i %in% eLim)==F){
	qC = setdiff(grep(qcTag[i],qcTag),eLim);qC = qC[qC!=i]
	if(length(qC)>0){ qC0 = which(qcREF[i,-1]>0)#; cat(i,",")
		for(i0 in qC){
			if(all(qcREF[i,-1]==qcREF[i0,-1])){
				eLim=c(eLim,i0)#; cat(i,",",i0,";")
			}else{
				qcREF[i0,qC0+1] = qcREF[i0,qC0+1]-1
}}}}};rm(i,qC)#;cat("\n")
qcREF = qcREF[-eLim,] # eliminate skippable textstring pattern

## limit text pattern signal range into [-1,+1]
x0 = colnames(qcREF)[1]
qcREFt = qcREF[,-1]
qcREFt[qcREFt<0] = -1
qcREF = cbind(qcREF[,1],qcREFt)
colnames(qcREF)[1] = x0;rm(x0,qcREFt)
write.csv(qcREF,"../data/otherSp-qcREF.csv", quote=F, row.names=F)

##### map direct record ##### 20220324
cat("map direct record:",date(),"\n")
for(i in 3:ncol(dirRef)){
	tK0 = which(gC$input==trimws(colnames(dirRef)[i], which="both"))
	tK = gC$med.microbe.other[tK0]
	if(tK!="other"){
		a0 = as.character(read.table(text=gC$QC.standardization[tK0], sep=";"))
		for(i0 in 1:length(a0)){
			if(tK=="medical"){
				a = which(colnames(mEdic)==a0[i0])
				mEdic[,a] = mEdic[,a] + as.numeric(dirRef[,i])
			}else if(tK=="microbe"){
				a = which(colnames(mIcro)==a0[i0])
				mIcro[,a] = mIcro[,a] + as.numeric(dirRef[,i])
}}}};rm(i,i0,a,a0,tK,tK0)

##### map medical indirect record ##### 20220324
cat("map medical indirect:",date(),"\n")
mediNum = which(gC$med.microbe.other=="medical")
for(i0 in 1:length(mediNum)){for(i1 in 3:ncol(indRef)){
	a = grep(gC$input[mediNum[i0]],tolower(indRef[,i1]))
	if(length(a)>0){
		a0 = which(colnames(mEdic) %in% as.character(read.table(text=gC$QC.standardization[mediNum[i0]], sep=";")))
		mEdic[a,a0] = mEdic[a,a0] +1
}}};rm(mediNum,i0,i1,a,a0)
for(i in 3:ncol(mEdic)){mEdic[,i] = ifelse(mEdic[,i]>0,1,0)};rm(i)
write.csv(mEdic,"../data/cf425Medic.csv", quote=F, row.names=F)

##### map microbial indirect record ##### 20220324
cat("standardise microbial indirect source columns:",date(),"\n")
for(i in 3:ncol(indRef)){
	indRef[,i] = gsub("[$]","@",gsub("[(]","@",gsub("[)]","@",gsub("[+]","@", gsub("\\s+", " ", gsub("[.*]"," ",indRef[,i]))))))
};rm(i)
cat("map microbial indirect:",date(),"\n")
qcTag = gsub("[$]","@",gsub("[(]","@",gsub("[)]","@",gsub("[+]","@",qcREF[,1])))) # filtered row
for(i in 1:nrow(qcREF)){
	wCol = which(qcREF[i,-1]!=0)+2 # column id for mIcro
	a0=0;for(i0 in 3:ncol(indRef)){
		a = grep(qcTag[i],indRef[,i0]);a0 = a0 + length(a)
		if(length(a)>0){
			if(length(wCol)>1){x = qcREF[i,wCol-1][rep(seq_len(1),each=length(a)),]}else{x = qcREF[i,wCol-1]}
			mIcro[a,wCol] = mIcro[a,wCol] + x*str_count(indRef[a,i0],qcTag[i])
	}};cat(i,":",date(),"(n =",a0,"-",qcREF$input[i],")\n")
};rm(i,wCol,a0,a,x)
cat("convert microbial indirect into presence/absence:",date(),"\n")
for(i in 3:ncol(mIcro)){mIcro[,i] = ifelse(mIcro[,i]>0,1,0)};rm(i)
write.csv(mIcro,"../data/cf425Micro.csv", quote=F, row.names=F)

save(mEdic,mIcro, file="../data/cf425MedMic.rda")
cat("data sorting completed:",date(),"\n")
#for(i in 2:ncol(qcREF)){if(all(qcREF[,i]<=0)){cat(i,",",colnames(qcREF)[i],"\n")}};rm(i)
#for(i in 3:ncol(mIcro)){if(all(mIcro[,i]<=0)){cat(i,",",colnames(mIcro)[i],"\n")}};rm(i)
#for(i in 3:ncol(dirRef)){if(all(dirRef[,i]<=0)){cat(i,",",colnames(dirRef)[i],"\n")}};rm(i)
