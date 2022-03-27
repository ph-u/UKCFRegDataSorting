#!/bin/env Rscript
# author: ph-u
# script: reArrange.r
# desc: rearrange columns to medical,microbe dataframes
# in: Rscript reArrange.r
# out: ../data/cf425MedMic.rda, ../data/{otherSp-qcREF,cf425Medic,cf425Micro}.csv
# arg: 0
# date: 20220315

##### env #####
cat("read reference data::",date(),"\n")
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
zZ = gsub("[$]","@",gsub("[(]","@",gsub("[)]","@",gsub("[+]","@",qcREF[,1]))))
for(i in 1:length(zZ)){if(any(qcREF[i,-1]!=0)){
	z9 = grep(zZ[i],zZ)
	z9 = z9[which(z9!=i)]
	if(sP$note[i]!=-1){z9 = z9[which(sP$note[z9]!=-1)]}else{z9 = z9[which(sP$note[z9]==-1)]} # only matching categories
## extract entries of identical categories (minimise duplicated id)
	z0 = z9[which(sP$QC.standardization[i]==sP$QC.standardization[z9])]
## preserve entries with non-identical categories
	z1 = setdiff(z9,z0)
	iNclude = grep(sP$QC.standardization[i],sP$QC.standardization[z1])
	if(length(iNclude)>0){z1 = z1[-iNclude]} # not mess up hybrids
## eliminate regex text grep artifact
	z2 = c(z0,z1)
	if(length(z2)>0){ # textstring is a subset of other textstrings
                for(i0 in 1:length(z2)){
			if(z2[i0] %in% z0){tX="<-"}else{tX="XX"}
                        cat(i,tX,z2[i0],":",zZ[i],tX,zZ[z2[i0]],"\n")
			i2 = which(qcREF[i,]!=0); i2 = i2[i2!=1]
			for(i3 in i2){
				qcREF[z2[i0],i3] = max(-1, qcREF[z2[i0],i3]-abs(qcREF[i,i3]))
}}}}};rm(i)
## remove replaceable textstrings
z0 = c();for(i in 1:nrow(qcREF)){if(sP$note[i]!=-1 & any(as.numeric(qcREF[i,-1])>0)==F){z0 = c(z0,i)}};rm(i)
qcREF = qcREF[-z0,]
## set text substrings to -1
zZ0 = gsub("[$]","@",gsub("[(]","@",gsub("[)]","@",gsub("[+]","@",qcREF[,1]))))
exHybrid = setdiff(zZ0,zZ[grep(";",sP$QC.standardization)])
for(i in 1:nrow(qcREF)){if(zZ0[i] %in% exHybrid){
	i0 = which(zZ0 %in% exHybrid[grep(zZ0[i],exHybrid)])
	if(length(i0)>1){i0 = i0[i0!=i]
		for(i1 in 1:length(i0)){
			i2 = which(qcREF[i,]!=0); i2 = i2[i2!=1]
			for(i3 in i2){
				qcREF[i0[i1],i3] = max(-1, qcREF[i0[i1],i3]-abs(qcREF[i,i3]))
}}}}}
write.csv(qcREF,"../data/otherSp-qcREF.csv", quote=F, row.names=F)

##### map direct record ##### 20220324
cat("map direct record:",date(),"\n")
for(i in 3:ncol(dirRef)){
	tK = gC$med.microbe.other[which(gC$input==colnames(dirRef)[i])]
	if(length(tK>0)){ # some colnames dropped due to duplicated info
		if(tK!="other"){
			a0 = gC$QC.standardization[which(gC$input==colnames(dirRef)[i])]
		};if(tK=="medical"){
			a = which(colnames(mEdic)==a0)
			mEdic[,a] = mEdic[,a] + as.numeric(dirRef[,i])
		}else if(tK=="microbe"){
			a = which(colnames(mIcro)==a0)
			mIcro[,a] = mIcro[,a] + as.numeric(dirRef[,i])
}}};rm(i,a,a0)

##### map medical indirect record ##### 20220324
cat("map medical indirect:",date(),"\n")
mediNum = which(gC$med.microbe.other=="medical")
for(i0 in 1:length(mediNum)){for(i1 in 3:ncol(indRef)){
	a = grep(gC$input[mediNum[i0]],tolower(indRef[,i1]))
	if(length(a)>0){
		a0 = which(colnames(mEdic)==gC$QC.standardization[mediNum[i0]])
		mEdic[a,a0] = mEdic[a,a0] +1
}}};rm(mediNum,i0,i1,a,a0)
for(i in 3:ncol(mEdic)){mEdic[,i] = ifelse(mEdic[,i]>0,1,0)};rm(i)
write.csv(mEdic,"../data/cf425Medic.csv", quote=F, row.names=F)

##### map microbial indirect record ##### 20220324
cat("standardise microbial indirect source columns:",date(),"\n")
for(i in 3:ncol(indRef)){
	indRef[,i] = gsub("[$]","@",gsub("[(]","@",gsub("[)]","@",gsub("[+]","@", gsub("\\s+", " ", gsub("[.*]"," ",indRef[,i]))))))
};rm(i)

cat("map microbial indirect:",date(),"\n") # zZ0 as input index
for(i in 1:nrow(qcREF)){
	wCol = which(qcREF[i,-1]!=0)+2 # column id for mIcro
	a = c();for(i0 in 3:ncol(indRef)){a = c(a,grep(zZ0[i],indRef[,i0]))};a = unique(a)
	cat(i,":",date(),"(n =",length(a),"-",qcREF$input[i],")\n")
	for(i0 in a){
		mIcro[i0,wCol] = mIcro[i0,wCol] + qcREF[i,wCol-1]
}}
cat("convert microbial indirect into presence/absence:",date(),"\n")
for(i in 3:ncol(mIcro)){mIcro[,i] = ifelse(mIcro[,i]>0,1,0)};rm(i)
write.csv(mIcro,"../data/cf425Micro.csv", quote=F, row.names=F)

save(mEdic,mIcro, file="../data/cf425MedMic.rda")
cat("data sorting completed:",date(),"\n")
