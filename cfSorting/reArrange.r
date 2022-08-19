#!/bin/env Rscript
# author: ph-u
# script: reArrange.r
# desc: rearrange columns to medical,microbe dataframes
# in: Rscript reArrange.r 1> ../data/reArrangeRec.txt
# out: ../data/cf425MedMic.rda, ../data/{otherSp,medName}-qcREF.csv, ../data/{cf425Medic,cf425Micro,cf425_noInfoList}.csv
# arg: 0
# date: 20220315,20220404,20220802

##### env #####
cat("read reference data:",date(),"\n")
library(stringr);
load("../data/cf425FULL.rda");
sP = read.csv("../data/otherSp-qcList.csv", header=T, stringsAsFactors=F);
gC = read.csv("../data/genCols-qcList.csv", header=T, stringsAsFactors=F);
mD = read.csv("../data/medName-qcList.csv", header=T, stringsAsFactors=F);
mdNA = mD[which(is.na(mD$QC.standardization) | is.na(mD$func.component)),1]
mD = mD[which(!is.na(mD$QC.standardization) | !is.na(mD$func.component)),]
miCol0 = c("s05culturespeciesresistotherspec","cult_specify","s05culturespeciesfungalotherspec","s05culturespeciesviralotherspeci")
miCol = which(colnames(indRef) %in% miCol0) # columns of microbial records

##### construct medical,microbe dataframes #####
cat("construct result dataframes:",date(),"\n")
m0 = str_split(c(gC$QC.standardization[which(gC$med.microbe.other=="medical")],mD$QC.standardization,mD$func.component),";")
namMedi = c();for(i in 1:length(m0)){namMedi = c(namMedi, m0[[i]])};namMedi = unique(trimws(tolower(namMedi), which="both"));
namMedi = namMedi[order(namMedi)];namMedi = namMedi[!is.na(namMedi)]
mEdic = as.data.frame(matrix(0, nr=nrow(rEc), nc=2+length(namMedi)))
colnames(mEdic) = c(colnames(rEc)[1:2], namMedi)

namMicro = unique(c(as.character(read.table(text=unique(gC$QC.standardization[which(gC$med.microbe.other=="microbe")]),sep=";")[,1]), as.character(read.table(text=unique(sP$QC.standardization),sep=";")[,1])))
namMicro = namMicro[order(namMicro)]
mIcro = as.data.frame(matrix(0, nr=nrow(rEc), nc=2+length(namMicro)))
colnames(mIcro) = c(colnames(rEc)[1:2], namMicro)

mEdic[,1:2] = mIcro[,1:2] = rEc[,1:2]

##### microbe reference dataframe ##### 20220323,20220327,20220814
cat("microbe reference dataframe:",date(),"\n")
eLim = c();for(i in 1:nrow(sP)){if((i %in% eLim)==F){ # eliminate repeated longer text
	a = setdiff(grep(paste0("\\Q",sP[i,1]),sP[,1]),i)
	eLim = c(eLim,a[which(sP$QC.standardization[a]==sP$QC.standardization[i] & sP$note[a]==sP$note[i])])
}};sP = sP[-unique(eLim),]

qcREF = as.data.frame(matrix(0,nr=nrow(sP),nc=length(namMicro)+1))
colnames(qcREF) = c(colnames(sP)[1],namMicro);qcREF[,1] = sP[,1]
a = str_split(sP$QC.standardization,";")
for(i in 1:length(a)){
	qcREF[i,a[[i]]] = ifelse(sP$note[i]==-1,-1,1)
	a0 = setdiff(grep(paste0("\\Q",qcREF[i,1]),qcREF[,1]),which(qcREF[,1]==qcREF[i,1])) # tag no extra info longer text
	if(length(a0)>0){qcREF[a0,a[[i]]] = qcREF[a0,a[[i]]]-ifelse(qcREF[a0,a[[i]]]>=0,str_count(qcREF[a0,1],paste0("\\Q",qcREF[i,1])),0)}
};qcREF = qcREF[-which(rowSums(abs(qcREF[,-1]))==0),] # eliminate longer text with no extra info

write.csv(qcREF,"../data/otherSp-qcREF.csv", quote=F, row.names=F)

##### map medical reference dataframe ##### 20220802
cat("medication reference dataframe:",date(),"\n")
mdNA = c(" mc",mdNA," 2units before breakfast and lunch sc", "chronic", "half ", " od ", "rescue", "unknown", " and ", "slow", "normal", "adult", "liquid", "alt days", "oral ", "as required", "asdirected", "aviva strips", "bg star strips", "cassette", "complan", "course complete", "dario", "mart up to", "microscoop", "not known", "l capsules", "ml eye drops", "po 30-60ml od", "for ", "day", "freestyle", "test strip", "testing strip", "butter", "fast clix", "fastclix", "bd", "may", "week", "other", "thick", "easy", "testing strips", "discomfort", "from feb to april", "cold", "for enuresis", "accu chek performa", "feed", "well teen tablets", "during hayfever season", "wk", "wfi", "pollen allergy", "lp", "alt months", "l", "r", "for 3 months", "nil", "wexham prescribed", "release", "2.5ml bd", "1.5 sachets per day", "not taking", "no current medications", "nil", "fluid thickening agent", "blood glucose stripes", "bone health", "blinded", "as req", "1.5 sachets per day", "2.5ml bd", "od", "dose being titrated up", "minute", "null", "rarely needed", "retired medication", "smart")
mdRef = as.data.frame(matrix(0,nr=nrow(mD)+length(mdNA),nc=length(namMedi)+1))
colnames(mdRef) = c(colnames(mD)[1],namMedi)
mdRef[,1] = c(mD[,1],mdNA)
for(i in 1:nrow(mdRef)){if(i<=nrow(mD)){
	i0 = c(trimws(str_split(mD$QC.standardization[i],";")[[1]], which="both"),str_split(mD$func.component[i],";")[[1]])
	mdRef[i,unique(tolower(i0))] = ifelse(mD$note[i]==-1,-1,1)
}}
cat("medication reference data export:",date(),"\n")
write.csv(mdRef,"../data/medName-qcREF.csv", quote=F, row.names=F)

##### map medical indirect record ##### 20220324,20220802,20220810
i0 = setdiff(3:ncol(indRef),miCol) # vitamin record inaccurate
m0 = indRef[,i0[1]]
for(i in i0[-1]){m0 = paste(m0,indRef[,i],sep=";")}
m0 = tolower(gsub(";","@",gsub("\"","",gsub(";NA$","",gsub("NA;","",m0))))) #tolower(gsub("\'","!",m0))
cat("map medical indirect:",date(),"\n")
##### map each data row 20220813
m1 = str_split(m0,"@")
for(i in 1:length(m1)){
	m2 = which(mdRef[,1] %in% m1[[i]]) # name match cols
	if(length(grep("a&d",gsub(" ","",gsub(",","&",m1[[i]]))))>0){m2 = c(m2,which(mdRef[,1] %in% "a&d"))} # non-matching exception
	if(i%%1000==0){cat(i,"(",round(i/length(m1)*100,2),"% )",date(),"\n")}
	m3 = setdiff(which(colSums(mdRef[m2,-1])>0),1) # count items
	if(length(m3)>0){mEdic[i,m3+2] = 1}
}
cat("convert medical indirect into presence/absence:",date(),"\n")
for(i in 3:ncol(mEdic)){mEdic[,i] = ifelse(mEdic[,i]>0,1,0)};rm(i)

##### map microbial indirect record ##### 20220324
cat("standardise microbial indirect source columns:",date(),"\n")
m0 = indRef[,miCol[1]]
for(i in miCol[-1]){m0 = paste(m0,indRef[,i],sep="@")};rm(i)
for(i in c("@NA@","NA@","@NA")){m0 = gsub(i,"@",m0)}
for(i in c("[.]","\\s+")){m0 = gsub(i," ",m0)}
m0 = tolower(m0)
cat("map microbial indirect:",date(),"\n")
for(i in 1:nrow(qcREF)){
	a = grep(paste0("\\Q",qcREF[i,1]),m0);#if(length(a)<1){break}
	a0 = setdiff(which(qcREF[i,]!=0),1)+1
	mIcro[a,a0] = mIcro[a,a0] + qcREF[rep(i,length(a)),a0-1]*str_count(m0[a],paste0("\\Q",qcREF[i,1]))
	cat(i,"(",round(i/nrow(qcREF)*100,2),"% ): \"",qcREF[i,1],"\" ( n =",length(a),")",date(),"\n")
}
#	wCol = which(qcREF[i,-1]!=0)+2 # column id for mIcro
#	a0=0;for(i0 in miCol){
#		a = grep(qcTag[i],indRef[,i0]);a0 = a0 + length(a)
#		if(length(a)>0){
#			if(length(wCol)>1){x = qcREF[i,wCol-1][rep(seq_len(1),each=length(a)),]}else{x = qcREF[i,wCol-1]}
#			mIcro[a,wCol] = mIcro[a,wCol] + x*str_count(indRef[a,i0],qcTag[i])
#	}};cat(i,":",date(),"(n =",a0,"-",qcREF$input[i],")\n")
#};rm(i,wCol,a0,a,x)
cat("convert microbial indirect into presence/absence:",date(),"\n")
for(i in 3:ncol(mIcro)){mIcro[,i] = ifelse(mIcro[,i]>0,1,0)};rm(i)

##### map direct record ##### 20220324,20220804
cat("map direct record (medical):",date(),"\n")
mD0 = gC[which(gC$med.microbe.other=="medical"),]
for(i in 1:nrow(mD0)){
        a = which(colnames(dirRef)==mD0$input[i])
        a0 = str_split(tolower(mD0$QC.standardization[i]),";")[[1]]
	i1 = c();for(i0 in 1:length(a0)){i1 = c(i1,which(colnames(mEdic)==a0[i0]))} # exact name
	mEdic[,i1] = mEdic[,i1] + as.numeric(dirRef[,a])
}

cat("map direct record (microbial):",date(),"\n")
#mC0 = gC[which(gC$med.microbe.other=="microbe"),]
#for(i in 1:nrow(mC0)){
#	a = which(colnames(dirRef)==mC0$input[i])
#	a0 = str_split(mC0$QC.standardization[i],";")[[1]]
#	i1 = c();for(i0 in 1:length(a0)){i1 = c(i1,which(colnames(mIcro)==a0[i0]))} # exact name
#	mIcro[,i1] = mIcro[,i1] + as.numeric(dirRef[,a])
#}
for(i in 3:ncol(dirRef)){
	tK0 = which(gC$input==trimws(colnames(dirRef)[i], which="both"))
	tK = gC$med.microbe.other[tK0]
	if(tK=="microbe"){
		a0 = as.character(read.table(text=gC$QC.standardization[tK0], sep=";"))
		for(i0 in 1:length(a0)){
#			if(tK=="medical"){
#				a = which(colnames(mEdic)==tolower(a0[i0]))
#				mEdic[,a] = mEdic[,a] + as.numeric(dirRef[,i])
#			}else if(tK=="microbe"){
				a = which(colnames(mIcro)==a0[i0])
				mIcro[,a] = mIcro[,a] + as.numeric(dirRef[,i])
}}};rm(i,i0,a,a0,tK,tK0)

##### double check each person has one entry at each data available year ##### 20220816
cat("check misrecognized entries:",date(),"\n")

## eliminate non-informative data -- no data input on both medication & microbes
eLim = which(rowSums(mIcro[,-c(1:2)])==0) # which(rowSums(mIcro[,-c(1:2)])==0 & rowSums(mEdic[,-c(1:2)])==0)
write.csv(mIcro[eLim,1:2], "../data/cf425_noInfoList.csv", quote=F, row.names=F)
mEdic = mEdic[-eLim,];mIcro = mIcro[-eLim,]

## fuse repeated patients occurrance in the same year
mD00 = as.data.frame(matrix(,nr=0,nc=ncol(mEdic))); colnames(mD00) = as.character(colnames(mEdic))
mC00 = as.data.frame(matrix(,nr=0,nc=ncol(mIcro))); colnames(mC00) = as.character(colnames(mIcro))
eLim = c();nR = nrow(mEdic);for(i in 1:nR){ if(!(i %in% eLim)){
	i0 = which(mEdic[,1]==mEdic[i,1] & mEdic[,2]==mEdic[i,2]);
	if(i%%1000==0){cat(round(i/nR*100,2),"% scanned on repeated patients occurrance:",date(),"\n")}
	if(length(i0)>1){
		cat(i,":",paste(mEdic[i,1:2], collapse=" "),"n =",length(i0),"-",date(),"\n");
		mD00[nrow(mD00)+1,] = c(mEdic[i,1:2],ifelse(colSums(mEdic[i0,-c(1:2)])>0,1,0));
		mC00[nrow(mC00)+1,] = c(mIcro[i,1:2],ifelse(colSums(mIcro[i0,-c(1:2)])>0,1,0));
		eLim = c(eLim,i0);
}}};mEdic = rbind(mEdic[-eLim,],mD00);mIcro = rbind(mIcro[-eLim,],mC00);

##### export #####
cat("record notation final conversion and export:",date(),"\n")
for(i in 3:ncol(mEdic)){mEdic[,i] = ifelse(mEdic[,i]>0,1,0)};rm(i)
write.csv(mEdic,"../data/cf425Medic.csv", quote=F, row.names=F)
for(i in 3:ncol(mIcro)){mIcro[,i] = ifelse(mIcro[,i]>0,1,0)};rm(i)
write.csv(mIcro,"../data/cf425Micro.csv", quote=F, row.names=F)
save(mEdic,mIcro, file="../data/cf425MedMic.rda")
cat("data sorting completed:",date(),"\n")
