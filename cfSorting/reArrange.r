#!/bin/env Rscript
# author: ph-u
# script: reArrange.r
# desc: rearrange columns to medical,microbe dataframes
# in: Rscript reArrange.r 1> ../data/reArrangeRec.txt
# out: ../data/cf425MedMic.rda, ../data/{otherSp-qcREF,cf425Medic,cf425Micro}.csv
# arg: 0
# date: 20220315,20220404,20220802

##### env #####
cat("read reference data:",date(),"\n")
library(stringr);
load("../data/cf425FULL.rda");
sP = read.csv("../data/otherSp-qcList.csv", header=T, stringsAsFactors=F);
gC = read.csv("../data/genCols-qcList.csv", header=T, stringsAsFactors=F);
mD = read.csv("../data/medName-qcList.csv", header=T, stringsAsFactors=F);
mD = mD[which(!is.na(mD$QC.standardization) | !is.na(mD$func.component)),]
miCol = c(3,7,20,21) # columns of microbial records

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

##### map medical reference dataframe ##### 20220802
cat("medication reference dataframe:",date(),"\n")
mdRef = as.data.frame(matrix(0,nr=nrow(mD),nc=length(namMedi)+1))
colnames(mdRef) = c(colnames(mD)[1],namMedi)
mdRef[,1] = mD[,1]
for(i in 1:nrow(mdRef)){
	i0 = c(trimws(str_split(mD$QC.standardization[i],";")[[1]], which="both"),str_split(mD$func.component[i],";")[[1]])
	mdRef[i,unique(tolower(i0))] = ifelse(mD$note[i]==-1,-1,1)
}
## at least 1 med similar name diff function
for(i in 1:nrow(mdRef)){i0 = setdiff(grep(paste0("\\Q",mdRef[i,1]),mdRef[,1]),i)
	if(length(i0)>0){cat(i,"included in",length(i0),"entries:",mdRef[i,1],"<--",paste(mdRef[i0,1],collapse=","),";\n")
		for(i1 in 1:length(i0)){
			i2 = setdiff(which(mdRef[i0[i1],]!=0),1)
			mdRef[i0[i1],i2] = mdRef[i0[i1],i2] - abs(mdRef[i,i2]) ## secure longer texts info if a shorter text is available
	}}
};mdRef = mdRef[which(rowSums(mdRef[,-1])!=0),] # eliminate physical lines, retain keyphrases consideration in categorisation step (ln2963 -> 2088)
write.csv(mdRef,"../data/medName-qcREF.csv", quote=F, row.names=F)

##### map medical indirect record ##### 20220324,20220802
cat("map medical indirect:",date(),"\n")
i0 = setdiff(3:ncol(indRef),miCol)
m0 = indRef[,i0[1]]
for(i in i0[-1]){m0 = paste(m0,indRef[,i],sep=";")}
m0 = tolower(gsub("\"","",m0)) #tolower(gsub("\'","!",m0))
for(i in 1:nrow(mdRef)){
	i1 = grep(paste0("\\Q",mdRef[i,1]),m0) # grep literal string
	i2 = which(mdRef[i,]!=0)
	cat(i,"(",round(i/nrow(mdRef)*100,2),"% ): n =",length(i1),"(",mdRef[i,1],")",date(),"\n")
	if(length(i1)==0){break}
	mEdic[i1,i2[-1]+1] = mdRef[i,i2[-1]]
}
cat("convert medical indirect into presence/absence:",date(),"\n")
for(i in 3:ncol(mEdic)){mEdic[,i] = ifelse(mEdic[,i]>0,1,0)};rm(i)

##### map microbial indirect record ##### 20220324
cat("standardise microbial indirect source columns:",date(),"\n")
for(i in miCol){
	indRef[,i] = tolower(gsub("[$]","@",gsub("[(]","@",gsub("[)]","@",gsub("[+]","@", gsub("\\s+", " ", gsub("[.*]"," ",indRef[,i])))))))
};rm(i)
cat("map microbial indirect:",date(),"\n")
qcTag = gsub("[$]","@",gsub("[(]","@",gsub("[)]","@",gsub("[+]","@",qcREF[,1])))) # filtered row
for(i in 1:nrow(qcREF)){
	wCol = which(qcREF[i,-1]!=0)+2 # column id for mIcro
	a0=0;for(i0 in miCol){
		a = grep(qcTag[i],indRef[,i0]);a0 = a0 + length(a)
		if(length(a)>0){
			if(length(wCol)>1){x = qcREF[i,wCol-1][rep(seq_len(1),each=length(a)),]}else{x = qcREF[i,wCol-1]}
			mIcro[a,wCol] = mIcro[a,wCol] + x*str_count(indRef[a,i0],qcTag[i])
	}};cat(i,":",date(),"(n =",a0,"-",qcREF$input[i],")\n")
};rm(i,wCol,a0,a,x)
cat("convert microbial indirect into presence/absence:",date(),"\n")
for(i in 3:ncol(mIcro)){mIcro[,i] = ifelse(mIcro[,i]>0,1,0)};rm(i)

##### map direct record ##### 20220324,20220804
cat("map direct record (medical):",date(),"\n")
mD0 = gC[which(gC$med.microbe.other=="medical"),]
for(i in 1:nrow(mD0)){
        a = which(colnames(dirRef)==mD0$input[i])
        a0 = str_split(tolower(mD0$QC.standardization[i]),";")[[1]]
#        i1 = c();for(i0 in 1:length(a0)){i1 = c(i1,grep(tolower(a0[i0]),colnames(mEdic)))} # fuzzy name identification
	i1 = c();for(i0 in 1:length(a0)){i1 = c(i1,which(colnames(mEdic)==a0[i0]))} # exact name
	mEdic[,i1] = mEdic[,i1] + as.numeric(dirRef[,a])
}

cat("map direct record (microbial):",date(),"\n")
mC0 = gC[which(gC$med.microbe.other=="microbe"),]
for(i in 1:nrow(mC0)){
	a = which(colnames(dirRef)==mC0$input[i])
	a0 = str_split(mC0$QC.standardization[i],";")[[1]]
	i1 = c();for(i0 in 1:length(a0)){i1 = c(i1,which(colnames(mIcro)==a0[i0]))} # exact name
	mIcro[,i1] = mIcro[,i1] + as.numeric(dirRef[,a])
}
#for(i in 3:ncol(dirRef)){
#	tK0 = which(gC$input==trimws(colnames(dirRef)[i], which="both"))
#	tK = gC$med.microbe.other[tK0]
#	if(tK!="other"){
#		a0 = as.character(read.table(text=gC$QC.standardization[tK0], sep=";"))
#		for(i0 in 1:length(a0)){
#			if(tK=="medical"){
#				a = which(colnames(mEdic)==a0[i0]) ##!
#				mEdic[,a] = mEdic[,a] + as.numeric(dirRef[,i])
#			}else if(tK=="microbe"){
#				a = which(colnames(mIcro)==a0[i0])
#				mIcro[,a] = mIcro[,a] + as.numeric(dirRef[,i])
#}}}};rm(i,i0,a,a0,tK,tK0)

##### export #####
cat("record notation final conversion and export:",date(),"\n")
for(i in 3:ncol(mEdic)){mEdic[,i] = ifelse(mEdic[,i]>0,1,0)};rm(i)
write.csv(mEdic,"../data/cf425Medic.csv", quote=F, row.names=F)
for(i in 3:ncol(mIcro)){mIcro[,i] = ifelse(mIcro[,i]>0,1,0)};rm(i)
write.csv(mIcro,"../data/cf425Micro.csv", quote=F, row.names=F)
save(mEdic,mIcro, file="../data/cf425MedMic.rda")
cat("data sorting completed:",date(),"\n")
#for(i in 2:ncol(qcREF)){if(all(qcREF[,i]<=0)){cat(i,",",colnames(qcREF)[i],"\n")}};rm(i)
#for(i in 3:ncol(mIcro)){if(all(mIcro[,i]<=0)){cat(i,",",colnames(mIcro)[i],"\n")}};rm(i)
#for(i in 3:ncol(dirRef)){if(all(dirRef[,i]<=0)){cat(i,",",colnames(dirRef)[i],"\n")}};rm(i)
