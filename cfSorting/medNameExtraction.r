#!/bin/env Rscript
# author: ph-u
# script: medNameExtraction.r
# desc: get standardized medication and active ingredient names
# in: Rscript medNameExtraction.r
# out: ../data/drug.csv
# arg: 0
# date: 20220726

##### env #####
p="../data/";
library(stringr)
a = read.csv(paste0(p,"medName-qcList.csv"), header=T)[,-c(2,3)];
a = a[which(!is.na(a$QC.standardization)),]
wEb = c("www.drugs.com/","www.medicines.org.uk/emc/product/","go.drugbank.com/drugs/","www.ndrugs.com/?s=","www.sdrugs.com/?c=drug&s=","pillintrip.com/medicine/","medlineplus.gov/");

##### identify where the url should be placed #####
idURL = function(URL,u=wEb){
	u0=rep(0,length(u))
	for(q in 1:length(u)){
		u0[q] = length(grep(paste0(c("\\Q",u[q]),collapse=""),URL))
	}
	if(any(u0>0)){
		u1 = which(u0>0)
		return(c(u1+1,strsplit(paste0(c("\\Q",URL),collapse=""),split=paste0(c("\\Q",u[u1]),collapse=""))[[1]][2]))
	}else{
		return(c(length(u)+2, URL))
	}
}

##### drug/active ingredient identification #####
nAm = read.table(text=unique(tolower(c(a$QC.standardization,a$func.component))), sep=";", fill=T); # ln2318
n0 = c();for(i in 1:ncol(nAm)){n0 = c(n0,nAm[,i])};n0 = unique(trimws(n0,which="both")); # ln1474

##### data recording frame #####
kCol = c("id","Drugs","Medicines","DrugBank","nDrugs","sDrugs","pillintrip","MedLinePlus","otherURL","purpose","class","function","diseases","medicationForm")
rEc = as.data.frame(matrix("",nr = length(n0),nc = length(kCol)))
colnames(rEc) = kCol
rEc$id = n0[order(n0)]
cat("Mapping url: ")
for(i in 1:nrow(rEc)){
	i0 = paste0(unique(gsub(" ","",a[unique(c(grep(rEc[i,1],tolower(a$QC.standardization)),grep(rEc[i,1],a$func.component))),"url"])), collapse=";")
	if(length(grep(";",i0))>0){
		p0 = strsplit(i0,split=";")
		p1 = c();for(q in 1:length(p0)){p1 = c(p1,p0[[q]])};
		p1 = unique(trimws(p1, which="both"))
		for(q in 1:length(p1)){
			i1 = idURL(p1[q])
			rEc[i,as.numeric(i1[1])] = paste0(c(rEc[i,as.numeric(i1[1])],i1[2]),collapse=";")
			rEc[i,as.numeric(i1[1])] = sub("^;","",rEc[i,as.numeric(i1[1])])
		}
	}else{
		i1 = idURL(i0)
		rEc[i,as.numeric(i1[1])] = i1[2]
	}
	if(i%%100==0){cat(i,", ")}
};cat("\n")

##### export #####
write.csv(rEc, paste0(p,"drug.csv"), quote=F, row.names=F);
