#!/bin/env Rscript
# author: ph-u
# script: colRefCF425.r
# desc: CF-Registry column dictionary xlsx reorganize
# in: Rscript colRefCF425.r
# out: ../data/cf425Dict.csv
# arg: 0
# date: 20220203

##### in #####
library(xlsx2dfs)
rEf = xlsx2dfs("../p_cf425/UK CF Registry Data dictionary - v4.2 Oct 2021.xlsx", rowNames=F, detectDates = T)
a = read.csv("../data/attributesTotal.csv", row.names=1, header=T)
rEf1 = read.csv("../data/cf425_extraDataGuide.csv", header=T)

##### ref df important columns #####
cOmmon = c("Field.name","Field.question","Mapping.from.PortCF")
iMp = vector(mode="list")
iMp$Demographic = iMp$"Annual Review" = iMp$"Chronic Medications" = cOmmon
iMp$Medications = c("DrugName","DrugType")
iMp$"PortCF medication mapping" = c(1,6:16)
iMp$Vaccinate = iMp$Culture = iMp$AcuteCOVID = iMp$"Tendon Rupture" = iMp$NTMCulture = cOmmon[-3]
iMp$Transplant = c(cOmmon[-3], "remap.from.PortCF")

##### record df ##### 572 attributes
a0Nam = c("attribute", "equivalent", "meaning", "spreadsheet", "row")
a0 = as.data.frame(matrix(NA, nr=length(colnames(a)), nc=length(a0Nam)))
colnames(a0) = a0Nam
a0[,1] = colnames(a)

##### refined ref dictionary #####
rEf0 = vector(mode="list")
for(i in 1:length(iMp)){
	rEf0[[i]] = rEf[[names(iMp)[i]]][,iMp[[i]]]
}
names(rEf0) = names(iMp)

##### dictionary/questions mapping (572 attributes) #####
pS = function(tX){return(paste(tX, collapse="/"))}
for(i in 1:nrow(a0)){
	sS = grep(tolower(a0[i,1]), tolower(rEf0))
	if(length(sS)>1 & any(sS)==1){sS=1} # 109 attributes
	if(length(sS)==1){ # 390 attributes
		cL = grep(tolower(a0[i,1]), tolower(rEf0[[sS]]))
		if(length(cL)>1){cL = grep(paste0(tolower(a0[i,1])," "), tolower(rEf0[[sS]]))}
		rW = grep(tolower(a0[i,1]), tolower(rEf0[[sS]][,cL]))
		if(names(rEf0)[sS]=="PortCF medication mapping"){
			if(length(rW)>1){
				eQ = c();for(i0 in 1:length(rW)){
					eQ = c(eQ,as.character(rEf0[[sS]][rW[i0],-1]))
				}}else{
				eQ = as.character(rEf0[[sS]][rW,-1])
			}
			eQ = unique(eQ[which(eQ!="NULL" & eQ!="n/a" & !is.na(eQ))])
			mE = unique(rEf0[[sS]][rW,iMp[[sS]][1]])
		}else{
			eQ = rEf0[[sS]][rW,iMp[[sS]][1]]
			mE = rEf0[[sS]][rW,iMp[[sS]][2]]
		};nAm = names(rEf0)[sS]
	}else{ # 73 attributes
		rW = which(rEf1==a0[i,1])
		if(length(rW)<1){eQ = mE = nAm = rW = "no dictionary entry"}else{
			eQ = unique(rEf1[rW,1])
			mE = unique(rEf1[rW,2])
			nAm = unique(rEf1[rW,3])
		}
	}
	a0[i,-1] = c(pS(eQ), pS(mE), pS(nAm), pS(rW))
}

##### replace newline #####
for(i in 1:ncol(a0)){
	nL = grep("\n", a0[,i])
	if(length(nL)>0){ for(i0 in 1:length(nL)){
		a0[nL[i0],i] = gsub("\n"," ", a0[nL[i0],i])
	}}}

write.table(a0, "../data/colDict-ExcSep.csv", sep="!", row.names=F, quote=F, col.names=T)
