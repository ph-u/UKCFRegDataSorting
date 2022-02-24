#!/bin/env Rscript
# author: ph-u
# script: medReconstruct.r
# desc: cf425 medication time-series reconstruction
# in: Rscript medReconstruct.r
# out: none
# arg: 0
# date: 20220208

##### env #####
load("../data/cf425.rda") # 5 sec - cf425, uniqNam, yEar
clREF = read.table("../data/colDict-ExcSep.csv", sep="!", header=T, stringsAsFactors=F)
#clSUM = read.csv("../data/attributesTotal.csv", header=T, stringsAsFactors=F, row.names=1)

##### patients id #####
#pG = 1:length(cf425) # every info: 12900 patients
yR = as.numeric(substr(names(cf425)[yEar], 1, 4))
regID = yRec = c();for(i in 1:length(yEar)){
	regID = c(regID, cf425[[yEar[i]]][,grep("regid", colnames(cf425[[yEar[i]]]))])
	yRec = c(yRec, rep(yR[i],nrow(cf425[[yEar[i]]])))
}#;regID = unique(regID)

##### modify column names #####
for(i in 1:nrow(clREF)){if(clREF[i,2]=="no dictionary entry" | clREF[i,2]==""){clREF[i,2]=clREF[i,1]}} # standardise colname source
clREF[,2] = tolower(clREF[,2])
namStroke = grep("/",clREF[,2])
namCol = c();for(i in 1:nrow(clREF)){if(i %in% namStroke){
	namCol = c(namCol, as.character(read.table(text=clREF[i,2], sep="/")))
}else{namCol = c(namCol, clREF[i,2])}
};namCol0 = namCol = unique(namCol) # back-up as df record colnames (n=435)
for(i in 1:length(namCol)){if(nchar(namCol[i])>32){namCol[i]=substr(namCol[i],1,32)}} # text length limit to 32 (observed from cf425 data)

##### df record ##### ? hr
#p0=date()
pAtients = as.data.frame(matrix(NA, nr=length(regID), nc=length(namCol)+1))
colnames(pAtients) = c("year", namCol0)
pAtients$year = yRec
pAtients[,2] = regID
for(i0 in 1:nrow(pAtients)){ #cat("[",i0,": ")
	if(i0==1){ # get data source df
		refData = grep(pAtients[i0,1], names(cf425))
	}else if(pAtients[i0,1]>pAtients[i0-1,1]){
		refData = grep(pAtients[i0,1], names(cf425))
	};for(i1 in 3:ncol(pAtients)){ #cat(i1,") ")
## short-list alternative headings to map
		sH = grep(namCol[i1-1], clREF[,2])
		sH = unique(c(namCol0[i1-1],namCol[i1-1],clREF[sH,1]))
## get data (find sheet, column, row)
		i2 = 1;cL = NA;repeat{
			if(!is.na(cL) | i2>length(refData)){break}
			sS = refData[i2]
			cL = which(colnames(cf425[[sS]]) %in% sH)
			cL = ifelse(length(cL)==0,NA,cL)
			i2 = i2+1
		}
		pAtients[i0,i1] = ifelse(is.na(cL),NA,paste(cf425[[sS]][which(cf425[[sS]][,1]==pAtients[i0,2]),cL], collapse="/"))
}} #;cat("\n")
#p0;date()

save(pAtients, file="../data/remapCF425.rda")

##### write colnames csv for printing #####
p = function(i){i1 = i+1
	sH = grep(namCol[i1-1], clREF[,2])
	sH = unique(c(namCol0[i1-1],namCol[i1-1],clREF[sH,1]))
	return(sH)}
kMc = c("No.", "Standardized column names", "Remapped equivalents (separated by ``/\")")
kM = as.data.frame(matrix(NA,nr=length(namCol0), nc=length(kMc)))
colnames(kM) = kMc
for(i in 1:nrow(kM)){
	kM[i,] = read.table(text=paste(c(i,gsub("_","\\\\_",sub(" / ","&",paste(p(i), collapse=" / ")))), collapse="&"), sep="&")
	if(nchar(kM[i,3])==1){kM[i,3]="x"}
}
write.table(kM,"../../thesis/fig/cf425Cols.txt", sep="&", quote=F, row.names=F)
