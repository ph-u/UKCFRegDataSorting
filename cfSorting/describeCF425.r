#!/bin/env Rscript
# author: ph-u
# script: describeCF425.r
# desc: cf425 population time-series & data collection summary
# in: Rscript describeCF425.r
# out: ../../thesis/fig/tsPatients.pdf ../data/attributesTotal.csv
# arg: 0
# date: 20220201

##### env #####
load("../data/cf425.rda") # 5 sec - cf425, uniqNam, yEar

## year as number
yNum = read.table(text=rev(names(cf425)), sep=" ", row.names=NULL, header=F, fill=T, stringsAsFactors=F)[,1]
yNum = unique(yNum[which(nchar(yNum)==4)])
yNum = as.numeric(yNum[order(yNum)])

## patients
regID=c();for(i in 1:length(yEar)){
	regID = c(regID,cf425[[yEar[i]]][,"regid_anon"])
};regID = unique(regID) # 12727 patients
cOl = c("regid", "startYear", "lastCheck","timepoints")
stEnd = as.data.frame(matrix(NA, nr=length(regID), nc=length(cOl)))
colnames(stEnd) = cOl
stEnd$regid = regID
stEnd$timepoints = 0

## attributes
dtAttri=c();for(i in 1:length(cf425)){
	dtAttri = c(dtAttri, colnames(cf425[[i]]))
};dtAttri = unique(dtAttri) # 572 attributes
attrCol = as.data.frame(matrix(NA,nr=length(yEar), nc=length(dtAttri)))
colnames(attrCol) = dtAttri
row.names(attrCol) = yNum

##### time-series details #####
for(i in 1:length(yEar)){
	stEnd[which(is.na(stEnd[,2]) & regID %in% cf425[[yEar[i]]]$regid_anon),2] = yNum[i] # start year (patient)
	stEnd[which(is.na(stEnd[,3]) &  regID %in% cf425[[yEar[length(yNum)-i+1]]]$regid_anon),3] = yNum[length(yNum)-i+1] # end year (patient)
	stEnd[which(regID %in% cf425[[yEar[i]]]$regid_anon),4] = stEnd[which(regID %in% cf425[[yEar[i]]]$regid_anon),4]+1 # timepoint count (patient)
}
stEnd$duration = paste0(stEnd[,2],"-",stEnd[,3])

##### patients time-series type #####
durType = as.data.frame(table(stEnd$duration))
durType = cbind(durType,read.table(text=as.character(durType[,1]), sep="-"))[,c(3,4,2)]
colnames(durType) = c("start","last","freq")
durType = durType[order(durType[,3]),]
cOl = ifelse(durType[,2]-durType[,1]>2,"#000000ff", ifelse(durType[,2]-durType[,1]<2, "#ff0000ff", "#0000ffff"))
durType0 = as.data.frame(t(durType))
yAx = durType0[-3,]
yAx[1,]=yAx[2,]=1:ncol(durType0)

##### patients time-series frequency plot #####
pdf("../../thesis/fig/tsPatients.pdf", width=7, height=14)
par(mar=c(5,0,0,0)+.1)
matplot(durType0[-3,],yAx, type="l", lty=1, yaxt="n", col="#00000077", ylab="", xlab="year")
text((durType0[1,]+durType0[2,])/2, 1:ncol(durType0)+.1, labels=durType0[3,], col=cOl)
invisible(dev.off())

##### Trim cf425 #####
seFilter = stEnd[which(stEnd[,4]>2),] # time-series = >2 timepoints; 11507 patients
#median(seFilter[,4]) # median 11 timepoints

##### attribute check #####
#for(i in 1:ncol(attrCol)){cat(i," -- ",colnames(attrCol)[i], ": ", length(which(!is.na(attrCol[,i]))),"\n")} # 13 years * 405 unique col
for(i0 in 1:nrow(attrCol)){
	i2 = grep(row.names(attrCol)[i0],names(cf425))
	for(i1 in 1:ncol(attrCol)){ i3=1;while(i3<=length(i2) & is.na(attrCol[i0,i1])){
		i4 = which(colnames(cf425[[i2[i3]]]) %in% colnames(attrCol)[i1])
			if(length(i4)>0){
				attrCol[i0,i1] = class(cf425[[i2[i3]]][,i4])
			}else{i3=i3+1}
	}}}
write.csv(attrCol, "../data/attributesTotal.csv", quote=F, row.names=T)
#i="s05culturespeciesotherpseudomona"
#i="cult_aeruginosa_mucoid"
#i="s05cultpseudoaeruginosasamples"
#i0=which(attrCol[,i]>0)
#boxplot(cf425[[yEar[i0[1]]]][,i])
