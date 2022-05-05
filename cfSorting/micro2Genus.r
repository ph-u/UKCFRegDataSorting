#!/bin/env Rscript
# author: ph-u
# script: micro2Genus.r
# desc: sort mIcro data in genus time-series
# in: Rscript micro2Genus.r [% presence threshold]
# out: genusCF425_gLV.csv, cf425Genus.pdf
# arg: 1
# date: 20220428

##### f: import csv with colnames converted in standard format #####
iPt = function(x){
	x0 = read.csv(paste0("../data/",x,".csv"), header=T, stringsAsFactors=F)
	colnames(x0)[grep("[.][.]",colnames(x0))] = gsub("[.]$",")",gsub("[.][.]"," (",colnames(x0)[grep("[.][.]",colnames(x0))]))
	colnames(x0) = gsub("[.]"," ",colnames(x0))
	return(x0)
}

##### env #####
argv = (commandArgs(T))
library(stringr)
mIcro = iPt("cf425Micro")#;mEdic = iPt("cf425Medic")
#load("../data/cf425MedMic.rda")

##### resort according to the list of genus #####
gEn = word(colnames(mIcro)[-c(1,2)],1)
uqGen = unique(gEn)
mcGen = as.data.frame(matrix(0,nr=nrow(mIcro),nc=2+length(uqGen)))
mcGen[,1:2] = mIcro[,1:2];colnames(mcGen) = c(colnames(mIcro)[1:2],uqGen)
for(i in 1:length(gEn)){
	tMp = grep(gEn[i],colnames(mcGen))
	mcGen[,tMp] = mcGen[,tMp] + mIcro[,2+i]
};for(i in 3:ncol(mcGen)){mcGen[,i] = ifelse(mcGen[,i]>0,1,0)}
## info check
#tMp = sample(1:nrow(mIcro),30,replace=F);i=1
#mIcro[tMp[i],which(mIcro[tMp[i],]!=0)];mcGen[tMp[i],which(mcGen[tMp[i],]!=0)];i = min(i+1,length(tMp))

##### time-series conversion (unit: fraction) #####
yR = unique(mIcro$year);yR = yR[order(yR)]
mcGenTS = as.data.frame(matrix(0,nr=length(yR),nc=ncol(mcGen)-1))
colnames(mcGenTS) = colnames(mcGen)[-1]
mcGenTS$year = yR
for(i0 in 1:nrow(mcGenTS)){
	i1 = mcGen[which(mcGen$year==mcGenTS$year[i0]),-c(1:2)]
	mcGenTS[i0,-1] = apply(i1,2,sum)/nrow(i1)
};rm(i0,i1)

##### plot presence fractions #####
tHs = as.numeric(argv[1])/100 ## 1% presence threshold (arbitrary from data observation)
pdf("../../thesis/fig/cf425Genus.pdf", width=21, height=7)
matplot(1:(ncol(mcGenTS)-1),t(mcGenTS)[-1,], type="l", xlab="Microbial genus presence in patients (alphabetical order)", ylab="proportion of patients", ylim=c(-.25,.7), cex.lab=1.4)
abline(h=tHs); text(0,tHs+.05,paste("threshold =",tHs))
legend("topleft",legend = mcGenTS$year, col = 1:nrow(mcGenTS), lty=1:nrow(mcGenTS), lwd=2)
xRef = colnames(mcGenTS)[-1]
for(i in 2:ncol(mcGenTS)){if(all(mcGenTS[,i]<tHs)){
	xRef[i-1] = "others"
}else{
	text(i-1,-.15,colnames(mcGenTS)[i],srt=90)
}}
invisible(dev.off())

##### combine non-common species per patient per year #####
mcCmb = as.data.frame(matrix(0,nr=nrow(mcGen),nc=1+length(unique(xRef))))
colnames(mcCmb) = c(colnames(mcGen)[2],unique(xRef))
mcCmb[,1] = mcGen[,2]
for(i in 1:length(xRef)){
	mcCmb[,xRef[i]] = mcCmb[,xRef[i]] + mcGen[,i+2]
};for(i in 2:ncol(mcCmb)){mcCmb[,i] = ifelse(mcCmb[,i]>0,1,0)}

##### sorted time-series #####
mcGenTSf = as.data.frame(matrix(0,nr=nrow(mcGenTS),nc=ncol(mcCmb)))
colnames(mcGenTSf) = colnames(mcCmb)
mcGenTSf$year = mcGenTS$year
for(i0 in 1:nrow(mcGenTSf)){
        i1 = mcCmb[which(mcCmb$year==mcGenTSf$year[i0]),-1]
        mcGenTSf[i0,-1] = apply(i1,2,sum)/nrow(i1)
};rm(i0,i1)

##### export #####
write.csv(mcGenTSf,"../data/genusTimeSeries_gLV.csv", quote=F, row.names=F)
