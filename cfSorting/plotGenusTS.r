#!/bin/env Rscript
# author: ph-u
# script: plotGenusTS.r
# desc: plot time-series genus data with fitted simulations
# in: Rscript plotGenusTS.r
# out: gTS_overall.pdf
# arg: 0
# date: 20220502

##### env #####
pT = "../../00_biLVC/"
source(paste0(pT,"pipeline/src.r"))
library(deSolve)
tS = c("0820","0811","0813","1015","1113","1315","1619")
tsRaw = paRaw = vector(mode='list', length=length(tS))
for(i in 1:length(tS)){
	tsRaw[[i]] = read.csv(paste0(pT,"data/gTS_",tS[i],"_gLV-log.csv"), header=T)
	paRaw[[i]] = read.csv(paste0(pT,"data/gTS_",tS[i],"_gLV-sam.csv"), header=T)
}

##### colour #####
cBp = c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#e79f00", "#9ad0f3", "#F0E442", "#999999", "#cccccc", "#6633ff", "#00FFCC", "#0066cc", "#000000")
cBl = c("#E69F0028", "#56B4E944", "#009E7349", "#0072B233", "#D55E0033", "#CC79A733", "#e79f0033", "#9ad0f333", "#F0E44233", "#99999933", "#cccccc33", "#6633ff33", "#00FFCC33", "#0066cc33", "#00000033")
ptCol0 = rep(cBp,ceiling((ncol(tsRaw[[1]])-1)/length(cBp)))
lnCol0 = rep(cBl,ceiling((ncol(tsRaw[[1]])-1)/length(cBl)))

##### plot #####
pdf("../../thesis/fig/gTS_overall.pdf", width=14)
par(mar=c(5,5,1,12)+.1, xpd=T)

## overall data
matplot(tsRaw[[1]][,1],tsRaw[[1]][,-1], type="p", pch=(1:(ncol(tsRaw[[1]])-1))%%25, cex=1.2, col=ptCol0,
        xlab=paste0(gsub("_"," (",colnames(tsRaw[[1]])[1]),ifelse(length(grep("_",colnames(tsRaw[[1]])[1]))>0,")","")),
        ylab="percentage presence [%]", cex.axis=2, cex.lab=2)
legend("topright", inset=c(-.19,0), legend = colnames(tsRaw[[1]])[-1], pch = (1:(ncol(tsRaw[[1]])-1))%%25, lty=(1:(ncol(tsRaw[[1]])-1))%%5+1, lwd=2, col = ptCol0)

## plot each manual-checked fit simulation set
for(i in 2:length(tsRaw)){
## identified Lotka-Volterra ecological period
	segments(x0=as.numeric(paste0("20",substr(tS[i],1,2))), x1=as.numeric(paste0("20",substr(tS[i],3,4))), y0=-(length(tsRaw)-i)/2.5)

## simulations
	for(i0 in 1:nrow(paRaw[[i]])){
		a = solveLV(as.numeric(tsRaw[[i]][1,-1]), as.numeric(paRaw[[i]][i0,]), range(tsRaw[[i]][,1]), "g")
		for(i1 in 2:ncol(a)){a[,i1] = ifelse(a[,i1]>150 | a[,i1]<0,-100,a[,i1])}
		a[is.na(a)] = -100
		if(all(a!=-100)){matplot(a[,1],a[,-1], type="l", add=T, col=lnCol0, lty=(1:(ncol(tsRaw[[i]])-1))%%5+1)}
	}
}

invisible(dev.off())
