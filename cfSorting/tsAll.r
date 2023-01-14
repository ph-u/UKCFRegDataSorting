#!/bin/env Rscript
# author: ph-u
# script: tsAll.r
# desc: plot time-series genus data with fitted simulations
# in: Rscript tsAll.r [medication group code]
# out: cf425_overall.pdf
# arg: 0 [for i in `ls ../F508del_data/*-eco* | cut -f 3 -d "_" | uniq`;do Rscript tsAll.r ${i};done]
# date: 20230113

##### env #####
argv = (commandArgs(T))
source("../../00_biLVC/src/src.r")
library(deSolve)
nRep = 7; simO = 500
ptIN = "../F508del_data/"
ptOT = "../graph/"
acRatio = .25 # acceptance ratio in plotting
cBp = palette.colors(palette = "Okabe-Ito", alpha=1, recycle = T)
cBl = palette.colors(palette = "Okabe-Ito", alpha=.01, recycle = T)
rEf = data.frame(code=c("0000","0010","0100","0110","1010","1011","1110","1111"), group=c("No medications", "Drugs/supplements", "Antimicrobials", "Antimicrobials + Drugs/supplements", "CFTR modulators + Drugs/supplements", "CFTR modulators + Drugs/supplements + Drug-drug interaction", "CFTR modulators + Antimicrobials + Drugs/supplements", "CFTR modulators + Antimicrobials + Drugs/supplements + Drug-drug interaction"))

##### f: capitalise first letter #####
capFirst = function(x){return(paste0(toupper(substr(x,1,1)),substr(x,2,nchar(x))))}

##### f: plot legend #####
legPlot = function(x,nDim=3){
        nDim = c(nDim, ceiling(length(x)/nDim))
        legMx = matrix(1:prod(nDim), nrow=nDim[1], ncol=nDim[2], byrow=F)
        legBd = rep("#000000ff", prod(nDim))
        legBd[legMx>length(x)] = legMx[legMx>length(x)] = NA
        return(list(legMx,nDim[2]))
}

##### treatment distribution along time #####
cat("Sample size [",argv[1],"] calculation: ",date(),"\n")
samSize = read.csv("../graph/cftrm_samFreq.csv", header=T)
x = c("Year","total","CFTR_modulators","antimicrobials","drugs/supplements","drug_interaction")
tReat = as.data.frame(matrix(nr=nrow(samSize), nc=length(x)))
colnames(tReat) = x;rm(x)
tReat$Year = samSize$year
tReat$total = apply(samSize[,-1],1,sum)
for(i in 3:ncol(tReat)){
        tReat[,i] = apply(samSize[,which(substr(colnames(samSize),i-1,i-1)==1)],1,sum)
};rm(i)

##### time-series data #####
cat("Time-series data [",argv[1],"]: ",date(),"\n")
f = list.files(ptIN,"log")
f = f[grep(argv[1],f)]
for(i in 1:length(f)){
	if(i==1){f0 = read.csv(paste0(ptIN,f[i]), header=T)}else{f0 = merge(f0,read.csv(paste0(ptIN,f[i]), header=T), all=T)}
};rm(i)

## plot
pdf(paste0(ptOT,"cf425_",argv[1],"_overall.pdf"), width=14, height=14)
#pdf(paste0(ptOT,"cf425_overall.pdf"), width=14, height=14)
par(mar=c(5,5,1,0)+.1, mfrow = c(2,1), cex.axis=1.4, xpd=F)
matplot(f0[,1], f0[,-1], "p", pch=4+1:(ncol(f0)-1), xaxt="n", xlab="", ylab="Prevalence (%)", cex=5, col=cBp, cex.axis=2.1, cex.lab=2.1, ylim = c(0,100), xlim = c(2008,2020))
abline(h=0, col="#000000ff")
axis(1, at=seq(2008,2020), padj=.7, labels=paste0(seq(2008,2020),"\n(",samSize[,which(colnames(samSize)==paste0("g",argv[1]))],")\n[",tReat[,2],"]"))
mtext(paste0("Year (Group Sample Size) [Total Sample Size]\n",rEf$group[which(rEf$code==argv[1])]),side=1,padj=2.8,cex=2.1)
text(2007.75, 93, labels="n =", cex=2)

##### Simulation Time-series #####
cat("Simulation [",argv[1],"] plotting: ",date(),"\n")
p = list.files(ptIN,"filter")
p = p[grep(argv[1],p)]

for(i in 1:length(p)){
	iNfo = c(strsplit(p[i],"_")[[1]][3],strsplit(p[i],"-")[[1]][2])
	t0 = as.numeric(c(substr(iNfo[1],1,2), substr(iNfo[1],3,4)))+2000
	x0 = as.numeric(f0[which(f0[,1]==paste0("20",substr(iNfo[1],1,2))),-1])
	oDe = strsplit(strsplit(p[i],"_")[[1]][4],"-")[[1]][1]
## simulation
	a = read.csv(paste0(ptIN,p[i]), header=T)
	if(nrow(a)>0){
		s = read.csv(paste0(ptIN,sub("filter","seed",p[i])), header=T)
		for(i0 in 1:nrow(a)){
			set.seed(s[a[i0,1],2])
			a0 = solveLV(x0, as.numeric(a[i0,-1]), t0, oDe)
			matplot(a0[,1], a0[,-1], type="l", add=T, col=cBl, lty = (1:(ncol(a0)-1))%%5+1)
		};rm(i0)
	}; text(t0[1]+.1,ifelse(t0[1]%%2==0, 97,95), labels=nrow(a), cex=2)
};rm(i,iNfo)

plot.new()
lPt = legPlot(colnames(f0)[-1], nDim = 3)
legend(x=c(-.12,1.04),y=c(.85,.45), legend = capFirst(colnames(f0)[-1]), title = paste("Microbial categories:",nRep,"replicates,",simO,"simulations each; N =",nRep*simO), border=NA, xpd=T, ncol=lPt[[2]], pch = 4+1:(ncol(f0)-1), col = cBp, cex=2.1)
invisible(dev.off())
