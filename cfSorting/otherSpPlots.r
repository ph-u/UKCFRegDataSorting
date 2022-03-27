#!/bin/env Rscript
# author: ph-u
# script: otherSpPlots.r
# desc: make preliminary reference frame for other species
# in: Rscript otherSpPlots.r
# out: otherSp-nameList.csv, GoogleSearchEfficiency.pdf, DataIrregularity.pdf
# arg: 0
# date: 20220226

library(stringdist)
a = read.csv("../data/otherSp-qcList.csv", header=T, stringsAsFactors=F)
lAbel = c("automated", "hybrid", "manual")
cBp = c("#E69F00ff", "#56B4E9ff", "#009E73ff", "#0072B277", "#D55E0077", "#CC79A777", "#e79f0077", "#9ad0f333", "#F0E44233", "#99999933", "#cccccc33", "#6633ff33", "#00FFCC33", "#0066cc33", "#00000033")

##### automated Google search efficiency #####
# categorize Google help
a$similar = stringsim(a[,2],a[,4])
a$auto.hyb.man = ifelse(a$similar==1,"a",ifelse(a$similar<.3,"m","h")) # 0.3 threshold set by observation

pdf("../../thesis/fig/GoogleSearchEfficiency.pdf")
barplot(table(a$auto.hyb.man)/nrow(a), col="#00000000", ylab="proportion", xlab="Google search quality check", xaxt="n", cex.axis=1.2, cex.lab=1.2, ylim=c(0,.4))
axis(1,at=c(.7,1.9,3.1),labels=lAbel)
invisible(dev.off())

##### data typos/irregularity quantification #####
sP = unique(read.table(text=a$QC.standardization, sep=";")[,1])
txUq = as.data.frame(matrix(0,nc=length(lAbel),nr=length(sP)))
colnames(txUq) = c("a","h","m");row.names(txUq) = sP
for(i in 1:nrow(a)){
	tMp0 = as.character(read.table(text=a[i,4],sep=";")[1,])
	txUq[tMp0,a[i,7]] = txUq[tMp0,a[i,7]] + 1
};rm(i,tMp0)
lIm=20
sP0 = apply(txUq,1,sum)
txAhm = txUq[which(sP0>lIm),]
txAhm = txAhm[rev(order(apply(txAhm,1,sum))),]
colnames(txAhm) = lAbel

pdf(paste0("../../thesis/fig/DataIrregularity.pdf"), height=9,width=7)
par(mar=c(4.9,14,0,3)+.1)
barplot(as.matrix(t(txAhm)),col=cBp[1:ncol(txAhm)], border=F, space=.04, xlab="frequency", cex.axis=1.2, cex.lab=1.2, hori=T, las=1, xlim=c(0,ceiling(max(apply(txAhm,1,sum))/10)*10))
legend("topright", legend=lAbel, fill=cBp[1:ncol(txAhm)], title="Google search quality check", border=NA, xpd=T, cex=1.2)
#barplot(txUq[txUq>20], col="#000000ff", border=F, xlab="frequency", cex.axis=1.2, cex.lab=1.2, hori=T, las=1)
invisible(dev.off())

txFq = apply(txUq[which(sP0<=lIm),],1,sum)
txDF = as.data.frame(matrix("",nr=max(table(txFq)),nc=max(txFq)))
for(i in 1:ncol(txDF)){
	zP = names(txFq)[which(txFq==i)]
	txDF[,i] = c(zP, rep("",nrow(txDF)-length(zP)))
};rm(i,zP)
write.csv(txDF,"../data/otherSp-nameList.csv", quote=F, row.names=F)
