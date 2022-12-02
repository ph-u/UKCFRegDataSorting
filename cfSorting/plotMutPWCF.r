#!/bin/env Rscript
# author: ph-u
# script: plotMutPWCF.r
# desc: plot CF mutation distribution
# in: Rscript plotMutPWCF.r
# out: ../graph/mutPWCF.pdf
# arg: 0
# date: 20221202

##### import #####
a = read.csv("../data/mutPWCF.csv", header=T)

##### percentage distribution calculation #####
a0 = apply(a[,-c(1:2)],2,sum)/nrow(a)*100
pdf("../graph/mutPWCF.pdf", width=14)
par(mar=c(5,5,0,0)+.1)
plot(1:(ncol(a)-2),a0, xlab="Mutation category code", ylab="Prevalence (%)",pch=4,cex.axis=2, cex.lab=2)
invisible(dev.off())
