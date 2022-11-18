#!/bin/env Rscript
# author: ph-u 
# script: cultureMethod.r
# desc: cf425 culture method extraction
# in: Rscript cultureMethod.r
# out: kk
# arg: 0
# date: 20221117

##### env #####
load("../data/cf425RAW.rda")

##### extract culture methods #####
cultMet = rEc[,c("regid_anon","year","s05culturetypedone","cult_type_ctr","s05numberofsputumcultures","cult_sample_sputum","s05numberofcoughswabcultures","cult_sample_throat","s05numberofbronchoscopycultures","cult_sample_bronch")]
cmNA = !is.na(cultMet)
for(i in seq(3,ncol(cultMet),2)){
	cultMet[,ncol(cultMet)+1] = ((cmNA[,i] + cmNA[,i+1])>0)
};rm(i)
colnames(cultMet)[11:14] = c("sampleLoc","sputum","coughswab","bronchoscopy")

#### proportion of samples through time #####
yR = as.numeric(unique(cultMet[,2]))
cT = c(colnames(cultMet)[2:1],colnames(cultMet)[11:14])
tSamples = as.data.frame(matrix(nr=length(yR),nc=length(cT)))
colnames(tSamples) = cT
tSamples[,1] = yR[order(yR)]
for(i in 1:nrow(tSamples)){
	tmp = cultMet[which(cultMet[,2]==tSamples[i,1]),c(1,2,(ncol(cultMet)-3):ncol(cultMet))]
	tSamples[i,-1] = c(length(unique(tmp[,1])),apply(tmp[,-(1:2)],2,sum))/nrow(tmp)
};rm(i,tmp,yR,cT)

##### plot #####
cOl = palette.colors(palette = "Okabe-Ito", alpha=1, recycle = T)
pdf("../data/cultureMethod.pdf", width=14)
par(mar=c(5,4,0,10)+.1, xpd=T)
matplot(tSamples[,1],tSamples[,-1]*100, type="l", col=cOl[1:(ncol(tSamples)-1)], lty=1, lwd=3, xlab="Year", ylab="pwCF %")
legend("topright", inset=c(-.16,0), legend=colnames(tSamples)[-1], lty=1, lwd=3, col=cOl[1:(ncol(tSamples)-1)])
dev.off()
