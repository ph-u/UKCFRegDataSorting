#!/bin/env Rscript
# author: ph-u
# script: F508delBetaDiv.r
# desc: calculate diversity dissimilarity between microbial pairs among all medication groups
# in: Rscript F508delBetaDiv.r
# out: medPairPCA.pdf, medPLSDA.pdf
# arg: 0
# date: 20230111

##### env #####
library(ggbiplot)
nRep = 7
simO = 500
ptIN = "../F508del_data/"
ptOT = "../graph/F508del/"
cBp = palette.colors(palette = "Okabe-Ito", alpha=1, recycle = T)

##### files #####
cat("Import file list: ",date(),"\n")
f = list.files(ptIN,"-eco")
f0 = unlist(strsplit(f,"_"))
fNam = as.data.frame(matrix(nr=length(f),nc=2))
for(i in 1:nrow(fNam)){fNam[i,] = f0[(2:3)+length(f0)/length(f)*(i-1)]};rm(i,f0)
fNam$start = substr(fNam[,2],1,2); fNam$end = substr(fNam[,2],3,4)
f0 = unique(fNam[,1])
yR = as.numeric(unique(c(fNam$start,fNam$end)))+2000;yR = yR[order(yR)]

##### ecology count #####
cat("Summarizing simulations: ",date(),"\n")
for(i in 1:length(f)){
	f1 = read.csv(paste0(ptIN,f[i]), header=T, stringsAsFactors=F)
	x = which(f0==strsplit(f[i],"_")[[1]][2])
	if(i==1){
		eCo = unique(f1$c1_is)
		rEc = as.data.frame(matrix(0, nr=length(f0), nc=length(eCo)))
		pAir = length(sP <- unique(paste0(f1$category1,"-",f1$category2)))
		colnames(rEc) = eCo; rm(eCo)
		row.names(rEc) = c("No medications", "Drugs/supplements", "Antimicrobials", "Antimicrobials + Drugs/supplements", "CFTR modulators + Drugs/supplements", "CFTR modulators + Drugs/supplements +\nDrug-drug interaction", "CFTR modulators +\nAntimicrobials + Drugs/supplements", "CFTR modulators + Antimicrobials +\n     Drugs/supplements + Drug-drug interaction")
		#row.names(rEc) = paste0("g",f0)
		spP = as.data.frame(matrix(0, nr=pAir*length(f0), nc=ncol(rEc)+2))
		colnames(spP) = c("pairs","group",colnames(rEc))
		spP[,1] = sP; spP[,2] = rep(f0, each=pAir)
	}; for(i0 in 1:ncol(rEc)){
		rEc[x,i0] = rEc[x,i0] + sum(f1$count[which(f1$c1_is==colnames(rEc)[i0])])
		i2 = paste0(f1$category1,"-",f1$category2)
		for(i1 in 1:length(sP)){
			i3 = which(spP[,1]==sP[i1] & spP[,2]==strsplit(f[i],"_")[[1]][2])
			spP[i3,i0+2] = spP[i3,i0+2] + sum(f1$count[which(i2==sP[i1] & f1$c1_is==colnames(spP)[i0+2])])
		};rm(i1)
	};rm(i0)
};rm(i,x, f1)
rEc$unpredictable = nRep * simO * pAir * table(fNam$V1) - rowSums(rEc)
spP$unpredictable = nRep * simO * rep(table(fNam$V1),each=pAir) - rowSums(spP[,-c(1:2)])

##### dissimilarity #####
cat("Clustering: ",date(),"\n")
#hC = hclust(dist(rEc))
#pdf(paste0(ptOT,"medCluster.pdf"), width=14, height=21)
#plot(hC, main="", xlab="", xaxt="n", yaxt="n", ylab="", cex=2.8, lwd=10)
#segments(x0 = 1:2*4, y0 = rep(5e4,2), y1=rep(1.8e5,2), lty=1, lwd=21, col=cBp[7])
#polygon(x = rep(c(3.5,4.4),each=2), y = c(.54,-17,-17,.54)*1e5, lwd=21, border=cBp[7])
#polygon(x = rep(c(7.7,8.2),each=2), y = c(.54,-14,-14,.54)*1e5, lwd=21, border=cBp[7])
#invisible(dev.off())

##### PCA #####
rEf = data.frame(ori=colnames(spP)[-c(1:2)], mod=paste0("<-",c("Mutualism", "Commensalism", "Predatory/Parasitism", "Commensalism", "Neutral/No interaction", "Amensalism", "Predatory/Parasitism", "Amensalism", "Competition", "Unpredictable"),"-"))
sP0 = as.data.frame(matrix(0, nr=nrow(spP), nc=length(unique(rEf[,2]))+2))
colnames(sP0) = c(colnames(spP)[1:2],unique(rEf[,2]))
sP0[,1:2] = spP[,1:2]
for(i in 3:ncol(spP)){sP0[,rEf[which(rEf[,1]==colnames(spP)[i]),2]] = sP0[,rEf[which(rEf[,1]==colnames(spP)[i]),2]] + spP[,i]};rm(i)

rEf0 = data.frame(group = unique(sP0[,2]),name = c("No medication", "Drugs/supplements", "Antimicrobials", "Antimicrobials + Drugs/supplemtns", "CFTR modulators + Drugs/supplements", "CFTR modulators + Drugs/supplements + Drug-drug interaction", "CFTR modulators + Antimicrobials + Drugs/supplements", "CFTR modulators + Antimicrobials + Drugs/supplements + Drug-drug interaction"))
for(i in 1:nrow(rEf0)){sP0[which(sP0[,2]==rEf0[i,1]),2] = rEf0[i,2]};rm(i)

p = prcomp(sP0[,-c(1:2)], scale.=T)
pCa0 = ggbiplot(p, var.scale = 1,
	groups=sP0[,2], ellipse = TRUE, ellipse.prob = .95,
	labels = row.names(sP0), labels.size = 4,
	varname.size = 4, varname.adjust = c(5,3.5,3.3,3.5,4,4.5,4.5)-1
)+scale_color_manual(values=setNames(cBp[1:length(unique(sP0[,2]))], unique(sP0[,2])))+
guides(color = guide_legend(title="Treatment"))+theme_bw()+ylim(c(-10,7))+
coord_cartesian(xlim=c(-5,2))+
theme(legend.position = 'bottom', legend.direction = "vertical",
	legend.title = element_text(size=18),
	axis.text=element_text(size=16),axis.title=element_text(size=14),
	plot.margin=margin(t=0,r=0,b=0,l=1))
ggsave(paste0(ptOT,"medPairPCA.pdf"), plot=pCa0, width=6, height=7) # last_plot()

#p0 = summary(p)$importance[2,]
#pdf(paste0(ptOT,"medPairPCA.pdf"), width=14, height=14)
#biplot(p, cex=2.1, lwd=7, cex.axis = 2.1, cex.lab=2.1, xlab=paste0("PC1 (",round(p0[1],3)*100,"%)"), ylab=paste0("PC2 (",round(p0[2],3)*100,"%)"))
#invisible(dev.off())

##### PLS-DA #####
library(mixOmics) # citation('mixOmics'), v6.18.1
for(i in 1:nrow(rEf0)){sP0[which(sP0[,2]==rEf0[i,2]),2] = rEf0[i,1]};rm(i)
p0 = plsda(sP0[,-c(1:2)], as.factor(sP0[,2]))
#write.csv(sP0[c(252,232,52,8,452,163,150,184,376),1:2], paste0(ptOT, "plsda_outliers.csv"), quote=F, row.names=T)
pdf(paste0(ptOT,"medPLSDA.pdf"), width=14, height=14)
plotIndiv(p0, ellipse=T, ellipse.level=.95, col=cBp[1:nrow(rEf0)], size.xlabel=rel(3), size.ylabel=rel(3), size.axis=rel(2), legend=T, style="graphics", cex=2) # cBp[c(3,4,7,8,5,6,2,1)]
invisible(dev.off())
