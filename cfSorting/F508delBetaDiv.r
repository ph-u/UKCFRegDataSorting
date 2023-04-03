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
ptIN = "../p_raw_20230207/data/" #"../F508del_data/"
ptOT = "../p_raw_20230207/res/" #"../graph/F508del/"
cBp = palette.colors(palette = "Okabe-Ito", alpha=1, recycle = T)
m0 = read.csv(paste0(ptIN,"../medDist.csv"), header=T)
m1 = read.csv(paste0(ptIN,"../F508dd_samFreq.csv"), header=T)
x0 = c(colnames(m0)[grep("CFTR",colnames(m0))],colnames(m0)[grep("antim",colnames(m0))],colnames(m0)[grep("panc",colnames(m0))],"others")

##### f: capitalise first letter #####
capFirst = function(x){return(gsub("[.]"," ",paste0(toupper(substr(x,1,1)),substr(x,2,nchar(x)))))}

##### files #####
cat("Import file list: ",date(),"\n")
f = list.files(ptIN,"-eco"); f = f[-grep("1101",f)] # group 1101 has 1 small sample size, unmodellable time-series
f0 = unlist(strsplit(f,"_"))
fNam = as.data.frame(matrix(nr=length(f),nc=2))
for(i in 1:nrow(fNam)){fNam[i,] = f0[(2:3)+length(f0)/length(f)*(i-1)]};rm(i,f0)
fNam$start = substr(fNam[,2],1,2); fNam$end = substr(fNam[,2],3,4)
f0 = unique(fNam[,1])
mNam = rep(NA,length(f0)); for(i in 1:length(mNam)){
	i0 = rep(NA,length(x0))
	for(i1 in 1:length(x0)){i0[i1] = capFirst(ifelse(substr(f0[i],i1,i1)==0,"",x0[i1]))}
	mNam[i] = paste(i0[i0!=""], collapse=" + ")
};rm(i,i0,i1)

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
		#row.names(rEc) = c("No medications", "Drugs/supplements", "Antimicrobials", "Antimicrobials + Drugs/supplements", "CFTR modulators + Drugs/supplements", "CFTR modulators + Drugs/supplements +\nDrug-drug interaction", "CFTR modulators +\nAntimicrobials + Drugs/supplements", "CFTR modulators + Antimicrobials +\n     Drugs/supplements + Drug-drug interaction")
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
rEf = data.frame(ori=colnames(spP)[-c(1:2)], mod=paste0("<-",c("Mutualism", "Commensalism", "Predatory/Parasitism", "Commensalism", "Neutral/No interaction", "Amensalism", "Predatory/Parasitism", "Amensalism", "Competition", "Unpredictable"),"->"))
sP0 = as.data.frame(matrix(0, nr=nrow(spP), nc=length(unique(rEf[,2]))+2))
colnames(sP0) = c(colnames(spP)[1:2],unique(rEf[,2]))
sP0[,1:2] = spP[,1:2]
for(i in 3:ncol(spP)){sP0[,rEf[which(rEf[,1]==colnames(spP)[i]),2]] = sP0[,rEf[which(rEf[,1]==colnames(spP)[i]),2]] + spP[,i]};rm(i)

rEf0 = data.frame(group = unique(sP0[,2]),name = NA)#c("No medication", "Drugs/supplements", "Antimicrobials", "Antimicrobials + Drugs/supplements", "CFTR modulators + Drugs/supplements", "CFTR modulators + Drugs/supplements + Drug-drug interaction", "CFTR modulators + Antimicrobials + Drugs/supplements", "CFTR modulators + Antimicrobials + Drugs/supplements + Drug-drug interaction"))
for(i in 1:nrow(rEf0)){rEf0$name[i] = mNam[which(f0==rEf0$group[i])]};rm(i)
for(i in 1:nrow(rEf0)){sP0[which(sP0[,2]==rEf0[i,1]),2] = rEf0[i,2]};rm(i)

p = prcomp(sP0[,-c(1:2)], scale.=T)
pCa0 = ggbiplot(p, var.scale = 1,
	groups=sP0[,2], ellipse = TRUE, ellipse.prob = .95,
	labels.size = 4, #labels = row.names(sP0),
	varname.size = 4, varname.adjust = c(3,2.4,2,1.5,2.5,3.2,1.8) # mut,comm,,,ame,comp,unp
)+#xlab("")+ylab("")+
scale_color_manual(values=setNames(cBp[1:length(unique(sP0[,2]))], unique(sP0[,2])))+
guides(color = guide_legend(title="Treatment"))+theme_bw()+#ylim(c(-4,9))+xlim(c(-20,20))+
scale_y_continuous(breaks = seq(-3, 9, 3), limits = c(-3.5, 9))+
coord_cartesian(xlim=c(-3.5,5))+
theme(legend.position = 'bottom', legend.direction = "vertical",
	panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
	legend.title = element_text(size=18),
	axis.text=element_text(size=16),axis.title=element_text(size=14),
	plot.margin=margin(t=0,r=0,b=0,l=1))
ggsave(paste0(ptOT,"medPairPCA.pdf"), plot=pCa0, width=6, height=7) # last_plot()
#ggsave(paste0(ptOT,"medPairPCA.jpg"), plot=pCa0, width=6, height=7) # last_plot()

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
#jpeg(paste0(ptOT,"medPLSDA.jpg"), width=1000, height=1000)
plotIndiv(p0, pch=20, ellipse=T, ellipse.level=.95, col=cBp[1:nrow(rEf0)], size.xlabel=rel(1.4), size.ylabel=rel(4), size.axis=rel(2), legend=F, style="graphics", cex=3) # cBp[c(3,4,7,8,5,6,2,1)]
#plotIndiv(p0, ellipse=T, ellipse.level=.95, col=cBp[1:nrow(rEf0)], X.label ="", Y.label ="", size.ylabel=rel(4), size.axis=rel(2), legend=F, style="graphics", cex=3) # cBp[c(3,4,7,8,5,6,2,1)]
invisible(dev.off())

## Loading plots
#write.csv(p0$loadings$X, paste0(ptOT,"PLSDA_ecoLoadings.csv"), quote=F)
#write.csv(p0$loadings$Y, paste0(ptOT,"PLSDA_medLoadings.csv"), quote=F)
pLoadings = merge(p0$loadings$X,p0$loadings$Y, all=T)
pLoadings$src = c(substr(row.names(p0$loadings$X),3,nchar(row.names(p0$loadings$X))-1),row.names(p0$loadings$Y))
pLoadings$tYpe = c(rep("Ecology",nrow(p0$loadings$X)), rep("Medication",nrow(p0$loadings$Y)))
rEf = data.frame(code=f0, group=mNam)
#rEf = data.frame(code=c("0000","0010","0100","0110","1010","1011","1110","1111"), group=c("No medications", "Drugs/supplements", "Antimicrobials", "Antimicrobials + Drugs/supplements", "CFTR modulators + Drugs/supplements", "CFTR modulators + Drugs/supplements + Drug-drug interaction", "CFTR modulators + Antimicrobials + Drugs/supplements", "CFTR modulators + Antimicrobials + Drugs/supplements + Drug-drug interaction"))
for(i in which(pLoadings$tYpe=="Medication")){pLoadings$src[i] = rEf[which(rEf[,1]==pLoadings$src[i]),2]};rm(i)

pdf(paste0(ptOT,"loadPLSDA.pdf"), width=7, height=7)
#jpeg(paste0(ptOT,"loadPLSDA.jpg"), width=490, height=490)
plot(pLoadings$comp1,pLoadings$comp2, type="p", cex=as.numeric(gsub("Medication",5,gsub("Ecology",2,pLoadings$tYpe))), pch=as.numeric(gsub("Medication",20,gsub("Ecology",18,pLoadings$tYpe))), xlab="PC1 Loading", ylab="PC2 Loading", col=cBp[as.numeric(gsub("Medication",2,gsub("Ecology",3,pLoadings$tYpe)))], cex.axis=1.5, cex.lab = 1.5)
abline(h=0,v=0)
text(pLoadings$comp1,pLoadings$comp2, label=pLoadings$src, cex=.1)
invisible(dev.off())

