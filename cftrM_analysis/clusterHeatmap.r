#!/bin/env Rscript
# author: ph-u
# script: clusterHeatmap.r
# desc: cftrM data ecology comparisons
# in: Rscript clusterHeatmap.r
# out: cftrM_result/heatmap_*.pdf
# arg: 0
# date: 20220822

library(PMCMRplus) # v1.9.3
library(stringr)
library(lattice);library(latticeExtra)
nRep = 7 # number of replicates
pT = "../cftrM_data/"
nAm = list.files(pT,"-eco.csv")
#n = c();for(i in 1:length(nAm)){if(nrow(read.csv(paste0(pT,nAm[i]), header=T, stringsAsFactors=F))<1){n = c(n,i)}};nAm = nAm[-n]
n = data.frame(fName=nAm,gRoup=0,tGroup=0,sTart=0,eNd=0, stringsAsFactors=F)
nAm = str_split(n$fName,"_")
for(i in 1:nrow(n)){n[i,-1] = c(paste0("g",nAm[[i]][2:3]), paste0("20",c(substr(nAm[[i]][3],1,2),substr(nAm[[i]][3],3,4))))}
fQ = table(n$tGroup);fQ = fQ[which(fQ>1)]
ecoType = data.frame(eCo=unique(read.csv(paste0(pT,n[1,1]), header=T, stringsAsFactors=F)$c1_is),cAt=NA)
ecoType$cAt = ifelse(substr(ecoType$eCo,1,4)=="comm","commensalism",ifelse(substr(ecoType$eCo,1,3)=="pre","predation",ifelse(substr(ecoType$eCo,1,4)=="harm","amensalism",ecoType$eCo)))
## statistics bin
i = c(colnames(n)[2:3],"ecoType","replicate","ratio")
statBin = as.data.frame(matrix(nr=0,nc=length(i)))
colnames(statBin) = i

##### time-series ecology heatmap #####
for(i in 1:length(fQ)){
	a = vector(length=fQ[i], mode="list")
	i6 = rep(0,fQ[i]) # total valid simulation sets
	for(i0 in 1:fQ[i]){
		a[[i0]] = read.csv(paste0(pT,n[which(n$tGroup==names(fQ)[i]),1][i0]), header=T, stringsAsFactors=F)
		i7 = a[[i0]]$fit_sim[which(a[[i0]][,1] == a[[i0]][1,1] & a[[i0]][,2] == a[[i0]][1,2] & a[[i0]][,4] == a[[i0]][1,4])]
		i6[i0] = ifelse(length(i7)==0,0,sum(i7))
	}
## plot matrix
	heatMedian = matrix(0,nr=length(unique(ecoType$cAt)),nc=length(a))
	row.names(heatMedian) = unique(ecoType$cAt)
	colnames(heatMedian) = names(a) = n[which(n$tGroup==names(fQ)[i]),2]
## data gropuing
	for(i0 in 1:ncol(heatMedian)){for(i1 in 1:nrow(heatMedian)){ # [ecology, medication]
		i5 = ecoType$eCo[which(ecoType$cAt==row.names(heatMedian)[i1])]
		i2 = rep(0,nRep)
		for(i3 in unique(a[[i0]]$replicate)){
			i4 = a[[i0]][which(a[[i0]]$replicate==i3 & (a[[i0]]$c1_is %in% i5)),]
			if(length(i5)>1){i4$fit_sim[which(i4$c1_is %in% i5[-1])]=0} # avoid double-count possible max (permutation vs Combination)
			i2[i3] = sum(i4$count)/sum(i4$fit_sim) # ratio of specific type against possible max
		}
		statBin[nrow(statBin)+(1:length(i2)),] = data.frame(colnames(heatMedian)[i0],names(fQ)[i],row.names(heatMedian)[i1],1:7,i2)
		heatMedian[i1,i0] = mean(i2) # for plotting
	}}
	cOl = function(p,m=.5){return(gray(rev(seq(m,1,length.out=p))))}
## dendrogram
# https://stackoverflow.com/questions/6673162/reproducing-lattice-dendrogram-graph-with-ggplot2
	dEn = as.dendrogram(hclust(dist(t(heatMedian))))
	rOw = order.dendrogram(dEn)
## grouping matrix
	gMtx = matrix(0,nr=nchar(colnames(heatMedian)[1])-1,nc=ncol(heatMedian))
	colnames(gMtx) = colnames(heatMedian)
	row.names(gMtx) = c("CFTRm", "Antimicrobial", "Chemicals", "CFTRm interacts")
	for(i0 in 1:ncol(gMtx)){gMtx[which(str_split(colnames(gMtx)[i0],"")[[1]]==1)-1,i0] = 1}
	fCap = function(x){paste(toupper(substr(x,1,1)),substr(x,2,nchar(x)),sep="")}
	row.names(heatMedian) = fCap(row.names(heatMedian))
## plot
	pdf(paste0(gsub("data","result",pT),"heatMap_",names(fQ)[i],".pdf"), width=9)
## heatmap
# https://stackoverflow.com/questions/15587734/saving-levelplot-to-file-in-rs-lattice-package
	print(levelplot(t(rbind(gMtx[rev(1:nrow(gMtx)),rOw],heatMedian[,rOw])),
		col.regions=cOl(100),aspect="iso",xlab=list("Number of matching simulations",cex=2),ylab=list("Medication  |  Ecological roles",cex=2),
		scales=list(x=list(at=1:ncol(gMtx),labels=i6[rOw],cex=2, rot=90), y=list(cex=2)), colorkey=list(labels=list(cex=2)), # x=list(cex=2, rot=90)
		panel=function(...){panel.levelplot(...);panel.abline(h=nrow(gMtx)+.5, lwd=4)},
		legend=list(top=list(fun=dendrogramGrob,args=list(x=dEn,side="top",size=4)))))
	invisible(dev.off())
}
## medication grouping
# https://stackoverflow.com/questions/2540129/lattice-multiple-plots-in-one-window
#	print(levelplot(t(gMtx[rev(1:nrow(gMtx)),rOw]),col.regions=cOl(100,0), colorkey=F,aspect="iso", scales=list(x=list(at=NULL),y=list(cex=2)),xlab=list("",cex=0),ylab=list("Medication\ngroup",cex=2)), split=c(1,2,1,2), newpage=FALSE)
#	grid.arrange(pLt1,pLt2, nrow=2)
#	pL = c(1,1,1,1,1,1,2)
#	layout(matrix(pL,nrow=length(pL),ncol=1,byrow=T))
#	par(mar=c(3,0,0,0)+.1,cex=2)
#	heatmap(heatMedian, Rowv=NA, col=cOl(100), RowSideColors=cOl(nrow(heatMedian)),cexRow=2,cexCol=2)
#	mtext("Mean ratio across replicates",side=2,cex=2,adj=.25,padj=7)
#	mtext("Medication group",side=1,cex=2,adj=1,padj=1)
#	mtext("Ecological roles",side=4,cex=2,padj=-1.4)
#	text(rep(.2,2),c(.02,.77),label=round(range(heatMedian),1),cex=2)
