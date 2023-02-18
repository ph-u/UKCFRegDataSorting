#!/bin/env Rscript
# author: ph-u
# script: F508ddEcoCategory.r
# desc: rolling mean plots of 3 consecutive years for F508del mutation CFTRm segregated time-series (single category overall)
# in: Rscript F508ddEcoCategory.r [medication code] [microbial category]
# out: graph/F508del/*.pdf
# arg: 0
# date: 20230216

#argv=(commandArgs(T))
#argv = c("0111","Pseudomonas")
##### env #####
nRep = 7
simO = 500
mTx = .3
ptIN = "../p_raw_20230207/data/" #"../F508del_data/"
ptOT = "../p_raw_20230207/res/" #"../graph/F508del/"
sEq = c(2,9,3,4,5,6,7,1,8)
cBp = palette.colors(palette = "Okabe-Ito", alpha=1, recycle = T)[sEq]
#cBl = palette.colors(palette = "Okabe-Ito", alpha=.1, recycle = T)[sEq]
cBp = c(cBp,paste0(substr(cBp,1,nchar(cBp[1])-2),"aa"))

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

##### files #####
f = list.files(ptIN,"-eco"); f = f[-grep("1101",f)]
fNam = strsplit(f,"_"); fNam = matrix(unlist(fNam), nr=length(grep("_",strsplit(f[1],"")[[1]]))+1)
medGp = unique(fNam[2,])
spNam = unique(read.csv(paste0(ptIN,f[1]), header=T)[,2])
yR = unique(as.numeric(c(substr(fNam[3,],1,2),substr(fNam[3,],3,4)))+2000)
totSimu = nRep*simO*length(spNam)#(length(spNam)*2-1)

cat("Summarizing microbial taxa in each medication group:",date(),"\n")
for(i0 in 1:length(medGp)){ cat(medGp[i0],","); for(i1 in 1:length(spNam)){
	f0 = f[grep(medGp[i0],f)]
	for(i in 1:length(f0)){
		a0 = read.csv(paste0(ptIN,f0[i]), header=T)[c(1,2,4,5)]
		if(nrow(a0)>0){
			a0 = a0[which(a0$category1==spNam[i1] | a0$category2==spNam[i1]),]
			a1 = unlist(strsplit(f0[i], "_"))
			a0$start = as.numeric(substr(a1[3],1,2))+2000
			a0$end = as.numeric(substr(a1[3],3,4))+2000
			if(i==1){
				a = as.data.frame(matrix(nr=0, nc=ncol(a0)))
				colnames(a) = colnames(a0)
			}; a = rbind(a,a0)}; if(sum(a$count)>=(nRep*simO*length(spNam)*i)){break}
	};#rm(i,a0,a1)

##### set subject category #####
	eCo = data.frame(src=unique(a$c1_is), convert=unique(a$c1_is)[c(1,4,7,2,5,8,3,6,9)])
	for(i in which(a$category1!=spNam[i1] & a$category2==spNam[i1])){
		a[i,1:2] = a[i,c(2,1)]
		a$c1_is[i] = eCo[which(eCo$src==a$c1_is[i]),2]
	};rm(i)
	a$ratio = a$count/totSimu

##### summarize ecology #####
	a0 = unique(a[,c(5,6,3)])
	a0$ratio = 0
	for(i in 1:nrow(a0)){
		a0$ratio[i] = sum(a$ratio[which(a[,6]==a0[i,2] & a[,3]==a0[i,3])])
	};rm(i)
	ecoPlot = matrix(0,nr=nrow(eCo), nc=length(yR))
	colnames(ecoPlot) = yR; row.names(ecoPlot) = eCo$src
	for(i in 1:nrow(a0)){ecoPlot[a0$c1_is[i],which(colnames(ecoPlot)==a0$end[i])] = a0$ratio[i]};rm(i)
## collect overall record from taxonomic categories perspective
	a0$medication = medGp[i0]
	a0$subject = spNam[i1]
	if(i0==1 & i1==1){a0Spp = a0}else{a0Spp = rbind(a0Spp,a0)}

##### plot #####
	colnames(ecoPlot) = paste0(colnames(ecoPlot),"\n(",read.csv(paste0(ptOT,"../F508dd_samFreq.csv"), header=T)[,paste0("g",medGp[i0])],")\n[",read.csv(paste0(ptOT,"medication.csv"), header=T)[,2],"]")
	pdf(paste0(ptOT,"sp_",medGp[i0],"_",spNam[i1],".pdf"), width=14, height=11)
	par(mar=c(5,4.5,2,1)+.1, mfrow=c(2,1), cex.axis=1.4, xpd=T)
	barplot(ecoPlot*100, ylim=c(0,100), xaxt="n", ylab=paste(spNam[i1], "(%)"), xlab="", col=cBp, border="white", cex.axis=2.1, cex.lab=1.5)
	axis(1, at=-.5+1.2*(1:length(yR)), padj=.7, labels=colnames(ecoPlot))
	mtext("Year (Group Sample Size) [Total Sample Size]",side=1,padj=4.9,cex=2.1)

	plot.new()
	lPt = legPlot(eCo$src, nDim=length(eCo$src))
	legend("top", inset=c(0,0), legend = capFirst(sub(" of","",sub(" c2","",gsub("_", " ", eCo$src)))), title=paste("Ecological Relationship - 100% =",totSimu,"simulations"), border=NA, xpd=T, cex=2, ncol=lPt[[2]], pch = rep(19,length(eCo)), col = cBp)
	invisible(dev.off())
}};cat("\n")

##### PCA: ecology vector, taxa colour #####
cat("PCA (ecology vector, taxa colour):",date(),"\n")
library(ggbiplot)
eCo$ecoType = c("Mutualism", "Commensalism", "Predatory/Parasitism", "Commensalism", "Neutral/No interaction", "Amensalism", "Predatory/Parasitism", "Amensalism", "Competition")
p0 = unique(a0Spp[,5:6])
p0 = cbind(p0,as.data.frame(matrix(0,nr=nrow(p0),nc=length(unique(eCo$ecoType))+1)))
colnames(p0)[-c(1:2)] = paste0("<-",c(unique(eCo$ecoType),"Unpredictable"),"->")
for(i in 1:nrow(a0Spp)){
	corD = c(which(p0$subject==a0Spp$subject[i] & p0$medication==a0Spp$medication[i]),grep(eCo$ecoType[which(a0Spp$c1_is[i]==eCo$src)],colnames(p0)))
	p0[corD[1],corD[2]] = p0[corD[1],corD[2]] + a0Spp$ratio[i]/(length(yR)-2)
};rm(i)
# summarize medication group & year
p0[,ncol(p0)] = 1-rowSums(p0[,3:(ncol(p0)-1)])
row.names(p0) = 1:nrow(p0)

p = prcomp(p0[,-c(1:2)], scale.=T)
pCa0 = ggbiplot(p, var.scale = 1,
        groups=p0[,2], ellipse = F,#TRUE, ellipse.prob = .95,
        labels = 1:nrow(p0), labels.size = 4,
        varname.size = 2.8, varname.adjust = c(4.1,2,3.9,1.5,2.5,2,1.8) # mut,comm,pred,neut,ame,comp,unp
)+#xlab("")+ylab("")+
scale_color_manual(values=setNames(cBp[c(5,1:3,6:10,4)], unique(p0[,2])))+
guides(color = guide_legend(title="Taxonomic category"))+theme_bw()+#ylim(c(-4,9))+xlim(c(-20,20))+
scale_y_continuous(breaks = seq(-6, 4, 2), limits = c(-6, 4))+
coord_cartesian(xlim=c(-6,12))+
theme(legend.position = 'bottom', legend.direction = "vertical",
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        legend.title = element_text(size=18),
        axis.text=element_text(size=16),axis.title=element_text(size=14),
        plot.margin=margin(t=0,r=0,b=0,l=1))
ggsave(paste0(ptOT,"sp_PCA.pdf"), plot=pCa0, width=6, height=7) # last_plot()

##### PLS-DA: ecology vector, taxa colour #####
cat("PLS-DA (ecology vector, taxa colour):",date(),"\n")
library(mixOmics)
p1 = plsda(p0[,-c(1:2)], as.factor(p0[,2]))
pdf(paste0(ptOT,"sp_PLSDA.pdf"), width=14, height=14)
plotIndiv(p1, ellipse=T, ellipse.level=.95, col=cBp[1:length(spNam)], size.xlabel=rel(1.4), size.ylabel=rel(4), size.axis=rel(2), legend=T, style="graphics", cex=3) # cBp[c(3,4,7,8,5,6,2,1)]
invisible(dev.off())

##### PCA: medication vector, taxa colour (illogical, `.` both independent var) #####
cat("PCA (medication vector, taxa colour):",date(),"\n")
p2 = unique(a0Spp[,c(2,6)])
p2 = cbind(p2,as.data.frame(matrix(0,nr=nrow(p2),nc=length(unique(medGp)))))
colnames(p2)[-c(1:2)] = paste0("<-",medGp,"->")
for(i in 1:nrow(a0Spp)){
	corD = c(which(p2$end==a0Spp$end[i] & p2$subject==a0Spp$subject[i]),which(colnames(p2)==paste0("<-",a0Spp$medication[i],"->")))
	p2[corD[1],corD[2]] = p2[corD[1],corD[2]] + a0Spp$ratio[i]/(length(yR)-2)
};rm(i)
p = prcomp(p2[,-c(1:2)], scale.=T)
pCa0 = ggbiplot(p, var.scale = 1,
        groups=p2[,2], ellipse = TRUE, ellipse.prob = .95,
        labels = 1:nrow(p2), labels.size = 4,
        varname.size = 2.8#, varname.adjust = c(3,2,3,1.5,2.5,2,1.8) # mut,comm,pred,neut,ame,comp,unp
)+#xlab("")+ylab("")+
scale_color_manual(values=setNames(cBp[c(5,1:3,6:10,4)], unique(p0[,2])))+
guides(color = guide_legend(title="Taxonomic category"))+theme_bw()+#ylim(c(-4,9))+xlim(c(-20,20))+
#scale_y_continuous(breaks = seq(-6, 4, 2), limits = c(-6, 4))+
coord_cartesian(xlim=c(-2,2))+
theme(legend.position = 'bottom', legend.direction = "vertical",
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        legend.title = element_text(size=18),
        axis.text=element_text(size=16),axis.title=element_text(size=14),
        plot.margin=margin(t=0,r=0,b=0,l=1))
ggsave(paste0(ptOT,"med_PCA.pdf"), plot=pCa0, width=6, height=7) # last_plot()

##### cladistics: ecology base, medication tag #####
cat("cladistics (ecology base, medication tag):",date(),"\n")
row.names(p0) = paste0("(",row.names(p0),") ",p0[,2])
#row.names(p0) = paste0("(",row.names(p0),") ",p0[,1],"-",p0[,2])
hc = hclust(dist(p0[,-c(1:2)], method = "euclidean"), method = "centroid");
pdf(paste0(ptOT,"med_clus.pdf"), width=14, height=14);
plot(hc, lwd=3);
#rect.hclust(hc, k = 8, border=cBp[7]);
invisible(dev.off());
