#!/bin/env Rscript
# author: ph-u
# script: F508delTSPlots.r
# desc: rolling mean plots of 3 consecutive years for F508del mutation CFTRm segregated time-series
# in: Rscript F508delTSPlots.r
# out: graph/F508del/*.pdf
# arg: 0
# date: 20221204

##### env #####
nRep = 7
simO = 500
mTx = .3
ptIN = "../F508del_data/"
ptOT = "../graph/F508del/"
sEq = c(2,9,3,4,5,6,7,1,8)
cBp = palette.colors(palette = "Okabe-Ito", alpha=1, recycle = T)[sEq]
#cBl = palette.colors(palette = "Okabe-Ito", alpha=.1, recycle = T)[sEq]

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
samSize = read.csv("../graph/cftrm_samFreq.csv", header=T)
x = c("Year","total","CFTRm","antimicrobials","chemicals","drug_interaction")
tReat = as.data.frame(matrix(nr=nrow(samSize), nc=length(x)))
colnames(tReat) = x;rm(x)
tReat$Year = samSize$year
tReat$total = apply(samSize[,-1],1,sum)
for(i in 3:ncol(tReat)){
	tReat[,i] = apply(samSize[,which(substr(colnames(samSize),i-1,i-1)==1)],1,sum)
};rm(i)

pdf(paste0(ptOT,"medication.pdf"), width=14, height=11)
par(mar=c(5,4.5,1,1)+.1, mfrow=c(2,1), xpd=T)
matplot(tReat[,1], tReat[,-1], type="l", lty=1:(ncol(tReat)-1), col=1:(ncol(tReat)-1), lwd=3, xlab="Year", ylab="Samples", cex.axis=1.7, cex.lab=2.1)
plot.new()
lPt = legPlot(colnames(tReat)[-1], 2)
legend("top", legend=sub("_"," ",capFirst(colnames(tReat)[-1]))[lPt[[1]]], border=NA, ncol=lPt[[2]], lty=c(1:(ncol(tReat)-1),NA), title="Annual Review Medication Records", col=1:(ncol(tReat)-1), lwd=5, cex=2.1)
invisible(dev.off())

##### files #####
cat("Import file list: ",date(),"\n")
f = list.files(ptIN,"-eco")
f0 = unlist(strsplit(f,"_"))
fNam = as.data.frame(matrix(nr=length(f),nc=2))
for(i in 1:nrow(fNam)){fNam[i,] = f0[(2:3)+length(f0)/length(f)*(i-1)]};rm(i,f0)
fNam$start = substr(fNam[,2],1,2); fNam$end = substr(fNam[,2],3,4)
f0 = unique(fNam[,1])
yR = as.numeric(unique(c(fNam$start,fNam$end)))+2000;yR = yR[order(yR)]

##### data bin to summary rolling-time graphs #####
cat("Processing group:\n")
for(i in 1:length(f0)){z = 0; cat(i,"(",date(),"),\n")
	for(i0 in which(fNam[,1]==f0[i])){
		f1 = read.csv(paste0(ptIN,f[i0]), header=T, stringsAsFactors=F)
		t0 = as.numeric(fNam$end[i0])+2000
## set-up record data frame
		if(z==0){
			if(i==1){
				sPair = unique(paste0(f1$category1,"_",f1$category2))
				eCo = unique(f1$c1_is)
			}
			ecoTS = vector("list", length(sPair)) # list of pairwise time-series
			for(i1 in 1:length(ecoTS)){
				ecoTS[[i1]] = as.data.frame(matrix(0,nr=length(yR)-2, nc=length(eCo)+1))
				colnames(ecoTS[[i1]]) = c("year",eCo)
				ecoTS[[i1]]$year = yR[-c(1:2)]
			};rm(i1)
		z = 1}
## segregate data (one end year) into dataframes
		for(i1 in 1:length(sPair)){
			a0 = strsplit(sPair[i1],"_")[[1]]
			f2 = f1[which(f1$category1==a0[1] & f1$category2==a0[2]),]
			for(i2 in eCo){ecoTS[[i1]][which(ecoTS[[i1]]$year==t0),which(colnames(ecoTS[[i1]])==i2)] = sum(f2$count[which(f2$c1_is==i2)])} # count one eco relationship
		};rm(i1)
	};rm(i0)
## plot dataframes
	for(i0 in 1:length(sPair)){
#		for(i1 in 3:ncol(ecoTS[[i0]])){ecoTS[[i0]][,i1] = ecoTS[[i0]][,i1] + ecoTS[[i0]][,i1-1]};rm(i1)
		ecoTS[[i0]][,-1] = ecoTS[[i0]][,-1]/(simO*nRep)*100
		tS = c(ecoTS[[i0]]$year,rev(ecoTS[[i0]]$year))
### stacked barplot prep
		ecoPlot = ecoTS[[i0]][,-1]
		ecoPlot[nrow(ecoPlot)+c(1:2),] = 0
		ecoPlot = as.matrix(t(ecoPlot[c(nrow(ecoPlot)-1:0,1:(nrow(ecoPlot)-2)),]))
		colnames(ecoPlot) = paste0(yR,"\n(",samSize[,which(colnames(samSize)==paste0("g",f0[i]))],")\n(",tReat[,2],")")

### export
		pdf(paste0(ptOT,f0[i],"_",sPair[i0],".pdf"), width=14, height=11)
		par(mar=c(5,4.5,2,1)+.1, mfrow=c(2,1), cex.axis=1.4, xpd=T)
		barplot(ecoPlot, ylim=c(0,100), xaxt="n", ylab=paste(sPair[i0], "(%)"), xlab="", col=cBp, border="white", cex.axis=2.1, cex.lab=1.5)
		axis(1, at=-.5+1.2*(1:length(yR)), padj=.7, labels=colnames(ecoPlot))
		mtext("Year (Group Sample Size) (Total Sample Size)",side=1,padj=4.9,cex=2.1)
### legend plot
		plot.new()
		lPt = legPlot(eCo)
		legend("top", inset=c(0,0), legend = capFirst(gsub("_"," ",sub("c2","B",eCo[lPt[[1]]]))), title=paste("Ecological Relationship - 100% =",nRep*simO,"simulations"), border=NA, xpd=T, cex=2, ncol=lPt[[2]], pch = rep(19,length(eCo)), col = cBp)
		invisible(dev.off())
	};rm(i0)
};rm(i);cat("Done",date(),"\n")
