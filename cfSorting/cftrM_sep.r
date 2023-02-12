#!/bin/env Rscript
# author: ph-u
# script: cftrM_sep.r
# desc: categorize people with CF (pwCF) into CFTR-modifiers (cftrM) and interacting drugs categories + separate into different time-series
# in: Rscript cftrM_sep.r
# out: ../cftrM_raw/cftrM_{abs,mod,int}_[start][end]_gLV.csv
# arg: 0
# date: 20220814

##### import #####
cat("data reading, environment preparation:",date(),"\n")
load("../data/cf425MedMic.rda");
ptOT = "../p_raw_20230207/"
#library(stringr);
muT = read.csv("../data/mutPWCF.csv", header=T, stringsAsFactors=F);

##### generate list of filtered drugs 20230207 #####
#rEf = read.csv("../data/cftrm_interaction.csv", header=T,stringsAsFactors=F);
#cftrM = c();for(i in strsplit(unique(rEf$cftrm),";")){cftrM = c(cftrM,i)};rm(i);cftrM = tolower(unique(cftrM)); # caftors list

#anF = read.csv("../data/antiInfectives.csv", header=T, stringsAsFactors=F);
#antiF = c();for(i in 1:3){for(i0 in strsplit(unique(anF[,i]),";")){antiF = c(antiF,strsplit(i0," / ")[[1]])}};rm(i,i0);antiF = unique(tolower(antiF));antiF = unique(c(antiF,sub("s$","",antiF[grep("s$",antiF)]))) # anti-infectives list: brand, ingredients, class
#lsDrugs = data.frame(nAm=colnames(mEdic)[-c(1:2)],src="")
#lsDrugs[grep("lactobac", lsDrugs[,1]),2] = "probiotics"
#for(i in c("zinc","sodium","potassium","calc","iron","magnesium")){lsDrugs[grep(i, lsDrugs[,1]),2] = i};rm(i)
#lsDrugs[grep("mycobacterium", lsDrugs[,1]),2] = NA
#lsDrugs[,2] = ifelse(lsDrugs[,1] %in% antiF, "antimicrobials",ifelse(lsDrugs[,1] %in% cftrM, "CFTRm",lsDrugs[,2]))
#write.csv(lsDrugs, paste0(ptOT,"F508dd_drugs.csv"), quote=F, row.names=F)

##### remap medication history 20230210 #####
dRug = read.csv("../data/F508dd_drugs.csv", header=T)
d0 = unique(unlist(strsplit(unique(dRug$src), split = ";")), na.rm=T)
d0 = d0[order(d0)];d0 = d0[!is.na(d0)] # 210 medication groups
mEd0 = as.data.frame(matrix(0,nr=nrow(mEdic),nc=length(d0)+2))
colnames(mEd0) = c(colnames(mEdic)[1:2],d0)
for(i in which(!is.na(dRug[,2]))){
	d1 = unlist(strsplit(dRug[i,2],";"))
	mEd0[,d1] = mEd0[,d1] + mEdic[,dRug[i,1]]};rm(i,d1)
mEd0 = as.data.frame(mEd0>0)
mEd0[,1:2] = mEdic[,1:2]

##### env #####
lIm = .05; # presence threshold of a category to be considered as a separate one
yR = unique(mIcro$year);yR = yR[order(yR)]; # timeline
cOl = palette.colors(n =  length(yR), palette = "Okabe-Ito", .5, recycle=T)
#rEf = rEf[!is.na(rEf$affectedDrug),-ncol(rEf)]; # filter only confirmed interacting drugs
#tricaftor = cftrM[grep("tri",cftrM)] # tricaftor list
#intDr = c();for(i in strsplit(unique(rEf$drug),";")){intDr = c(intDr,i)};rm(i);intDr = tolower(unique(intDr)); # interacting drug list
#noMed = c("no_medication","exacerbations","aqua","starch","water",colnames(mEdic)[grep("intravenous",colnames(mEdic))]) # colnames representing no drugs
#dRugs = colnames(mEdic)[which(!(colnames(mEdic) %in% c(noMed,antiF,cftrM,"regid_anon","year")))]
# assume if data only contain drug class, the drug used is not interacting with cftrM

##### get only F508del pwCF #####
cat("get the F508dd pwCF group:",date(),"\n");
F508del = muT$regid_anon[which(muT[,3]==1 & rowSums(muT[,-c(1:2)])==1)] #F508 homozygote (6419/12727)
mIcro = mIcro[which(mIcro$regid_anon %in% F508del),] # 59422/112980
#mEdic = mEdic[which(mEdic$regid_anon %in% F508del),]
mEd0 = mEd0[which(mEd0$regid_anon %in% F508del),]

##### sort medication categories #####
cat("sort medication categories:",date(),"\n");
medPc = as.data.frame(matrix(0,nr=length(yR), nc=ncol(mEd0)))
colnames(medPc) = c(colnames(mEd0)[2],"total",colnames(mEd0)[-c(1:2)]); medPc[,1] = yR
for(i in 1:nrow(medPc)){
	d1 = mEd0[which(mEd0[,2]==medPc[i,1]),-c(1:2)]
	medPc[i,-c(1:2)] = colSums(d1)
	medPc[i,2] = nrow(d1)};rm(i,d1)

## plot medication distribution ratio
mP0 = t(medPc[,-c(1:2)]/medPc[,2]*100)
colnames(mP0) = yR
d1 = apply(mP0,1,function(s){return(any(s>50))})
write.csv(medPc,paste0(ptOT,"medDist.csv"), quote=F, row.names=F)
pdf(paste0(ptOT,"medDist.pdf"), width=21)
par(mar=c(5,5,1,0)+.1)
matplot(mP0, type="l", xlab = "Medication groups (alphabetical order)", ylab = "Presence in F508del homozygous (%)", ylim = c(-28,100), lty=1:length(yR), col=cOl, cex.axis = 2, cex.lab=2); abline(h=50)
text(which(d1>0),-10, labels=names(which(d1>0)),srt=90)
legend("topleft", legend=yR, box.col="#00000000", lty=1:length(yR), lwd=7, col=cOl, horiz=T, cex=1.4)
invisible(dev.off())

## simple medication categories
mD0 = c(d0[c(grep("CFTR",d0), grep("antim",d0), grep("panc",d0))],"others","grouping")
#mD0 = c("CFTRm","anti-infective","chemical","interaction","grouping");
mEdic0 = as.data.frame(matrix(0,nr=nrow(mEd0),nc=length(mD0)));
colnames(mEdic0) = mD0;
for(i in 1:3){mEdic0[,i] = as.numeric(mEd0[,mD0[i]])};rm(i)
#mEdic0[,1] = ifelse(rowSums(mEdic[,which(colnames(mEdic) %in% cftrM)])>0,1,0); # 0 = no CFTRm
#mEdic0[,2] = ifelse(rowSums(mEdic[,which(colnames(mEdic) %in% antiF)])>0,1,0); # 0 = no anti-infectives
#mEdic0[,3] = ifelse(rowSums(mEdic[,which(colnames(mEdic) %in% dRugs)])>0,1,0); # 0 = no chemicals (drugs, supplements)
#mEdic0[,4] = ifelse(rowSums(mEdic[,which(colnames(mEdic) %in% intDr)])>0 & mEdic0[,1]==1,1,0); # no CFTRm-prescription interactions
mEdic0[,4] = as.numeric(rowSums(mEd0[-c(1:2,which(colnames(mEd0)%in%mD0[1:3]))])>0)
mEdic0[,5] = paste0(mEdic0[,1],mEdic0[,2],mEdic0[,3],mEdic0[,4],sep=""); # grouping

##### genus categories sorting #####
cat("sort genus categories");
#gEn = word(colnames(mIcro)[-c(1,2)],1);
gEn0 = strsplit(colnames(mIcro)[-c(1,2)]," ");gEn = rep(NA,length(gEn0));
for(i in 1:length(gEn0)){gEn[i] = gEn0[[i]][1]};rm(i);
gEn0 = unique(gEn);
cat(" (",length(gEn0),"categories ):",date(),"\n")
yrFq = table(mIcro[,2]);
genProp = as.data.frame(matrix(0,nr=length(yR),nc=length(gEn0)+1));
colnames(genProp) = c("year",gEn0);
genProp$year = yR;
for(i in 1:nrow(genProp)){
	m0 = mIcro[which(mIcro$year==genProp[i,1]),];
	for(i0 in gEn0){genProp[i,i0] = sum(rowSums(cbind(m0[,which(gEn==i0)+2],rep(0,nrow(m0))))>0)} # count pwCF hosting genus
};rm(i,m0,i0);
for(i in 1:length(yrFq)){genProp[i,-1] = genProp[i,-1]>=yrFq[i]*lIm}
genProp = colSums(genProp[,-1])>0;
mIcro0 = as.data.frame(matrix(0,nr=nrow(mIcro),nc=sum(genProp)+3));
colnames(mIcro0) = c(colnames(mIcro)[1:2],"others",gEn0[which(genProp>0)]);
mIcro0[,1:2] = mIcro[,1:2];
# data.frame(spp=colnames(mIcro)[-c(1,2)],genus=gEn,map=ifelse(gEn %in% gEn0[which(genProp>0)],gEn,"others")) # listing column categorisation
gEn = ifelse(gEn %in% gEn0[which(genProp>0)],gEn,"others");
for(i in 3:ncol(mIcro0)){mIcro0[,i] = rowSums(cbind(mIcro[,2+which(gEn==colnames(mIcro0)[i])],rep(0,nrow(mIcro))))>0} # category presence T/F
mIcro0$Mycobacteria = ((mIcro0$Mycobacterium + mIcro0$Mycobacteroides)>0)
mIcro0$Mycobacterium = mIcro0$Mycobacteroides = NULL

##### medication category sorting #####
cat("sort genus per pwCF per year into medication groups:",date(),"\n");
## eliminate "no medication record" entries
mIcro0 = mIcro0[which(mEdic0$grouping!="0000"),]
mEdic0 = mEdic0[which(mEdic0$grouping!="0000"),]
## medication grouping
a = vector(length=length(unique(mEdic0$grouping)), mode="list");
names(a) = unique(mEdic0$grouping)[order(unique(mEdic0$grouping))];
for(i in 1:length(a)){a[[i]] = mIcro0[which(mEdic0$grouping==names(a)[i]),]};rm(i);

##### sample size per year output #####
cat("sample sizes counting:",date(),"\n")
samSize = as.data.frame(matrix(0,nr=length(yR),nc=length(a)+1));
colnames(samSize) = c("year",paste("g",names(a),sep=""));
samSize$year = yR
for(i in 2:ncol(samSize)){a0 = table(a[[i-1]][,2]);samSize[which(samSize[,1] %in% names(a0)),i] = a0}
## graph + counts
#write.csv(samSize,"../data/cftrm_samFreq.csv", row.names=F, quote=F);
write.csv(samSize,paste0(ptOT,"F508dd_samFreq.csv"), row.names=F, quote=F);
#pdf("../../thesis/fig/cftrm_samFreq.pdf", width=14);
#pdf("../graph/cftrm_samFreq.pdf", width=14);
pdf(paste0(ptOT,"F508dd_samFreq.pdf"), width=14);
par(mar=c(4,4,0,5)+.1, cex=1.5, xpd=T);
cOl[10:length(cOl)] = paste0(substr(cOl[10:length(cOl)],1,7),"FF",sep="")
#cOl = grey(seq(.1,.9,length=length(a))) # paste0("#000000",c("00","22","44"));
#cOl = cOl[c(seq(1,length(cOl),2),seq(2,length(cOl),2))]
p = t(samSize[,-1]);colnames(p) = samSize[,1];
yLim = max(rowSums(samSize[,-1]));
barplot(as.matrix(p), col=cOl, xlab="Year", ylab="Number of UK pwCF Annual Reviews", ylim=c(0,signif(yLim,2)+10^(round(log10(yLim)-1))));
legend("bottomright", inset=c(-.09,0),legend=rev(names(a)), fill=rev(cOl), box.col="#00000000", horiz=F, title="pwCF\nmedication\ngroups");
invisible(dev.off());

##### filter years with sufficient sample sizes #####
for(i in 2:ncol(samSize)){samSize[,i] = ifelse(samSize[,i]<5,0,samSize[,i])};rm(i);
sEp = which(colSums(samSize[,-1])>30);
for(i in sEp){cat(i,":",nrow(a[[i]])); a[[i]] = a[[i]][which(a[[i]][,2] %in% samSize[which(samSize[,i+1]>0),1]),];cat("->",nrow(a[[i]]),"\n")}

##### f: category presence proportion by year #####
cat("category presence fraction by year:",date(),"\n")
cProp = function(x,timeCol=2){
	t = unique(x[,timeCol]);
	x0 = as.data.frame(matrix(0,nr=length(t),nc=ncol(x)-1));
	colnames(x0) = colnames(x)[-1];
	x0[,1] = t[order(t)];
	for(i in 1:nrow(x0)){
		x1 = x[which(x[,timeCol]==x0[i,1]),-c(1:2)];
		x0[i,-1] = colSums(x1)/nrow(x1);
	}
	return(x0);
}
a0 = vector(length=length(sEp), mode="list");
names(a0) = names(a)[sEp];
for(i in 1:length(sEp)){a0[[i]] = cProp(a[[sEp[i]]])}
#noCFTRM0 = cProp(noCFTRM); olCFTRM0 = cProp(olCFTRM); itCFTRM0 = cProp(itCFTRM); # category presence proportion per year

##### plot overall percentage plots of medication groups #####
#pdf("../../thesis/fig/cftrm_perc.pdf", width=14);
#pdf("../graph/cftrm_perc.pdf", width=14);
pdf(paste0(ptOT,"F508dd_perc.pdf"), width=14);
par(mar=c(2,4,2,10)+.1, mfrow=c(3,3), cex=.7, xpd=T);
cOl = palette.colors(n =  length(a0), palette = "Okabe-Ito", .5, recycle=T);
for(i in 1:length(a0)){
	matplot(a0[[i]][,1],a0[[i]][,-1],type="b",pch=1:(ncol(a0[[i]])-1),xlim=range(yR),ylim=c(0,1),col=cOl,xlab="",ylab="pwCF ratio with genus presence",main=paste("Medication category:",names(a0)[i]))
	if(i==3){legend("topright", inset=c(-.5,0), legend=colnames(mIcro0)[-c(1:2)], pch=1:(ncol(mIcro0)-1), lty=1:(ncol(mIcro0)-1), col=cOl, box.col="#00000000", horiz=F, title="Microbial\ncategories");
}}
invisible(dev.off());

##### f: cut time-series #####
cat("cut data into different time-series:",date(),"\n")
tsCut = function(x,tYpe){ for(i0 in 1:(nrow(x)-2)){ for(i1 in 3:nrow(x)){ if(x[i1,1]-x[i0,1]==2){#>1){ # 3 consecutive years only
	#write.csv(x[i0:i1,],paste0("../cftrM_raw_20221202/cftrM_",tYpe,"_",paste0(substr(x[c(i0,i1),1],3,4), collapse=""),"_gLV.csv"), row.names=F, quote=F)
	write.csv(x[i0:i1,],paste0(ptOT,"ts/F508dd_",tYpe,"_",paste0(substr(x[c(i0,i1),1],3,4), collapse=""),"_gLV.csv"), row.names=F, quote=F)
}}}}
for(i in 1:length(a0)){tsCut(a0[[i]],names(a0)[i])}
#tsCut(noCFTRM0,"abs"); tsCut(olCFTRM0,"mod"); tsCut(itCFTRM0,"int") # export time-series of various time periods
cat("data preparation completed:",date(),"\n")
