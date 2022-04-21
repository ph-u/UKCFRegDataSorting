#!/bin/env Rscript
# author: ph-u 
# script: dataMapCF425.r
# desc: cf425 data remapping
# in: Rscript dataMapCF425.r
# out: {otherSp,genCols}.csv, cf425FULL.rda
# arg: 0
# date: 20220210

##### env #####
load("../data/cf425.rda") # 5 sec - cf425, uniqNam, yEar
clREF = read.table("../data/colDict-ExcSep.csv", sep="!", header=T, stringsAsFactors=F)
#clSUM = read.csv("../data/attributesTotal.csv", header=T, stringsAsFactors=F, row.names=1)
yR = as.numeric(substr(names(cf425)[yEar], 1, 4))

##### fuse data ##### 5 mins
#d0 = date()
rEc = cf425[[1]]; rEc$year = yR[1]
for(i in 2:length(yEar)){ #cat(yR[i],", ")
	if(i<7){ ## 2008-2013
		x = cf425[[yEar[i]]]
	}else if(i<9){
		x = merge(cf425[[yEar[i]]], cf425[[yEar[i]+1]], all=T)
	}else{
		p = names(cf425)[grep(yR[i], names(cf425))]
		mEd = grep("medi", p)
		x0 = cf425[[p[mEd[1]]]][,1:5]
		for(i0 in 2:length(mEd)){x0 = rbind(x0,cf425[[p[mEd[i0]]]][,1:5])}
		for(i0 in 2:ncol(x0)){x0[,i0] = tolower(x0[,i0])}
		i1 = c();for(i0 in c("vitamin", "minerals","biphosphanates")){i1 = c(i1,grep(i0,x0[,4]))}
		for(i0 in c("uppl", "sleep", "vitamin", "replacement", "anaphylaxis", "haliborange", "salt", "iron", "nurtini", "sugar", "skin", "epipen", "sea ", "eyedrops", "lotion", "nutri ", "juice", "feed", "energy", "magnesium", "nutrition", "leg ")){i1 = c(i1,grep(i0,x0[,5]))}
		x0 = x0[-unique(i1),]
		x0 = x0[which(!is.na(x0[,2])),]
## remap medication history to per patients basis
		medID = unique(x0$regid_anon)
		x0[grep("please",x0[,2]),2] = x0[which(is.na(x0[,3])),3] = ""
		for(i0 in c("vit", "shake", "juice", "water", "sodium chloride", "factor", "bath oil")){x0[grep(i0,x0[,2]),2] = x0[grep(i0,x0[,3]),3] = ""}
		x0$drugs = gsub("@$","",gsub("^@","",paste(x0[,2],x0[,3], sep="@")))
		x0$drugs[grep("[0-9]+/[0-9]+/",x0$drugs)] = x0[grep("[0-9]+/[0-9]+/",x0$drugs),4]
		x1 = as.data.frame(matrix(NA,nr=length(medID),nc=2))
		colnames(x1) = c(colnames(rEc)[1],"drugs")
		x1[,1] = medID
		for(i0 in 1:nrow(x1)){x1[i0,2] = paste(x0$drugs[which(x0[,1]==x1[i0,1])],collapse="@")}
		x1$drugs = gsub("!","e/a",gsub("/","@",gsub("e/a","!",x1$drugs)))
		for(i0 in c(",","%"," x"," & ","[(]","[)]","e[.]g[.]", " @", "@ ", rep("@@",2))){x1$drugs = gsub(i0,"@",x1$drugs)}
		x1$drugs = gsub("^@","",gsub("@$","",x1$drugs))
## merge reformatted medication histories to remaining data
		x0 = merge(x1, cf425[[p[grep("NTM", p)]]], all=T)
		x = merge(cf425[[p[grep("review", p)]]], x0, all=T)
		rm(x0, x1, p, mEd)
	};x$year = yR[i]
	colnames(x) = tolower(colnames(x))
	rEc = merge(rEc,x, all=T)
}#;cat("\n");d0;date()

##### extract annual review patients & useful data attributes ##### 20220212, 20220422
dRop = c();for(i in c("reason", "culturedate", "^is_pt", "_date", "freq$", "positivevalue", "dates$", "culturetype", "reastop")){dRop = c(dRop, grep(i,colnames(rEc)))};dRop = unique(dRop)
regID = c();for(i in yEar){regID = unique(c(regID, cf425[[i]]$regid_anon))};rm(i)
rEc = rEc[which(rEc$regid_anon %in% regID),-dRop]
# Only patients with annual review presence in data + question multichoices are "tick" or "not ticked" = unticked entries == negative response instead of no info (i.e. 0 instead of NA)

##### standardize data ##### 20220212, 20220315
for(i in 3:ncol(rEc)){rEc[,i] = tolower(rEc[,i])
	uQ = unique(rEc[,i]);uQ0 = c(NA,1,2)
	if(length(setdiff(uQ,uQ0))==length(setdiff(uQ0,uQ))){rEc[,i] = ifelse(rEc[,i]==1,1,0)}
};rm(i,uQ)
rEc[is.na(rEc)] = rEc[rEc=="no"] = rEc[rEc=="n"] = rEc[rEc=="nk"] = rEc[rEc=="unknown"] = 0
mUlti = c("cult_specify","s05culturespeciesresistotherspec","macrolidespec","enzymespecifyother","oralbroncho_bamed","oralbroncho_theophmed","inhbroncho_sabamed","inhbroncho_labamed","inhbroncho_saacmed","inhbroncho_laacmed","inhbroncho_baacmed","other_antibiotics_specify","s05culturespeciesfungalotherspec","s05culturespeciesviralotherspeci","corticocombomed","drugs")
inDir = which(colnames(rEc)%in%mUlti)
dirRef = rEc[,-c(1,2,inDir)];indRef = rEc[,c(1,2,inDir)] # sep (in)direct extractable info 20220323
dirRef[dirRef!=0] = 1 # restrict to presence/absence data per patient
dirRef = cbind(indRef[,1:2],dirRef)
colnames(dirRef) = trimws(colnames(dirRef), which="both")

##### extract column (medicine + standard species) names ##### 20220227
dEtails = which(colnames(indRef) %in% mUlti[c(3:12,15:length(mUlti))])
zZ = c();for(i in dEtails){zZ = c(zZ,unique(indRef[,i]))}
zZ = unique(zZ[-which(zZ=="0" | zZ=="nil" | zZ=="not taking")])
zZ1 = setdiff(grep(", ",zZ), grep(")",zZ))
zZ0 = read.table(text=zZ[zZ1], sep=",")
zZ = zZ[-zZ1]
for(i in 1:ncol(zZ0)){zZ = c(zZ, trimws(zZ0[,i], which="both"))}
zZ = zZ[-grep("^[[:digit:]]",zZ)]
zZ = c(zZ,colnames(dirRef)[-c(1,2)])
zZ = unique(zZ)
write.table(c("input",zZ[order(zZ)]),"../data/genCols.csv", sep="\t",quote=F,row.names=F, col.names=F)

save(rEc,yEar,yR,mUlti,dirRef,indRef, file="../data/cf425FULL.rda")
