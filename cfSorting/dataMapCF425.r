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
#clREF = read.table("../data/colDict-ExcSep.csv", sep="!", header=T, stringsAsFactors=F)
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
## remap medication history to per patients basis 20220421
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
dRop = c();for(i in c("reason", "culture[dt][ay]", "^is_pt", "_date", "freq$", "positivevalue", "dates$", "reastop", "s05numberof", "cult_[dt]", "result[ap]", "s05sectionstatus", "details$", "cult_sample_[stb]", "_specify_", "e_pulm", "s05ntm[hsnac][ateu]", "_validated", "rea[so][it]", "datesofstopping")){dRop = c(dRop, grep(i,colnames(rEc)))};dRop = unique(dRop)
regID = c();for(i in yEar){regID = unique(c(regID, cf425[[i]]$regid_anon))};rm(i)
rEc = rEc[which(rEc$regid_anon %in% regID),-dRop]
# Only patients with annual review presence in data + question multichoices are "tick" or "not ticked" = unticked entries == negative response instead of no info (i.e. 0 instead of NA)

##### separate direct/indirect reference columns ##### 20220421
mUlti = c("cult_specify","s05culturespeciesresistotherspec","macrolidespec","enzymespecifyother","oralbroncho_bamed","oralbroncho_theophmed","inhbroncho_sabamed","inhbroncho_labamed","inhbroncho_saacmed","inhbroncho_laacmed","inhbroncho_baacmed","other_antibiotics_specify","s05culturespeciesfungalotherspec","s05culturespeciesviralotherspeci","corticocombomed","drugs","corticoinhmed","leukotrienemodmed","antifungalsmed")
inDir = which(colnames(rEc)%in%mUlti)
dirRef = rEc[,-c(1,2,inDir)];indRef = rEc[,c(1,2,inDir)] # sep (in)direct extractable info 20220323

##### standardize data ##### 20220212, 20220315, 20220421
dRop = c(); toZero = c("no","n","nk","unknown")
for(i in 1:ncol(dirRef)){ if(length(grep("/",dirRef[,i]))>0){
	dirRef[,i] = ifelse(is.na(dirRef[,i]),0,1)
}else{
	dirRef[,i] = tolower(dirRef[,i])
	dirRef[is.na(dirRef[,i]),i] = 0
	for(i0 in toZero){dirRef[dirRef[,i]==i0,i] = 0}
	if(any(dirRef[,i]==2)){if(length(unique(dirRef[which(dirRef[,i]>2),i]))>0){
		dirRef[,i] = ifelse(dirRef[,i]!=0,1,0) # counts
	}else{
		dirRef[,i] = ifelse(dirRef[,i]==1,1,0)
		# 1 = checked; 2 = unchecked (sourced cf425 dictionary)
	}}
#	if(length(unique(dirRef[,i]))>2){cat(i,":",length(unique(dirRef[,i]))," - ",colnames(dirRef)[i]," - ",if(length(unique(dirRef[,i]))<10){unique(dirRef[,i])}else{head(unique(dirRef[,i]))},"\n")}
}};rm(i)
dirRef[dirRef!=0] = 1 # restrict to presence/absence data per patient
for(i in 1:ncol(dirRef)){if(length(unique(dirRef[,i]))<2){dRop = c(dRop,i)}} # non-info
dirRef = dirRef[,-dRop] # non-informative columns
dirRef = cbind(indRef[,1:2],dirRef)
colnames(dirRef) = trimws(colnames(dirRef), which="both")

##### extract direct ref column names ##### 20220227, 20220422
colNam = colnames(dirRef)[-c(1:2)]
write.table(c("input",colNam[order(colNam)]),"../data/genCols.csv", sep="\t",quote=F,row.names=F, col.names=F)

save(rEc,yEar,yR,mUlti,dirRef,indRef, file="../data/cf425FULL.rda")
