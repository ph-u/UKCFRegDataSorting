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

##### fuse data ##### ~6 mins
#d0 = date()
rEc = cf425[[1]]; rEc$year = yR[1]
for(i in 2:length(yEar)){ #cat(yR[i],", ")
	if(i<7){ ## 2008-2013
		x = cf425[[yEar[i]]]
	}else if(i==7 | i==8){
		x = merge(cf425[[yEar[i]]], cf425[[yEar[i]+1]], all=T)
	}else{
		p = names(cf425)[grep(yR[i], names(cf425))]
		mEd = grep("medi", p)
		x0 = cf425[[p[mEd[1]]]]
		for(i0 in 2:length(mEd)){x0 = rbind(x0,cf425[[p[mEd[i0]]]])}
		x0 = x0[which(!is.na(x0$drugname)),]
## remap medication history to per patients basis # 2 mins
		medID = unique(x0$regid_anon)
		#drugID = unique(as.character(apply(x0[,2:3],1, function(x){gsub("-NA","",paste(x, collapse="-"))})))
		drugID = unique(x0$drugname)
		x1 = as.data.frame(matrix(0, nr=length(medID), nc=length(drugID)+1))
		x1[,1] = medID
		colnames(x1) = c(colnames(x0)[1], drugID)
		#p01 = date()
		for(i0 in 1:nrow(x1)){ #cat(i0,",")
			#x2 = unique(as.character(apply(x0[which(x0$regid_anon==x1[i0,1]),2:3], 1, function(x){gsub("-NA","",paste(x, collapse="-"))}))) # >4000 Others medicine diversity; lack of memory
			x2 = unique(x0[which(x0$regid_anon==x1[i0,1]),"drugname"])
			x1[i0,x2] = 1
		}#;cat("\n");p01;date()
## merge reformatted medication histories to remaining data
		x0 = merge(x1, cf425[[p[grep("NTM", p)]]], all=T)
		x = merge(cf425[[p[grep("review", p)]]], x0, all=T)
		rm(x0, x1, x2, p, mEd)
	}
	x$year = yR[i]
	colnames(x) = tolower(colnames(x))
	rEc = merge(rEc,x, all=T)
}#;cat("\n");d0;date()

##### patch "other species" record #####
zZ = gsub("[|]+","@",gsub("[.*]"," ",gsub("[\n]","@",gsub("[;,] ","@",gsub("_x000d_","",gsub(" mar, jul, sept, nov 2019, jan 2020","",gsub("yeast[/,: ]","yeast@",gsub("yeasts[/,: ]","yeast@",gsub("yeasts and ","yeast@",gsub("yeasts &","yeast@",
	c(rEc$cult_specify,rEc[,"s05culturespeciesresistotherspec"],rEc[,"s05culturespeciesfungalotherspec"],rEc[,"s05culturespeciesviralotherspeci"])
	))))))))))
zZ[which(zZ=="?")] = zZ[which(zZ==0)] = NA
zZ = trimws(gsub("^[,?]|,$","",trimws(gsub("\\s+", " ", tolower(zZ)), which="both")), which="both")
zZ = zZ[-grep("sputum",zZ)] # no extra info & machine-confuse symbols
for(i in c("heavy ", "a moderate ", "scanty ", "light ", "")){zZ = gsub(paste0(i,"growth of "),"",zZ)}
for(i in c(" [/]+ ","\n")){zZ = gsub(i,"@",zZ)}
zZ = trimws(unique(zZ),which="both")
zZ = gsub("\'","!!",zZ)
for(i in grep("@",zZ)){zZ = c(zZ,as.character(read.table(text=zZ[i], sep="@")[1,]))}
zZ = zZ[-grep("@",zZ)]
zZ = gsub("!!","\'",unique(zZ))
write.table(zZ[order(zZ)], "../data/otherSp.csv", sep=",", quote=F, row.names=F, col.names=F)

#zZ0 = as.data.frame(table(zZ))
#write.csv(zZ0[,c(2,1)], "../data/otherSp.csv", quote=F, row.names=F)
#write.table(zZ0[,c(2,1)], "../data/otherSp-ExSep.csv", quote=F, row.names=F, sep="!")
#rEc[,"s05culturespeciesresistotherspec"] = zZ

##### extract annual review patients & useful data attributes ##### 20220212
dRop = c(13,15,17:22,25,33,43,seq(56,96,2),97:106,145,146,269,275,281,283,284,291,292,308:310,seq(318,324,2),328:330,334:339,366,369,373,376,377,381:394,402:409,411:414,416:453,seq(456,480,2),481:500,528:591,seq(598,616,2),617:619,seq(628,644,2),645:654)
regID = c();for(i in yEar){regID = unique(c(regID, cf425[[i]]$regid_anon))};rm(i)
rEc = rEc[which(rEc$regid_anon %in% regID),-dRop]
# Only patients with annual review presence in data + question multichoices are "tick" or "not ticked" = unticked entries == negative response instead of no info (i.e. 0 instead of NA)

##### standardize data ##### 20220212, 20220315
for(i in 3:ncol(rEc)){rEc[,i] = tolower(rEc[,i])}
rEc[is.na(rEc)] = rEc[rEc=="no"] = rEc[rEc=="n"] = rEc[rEc=="nk"] = rEc[rEc=="unknown"] = 0
mUlti = c("cult_specify","s05culturespeciesresistotherspec","macrolidespec","enzymespecifyother","oralbroncho_bamed","oralbroncho_theophmed","inhbroncho_sabamed","inhbroncho_labamed","inhbroncho_saacmed","inhbroncho_laacmed","inhbroncho_baacmed","other_antibiotics_specify","s05culturespeciesfungalotherspec","s05culturespeciesviralotherspeci")
inDir = which(colnames(rEc)%in%mUlti)
dirRef = rEc[,-c(1,2,inDir)];indRef = rEc[,c(1,2,inDir)] # sep (in)direct extractable info 20220323
dirRef[dirRef!=0] = 1 # restrict to presence/absence data per patient
dirRef = cbind(indRef[,1:2],dirRef)
colnames(dirRef) = trimws(colnames(dirRef), which="both")

##### f: listing column details ##### 20220212
dTl = function(x){
        x0 = unique(rEc[,x])
        cat(colnames(rEc)[x], ",", length(x0),",",class(rEc[,x]), "\n")
        if(length(x0)<10){cat(x0, "\n");print(table(rEc[,x]))}}
#gWeird = function(x){return(grep(x, colnames(rEc0)))}
#for(i in 3:ncol(rEc)){if(length(unique(rEc[,i]))>10){dTl(i)}};rm(i)

##### extract column (medicine + standard species) names ##### 20220227
dEtails = which(colnames(rEc) %in% mUlti[c(3:12)])
zZ = colnames(rEc)[-c(1,2,245:252,374,dEtails)]
for(i in dEtails){zZ = c(zZ,unique(rEc[,i]))}
zZ = unique(zZ[-which(zZ=="0" | zZ=="nil" | zZ=="not taking")])
zZ1 = setdiff(grep(", ",zZ), grep(")",zZ))
zZ0 = read.table(text=zZ[zZ1], sep=",")
zZ = zZ[-zZ1]
for(i in 1:ncol(zZ0)){zZ = c(zZ, trimws(zZ0[,i], which="both"))}
zZ = zZ[-grep("^[[:digit:]]",zZ)]
zZ = unique(zZ)
#zZ = unique(c(colnames(rEc)[-c(1,2,245:250,374,375)],sub("^ ","",read.table(text=zZ[-which(zZ=="0" | zZ=="nil")],sep=",")[,1])))
write.table(c("input",zZ[order(zZ)]),"../data/genCols.csv", sep="\t",quote=F,row.names=F, col.names=F)

save(rEc,yEar,yR,mUlti,dirRef,indRef, file="../data/cf425FULL.rda")
