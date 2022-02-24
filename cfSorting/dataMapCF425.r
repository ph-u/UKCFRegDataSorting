#!/bin/env Rscript
# author: ph-u 
# script: dataMapCF425.r
# desc: cf425 data remapping
# in: Rscript dataMapCF425.r
# out: none
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
zZ = gsub("[|]+",",",gsub("[.*]"," ",gsub("[:;\n]",",",rEc[,"s05culturespeciesresistotherspec"])))
zZ[which(zZ=="?")] = NA
zZ = trimws(gsub("^[,?]|,$","",trimws(gsub("\\s+", " ", tolower(zZ)), which="both")), which="both")
zZ0 = as.data.frame(table(zZ))
write.csv(zZ0[,c(2,1)], "../data/otherSp.csv", quote=F, row.names=F)
#write.table(zZ0[,c(2,1)], "../data/otherSp-ExSep.csv", quote=F, row.names=F, sep="!")
rEc[,"s05culturespeciesresistotherspec"] = zZ
save(cf425,rEc,yEar,yR, file="../data/cf425FULL.rda")
