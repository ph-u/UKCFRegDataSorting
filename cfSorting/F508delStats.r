#!/bin/env Rscript
# author: ph-u
# script: F508delStats.r
# desc: statistical tests for F508del group
# in: Rscript F508delStats.r
# out: ?
# arg: 0
# date: 20221207

## f: time comparison for same pair
tComp = function(d1,d2,tAg,mAx=500,nRep=7){
	eCo = unique(d1$c1_is)
	tAr = c("Eco","Rep","d1","d2")
	d0 = as.data.frame(matrix(0,nr=length(eCo)*nRep, nc=length(tAr)))
	colnames(d0) = tAr
	d0$Eco = eCo; d0$Rep = rep(1:nRep,each=length(eCo))

	d10 = d1; d20 = d2
	for(i in 1:length(tAg)){
		d10 = d10[grep(tAg[i],d10[,i]),]
		d20 = d20[grep(tAg[i],d20[,i]),]
	};rm(i)

	for(i in 1:nrow(d10)){d0$d1[which(d0$Eco==d10$c1_is[i] & d0$Rep==d10$replicate[i])] = d10$count[i]}
	for(i in 1:nrow(d20)){d0$d2[which(d0$Eco==d20$c1_is[i] & d0$Rep==d20$replicate[i])] = d20$count[i]}
	d0[,3:4] = d0[,3:4]/mAx

	d00 = data.frame(Eco=eCo,stat=NA,p.val=NA,adj.p=NA)
	for(i in 1:nrow(d00)){
		x = d0[which(d0$Eco==d00$Eco[i]),]
		x0 = wilcox.test(x$d1,x$d2,p.adjust.methods="none")
		d00[i,c(2,3)] = c(unname(x0$statistic), x0$p.value)
	};rm(i,x)
	d00$adj.p = p.adjust(d00$p.val,method="BH")

	cat("Any adj p-value < 0.05:",any(d00$adj.p<.05, na.rm=T),"\n")
	return(list(data=d0,result=d00))
}

##### 2012-14 vs 2018-20 #####
q1 = read.csv("../F508del_data/cftrM_0110_1214_gLV-eco.csv", stringsAsFactors=F)
q2 = read.csv("../F508del_data/cftrM_0110_1820_gLV-eco.csv", stringsAsFactors=F)
q3 = read.csv("../F508del_data/cftrM_0110_1618_gLV-eco.csv", stringsAsFactors=F)
q4 = read.csv("../F508del_data/cftrM_0110_0810_gLV-eco.csv", stringsAsFactors=F)

cat("14 vs 20\n")
PSa14v20 = tComp(q1,q2,c("Pseu","Staph"))
PSe14v20 = tComp(q1,q2,c("Pseu","Steno"))

cat("14 vs 18\n")
PSa14v18 = tComp(q1,q3,c("Pseu","Staph"))
PSe14v18 = tComp(q1,q3,c("Pseu","Steno"))
PSa14v18$data[grep("prey", PSa14v18$data$Eco),] # 2012-14 data all 0
PSe14v18$data[grep("mutu", PSe14v18$data$Eco),] # 2012-14 data all 0

cat("18 vs 10\n")
PSa18v10 = tComp(q3,q4,c("Pseu","Staph"))
PSe18v10 = tComp(q3,q4,c("Pseu","Steno"))

##### drug interaction influence #####
p1 = read.csv("../F508del_data/cftrM_1110_1517_gLV-eco.csv", stringsAsFactors=F)
p2 = read.csv("../F508del_data/cftrM_1111_1517_gLV-eco.csv", stringsAsFactors=F)

cat("Drug interaction influence: without vs with\n")
drugInt = tComp(p1,p2,c("Pseu","Staph"))
drugInt$data[which(drugInt$data$Eco=="mutualism"),] # 1110 group all 0
drugInt$data[grep("preda", drugInt$data$Eco),] # 1111 group not all 0

##### Same pair, ecological diff wihtin same period ##### 20221221
tEst = c("0911","1012","1719")
x = c("W","p","adj.p")
wxBH = as.data.frame(matrix(nr=length(tEst), nc=length(x)))
colnames(wxBH) = x;rm(x)
for(i0 in 1:length(tEst)){
	d1 = read.csv(paste0("../F508del_data/cftrM_0110_",tEst[i0],"_gLV-eco.csv"), stringsAsFactors=F)
	d1 = d1[which(d1$category1=="Pseudomonas" & d1$category2=="Stenotrophomonas"),]
	eLim = c();for(i in unique(d1$c1_is)){if(any(d1$count[which(d1$c1_is==i)]!=0)){eLim = c(eLim,i)}};rm(i)
	d1 = d1[which(d1$c1_is %in% eLim),]
	w0 = wilcox.test(d1$ratio_in_rep[which(substr(d1$c1_is,1,4)=="comm")],d1$ratio_in_rep[which(substr(d1$c1_is,1,4)!="comm")])
	wxBH[i0,-ncol(wxBH)] = c(w0$statistic, w0$p.value)
};rm(i0,w0)
wxBH[,ncol(wxBH)] = p.adjust(wxBH[,2],method="BH")

d1 = read.csv(paste0("../F508del_data/cftrM_1111_1315_gLV-eco.csv"), stringsAsFactors=F)
d1 = d1[which(d1$category1=="Mycobacteroides" & d1$category2=="Pseudomonas"),]
eLim = c();for(i in unique(d1$c1_is)){if(any(d1$count[which(d1$c1_is==i)]!=0)){eLim = c(eLim,i)}};rm(i)
d1 = d1[which(d1$c1_is %in% eLim),]
x0 = d1$ratio_in_rep[which(substr(d1$c1_is,1,4)=="comm")]
x1 = d1$ratio_in_rep[which(substr(d1$c1_is,1,4)!="comm")]
wilcox.test(c(x0,rep(0,7-length(x0))),c(x1,rep(0,7-length(x1))))
