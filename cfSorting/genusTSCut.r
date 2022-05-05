#!/bin/env Rscript
# author: ph-u
# script: genusTSCut.r
# desc: cut genus time-series data in multiple csv
# in: Rscript genusTSCut.r
# out: gTS_[startYr][endYr]-gLV.csv
# arg: 0
# date: 20220501

##### env #####
a = read.csv("../data/genusTimeSeries_gLV.csv", header=T)

##### export start/end year data segregations #####
for(i0 in 1:(nrow(a)-2)){for(i1 in 3:nrow(a)){if(i1-i0>=2){
	a0 = a[which(a$year>=a$year[i0] & a$year<=a$year[i1]),]
	write.csv(a0,paste0("../data/gTS_",substr(as.character(a$year[i0]),3,4),substr(as.character(a$year[i1]),3,4),"_gLV.csv"), row.names=F, quote=F)
}}}
