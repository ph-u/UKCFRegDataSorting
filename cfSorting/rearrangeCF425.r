#!/bin/env Rscript
# author: ph-u
# script: rearrangeCF425.r
# desc: CF-Registry data xlsx reorganize
# in: Rscript rearrangeCF425.r
# out: ../data/cf425.rda
# arg: 0
# date: 20220131

##### Convert xlsx into RData ##### 26 mins
library(xlsx2dfs)
cf425 = xlsx2dfs("../p_cf425/Data for 425 - final.xlsx", rowNames=F, detectDates = T) ## long runtime (31 sheets, 6000s row * 100s col mostly each)

##### unique patients in each data set #####
cOl = c("dfNam","uniqPatient","totalEntry","nCol")
uniqNam = as.data.frame(matrix(NA, nr=length(cf425), nc=length(cOl)))
colnames(uniqNam) = cOl
for(i in 1:length(cf425)){
        uniqNam[i,] = c(names(cf425)[i], length(unique(cf425[[i]][,"regid_anon"])), dim(cf425[[i]]))
}

##### unique patients & data attributes in time-series #####
yEar = c(which(nchar(names(cf425))==4), grep("annual", names(cf425)))

save(cf425, uniqNam, yEar, file="../data/cf425.rda")

