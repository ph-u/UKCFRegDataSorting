#!/bin/env Rscript
# author: ph-u
# script: otherSp.r
# desc: patch "other species" record
# in: Rscript otherSp.r
# out: ../data/otherSp.csv
# arg: 0
# date: 20220210, 20220414

load("../data/cf425FULL.rda")

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

