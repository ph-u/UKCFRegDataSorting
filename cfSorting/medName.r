#!/bin/env Rscript
# author: ph-u
# script: medName.r
# desc: extract medicine name
# in: Rscript medName.r
# out: medName.csv
# arg: 0
# date: 20220422

##### env #####
load("../data/cf425FULL.rda")

##### medicine names (anti-microbes only) #####
medNum = which(colnames(indRef) %in% mUlti[c(3:12,15:length(mUlti))])
medNam = c();for(i in medNum){medNam = c(medNam,unique(indRef[,i]))}
medNam = unique(medNam[-which(medNam=="0" | medNam=="nil" | medNam=="not taking" | medNam=="?")])
for(i in c(" [[:digit:]]+[m ]", " [[:digit:]]+[.][[:digit:]]+[m ]", "slow", " 2units", " 40g", " open ", " as ")){medNam = gsub(i,"@0",medNam)}
for(i in c(" and ", "20 mg", "[+]", "oral ", " or ", "[|]", ", ")){medNam = gsub(i,"@",medNam)}
medN = medNam[-grep("@",medNam)];medN0 = medNam[grep("@",medNam)]

## separate concatenated medicine names
i = read.table(text = gsub("[']","!",medN0), sep="@", fill=T)
i[i=="0"] = i[i=="null"] = i[i==""]= NA
mD = c();for(i0 in 1:ncol(i)){mD = c(mD,i[which(!is.na(i[,i0])),i0])}
medN = unique(trimws(c(medN,gsub("!","'",mD)), which="both"))

## remove non-medicine entries
i0 = c();for(i in c("^[[:digit:]]", "glucos", "thick", "being", "epipen", "discomfort", "drink", "gaviscon", "ointment", "energy", "oil", "w[efk]", "may", "from", "allergy", "cold", "zer", "needles", "nutrisan", "minute", "smart", "po 3", "nail", "ensure", "two cal", "twocal", "feed", "electrolytes", "testing strips", "nutri[ispl]", "current", "cream", "for ", "emolient", "just ", "day", "as ", "^trial", "freestyle", "during ", "^liquid ", "^powder ", "^rarely ", "reusable", "^accu")){
	i0 = c(i0,grep(i,medN))};medN = medN[-unique(i0)]
for(i in c("r", "d", "l", "bd", "tip", "patch", "nutrini", "null", "for", "", "-", "?", ".9", "easy", "bone health", "adult", "alt months", "mobile cassette", "scoops", "retired medication")){medN = medN[-which(medN==i)]}

##### export #####
write.table(c("input",medN[order(medN)]),"../data/medName.csv", sep="\t",quote=F,row.names=F, col.names=F)
