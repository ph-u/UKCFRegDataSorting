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

##### medication names #####
medNum = which(colnames(indRef) %in% mUlti[c(3:12,15:length(mUlti))])
medNam = c();for(i in medNum){medNam = c(medNam,unique(indRef[,i]))}
medMtx = read.table(text = gsub("[']","!",medNam), sep="@", fill=T)
mL = c();for(i in 1:ncol(medMtx)){mL = c(mL,medMtx[,i])} # ln969323
mL = unique(trimws(mL[!is.na(mL)], which="both")) # ln3913
mL = mL[-grep("^[^a-zA-Z]+$",mL)] # ln3866
mL = mL[mL!=""] # ln3865

##### rough rm entries with quantity indications #####
gM = grep("*[0-9][ mg][ g]+ *", mL);gM0 = mL[-gM];gM1 = mL[gM]
gM1 = unique(trimws(read.table(text = sub("[0-9]","@",gM1), sep="@", fill=T)[,1], which="both"))
gM0 = unique(c(gM0, gM1[gM1!=""])) # all duplicates
gM0 = gM0[-grep("^[0-9]",gM0)] # duplicate content / deducible; ln3536
gM = grep(" m[cg]", gM0);mL0 = gM0[-gM];mL1 = gM0[gM] # ln3213,323
mL0 = unique(c(mL0, trimws(read.table(text = sub("[0-9]","@",mL1), sep="@", fill=T)[,1], which="both"))) # ln3255
gM1 = read.table(text = mL0, sep=",", fill=T)[,1] # ln3258

##### meaningless phrases replacement #####
gM = c(" 2units before breakfast and lunch sc", "[|]", "chronic", "half ", " od ", "rescue", "unknown", " and ", "[+]", " - ", "slow", "normal", "adult", "liquid", "alt days", "oral ", "as required")
for(i in gM){gM1 = gsub(i,"@",gM1)}
gM = grep("@",gM1);mL0 = gM1[-gM];mL1 = unique(gM1[gM]) # ln3133,122
mL1 = read.table(text=mL1, sep="@", fill=T)
for(i in 1:ncol(mL1)){mL0 = c(mL0,mL1[,i])}
gM1 = unique(trimws(tolower(mL0), which="both")) # ln3148

##### meaningless entries removal #####
gM = c("asdirected", "aviva strips", "bg star strips", "cassette", "complan", "course complete", "dario", "mart up to", "microscoop", "not known", "l capsules", "ml eye drops", "po 30-60ml od", "^for ", "^day", "freestyle", "test strip", "testing strip", "butter", "fast clix", "fastclix")
for(i in gM){gM1 = gM1[-grep(i,gM1)]} # ln3109
gM = c("d", "bd", "may", "week", "other", "thick", "easy", "testing strips", "discomfort", "from feb to april", "", "cold", "for enuresis", "accu chek performa", "feed", "well teen tablets", "during hayfever season", "wk", "wfi", "pollen allergy", "lp", "days", "alt months", "day", "l", "r", "for 3 months", "nil", "wexham prescribed", "release", "2.5ml bd", "1.5 sachets per day", "not taking", "no current medications", "nil", "fluid thickening agent", "blood glucose stripes", "bone health", "blinded", "as req", "1.5 sachets per day", "2.5ml bd", "4", "0", "od", "dose being titrated up", "minute", "null", "rarely needed", "retired medication", "smart")
for(i in gM){gM1 = gM1[which(gM1 != i)]} # ln3066

##### remove quantities #####
gM = intersect(grep("[0-9]", gM1), unique(c(grep("m[lc]", gM1), grep(" per ", gM1), grep("unit", gM1), grep("sachet", gM1), grep(" ug", gM1))))
mL0 = gM1[-gM];mL1 = unique(gM1[gM]) # ln2952,114
mL0 = unique(c(mL0,trimws(read.table(text=sub(" [0-9]","@",mL1), sep="@", fill=T)[,1], which="both"))) # 2969

##### export #####
write.table(c("input",gsub("!","'",mL0[order(mL0)])),"../data/medName.csv", sep="\t",quote=F,row.names=F, col.names=F)

#medNam = unique(medNam[-which(medNam=="0" | medNam=="nil" | medNam=="not taking" | medNam=="?")])
#for(i in c(" [[:digit:]]+[m ]", " [[:digit:]]+[.][[:digit:]]+[m ]", "slow", " 2units", " 40g", " open ", " as ")){medNam = gsub(i,"@0",medNam)}
#for(i in c("breakfast and lunch sc", " and ", "20 mg", "[+]", "oral ", " or ", "[|]", ", ", "asdirected", "aviva strips", "berry", "bg star strips", "cassette", "chronic", "complan", "continuous", "course complete", "dario", "half ", "mart up to", "microscoop", "not known", "od ", "other", "rescue", "unknown")){medNam = gsub(i,"@",medNam)}
#medN = medNam[-grep("@",medNam)];medN0 = medNam[grep("@",medNam)]

## separate concatenated medicine names
#i = read.table(text = gsub("[']","!",medN0), sep="@", fill=T)
#i[i=="0"] = i[i=="null"] = i[i==""]= NA
#mD = c();for(i0 in 1:ncol(i)){mD = c(mD,i[which(!is.na(i[,i0])),i0])}
#medN = unique(trimws(c(medN,gsub("!","'",mD)), which="both"))

## remove non-medicine/nutrient entries (no macronutrient)
#i0 = c();for(i in c("^[[:digit:]]", "glucos", "thick", "being", "epipen", "discomfort", "drink", "gaviscon", "ointment", "energy", "oil", "w[efk]", "may", "from", "allergy", "cold", "zer", "needles", "nutrisan", "minute", "smart", "po 3", "nail", "ensure", "two cal", "twocal", "feed", "electrolytes", "testing strips", "nutri[ispl]", "current", "cream", "for ", "emolient", "just ", "day", "as ", "^trial", "freestyle", "during ", "^liquid ", "^powder ", "^rarely ", "reusable", "^accu")){
#	i0 = c(i0,grep(i,medN))};medN = medN[-unique(i0)]
#for(i in c("r", "d", "l", "bd", "tip", "patch", "nutrini", "null", "for", "", "-", "?", ".9", "easy", "bone health", "adult", "alt months", "mobile cassette", "scoops", "retired medication")){medN = medN[-which(medN==i)]}

##### export #####
#write.table(c("input",medN[order(medN)]),"../data/medName.csv", sep="\t",quote=F,row.names=F, col.names=F)
