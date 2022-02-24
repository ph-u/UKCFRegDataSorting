#!/bin/env Rscript
# author: ph-u
# script: dataStd.r
# desc: cf425 data standardization
# in: Rscript dataStd.r
# out: none
# arg: 0
# date: 20220212

##### env #####
library(stringr)
load("../data/cf425FULL.rda")

##### exclude non-useful columns & patients without annual reviews #####
dRop = c(13,15,17:22,25,33,43,seq(56,96,2),97:106,145,146,269,275,281,283,284,291,292,308:310,seq(318,324,2),328:330,334:339,366,369,373,376,377,381:394,402:409,411:414,416:453,seq(456,480,2),481:500,528:591,seq(598,616,2),617:619,seq(628,644,2),645:654)
regID = c();for(i in yEar){regID = unique(c(regID, cf425[[i]]$regid_anon))};rm(i)
rEc = rEc[which(rEc$regid_anon %in% regID),-dRop]
# Only patients with annual review presence in data + question multichoices are "tick" or "not ticked" = unticked entries == negative response instead of no info (i.e. 0 instead of NA)

##### standardize data recording method #####
rEc[is.na(rEc)] = 0
rEc[rEc=="SELECT"] = rEc[rEc=="seleted"] = rEc[rEc=="selected"] = rEc[rEc==2] = 1 # restrict to presence/absence data per patient
# ntm species treatment = species detected = presence in microbiome
for(i in 3:ncol(rEc)){rEc[,i] = tolower(rEc[,i])};rm(i)
nkRec = c(grep("nk", rEc),grep("unknown", rEc))
for(i in 1:length(nkRec)){if(length(unique(rEc[,nkRec[i]]))==4){
        x = rEc[,nkRec[i]]
        x[x=="nk"] = x[x=="unknown"] = x[x==0] = .5
        x[x=="n"] = x[x=="no"] = 0
        x[x=="y"] = x[x=="yes"] = 1
        rEc[,nkRec[i]] = x
        # unknown = no info; both have no clear indication
}}; rm(i,x)
x = rEc[,"colistimethate"];x[x==1]="yes";rEc[,"colistimethate"] = x;rm(x) # "SELECT" = drug presence in treatment = "yes"

##### f: listing column details #####
dTl = function(x){
        x0 = unique(rEc[,x])
        cat(colnames(rEc)[x], ",", length(x0),",",class(rEc[,x]), "\n")
        if(length(x0)<10){cat(x0, "\n");print(table(rEc[,x]))}}
#gWeird = function(x){return(grep(x, colnames(rEc0)))}
#for(i in 3:ncol(rEc)){if(length(unique(rEc[,i]))>10){dTl(i)}};rm(i)

##### f: add new named column #####
nwCol = function(nAm, tX="", eX="", df=rEc, rEf="s05culturespeciesresistotherspec"){
	if(tX[1]==""){tX=nAm}
	tX = unique(c(nAm, tX))
	df[,nAm] = 0; rC = c()
	for(i in tX){rC = c(rC,grep(i, df[,rEf]))}
	rC = unique(rC)
	if(eX[1]!=""){rC0=c() # exclude a subset of selected text
		for(i in eX){rC0 = c(rC0,grep(i, df[rC,rEf]))}
		rC0 = rC[unique(rC0)] # row potentially be excluded
		for(i in rC0){
			r1=c(); for(i0 in tX){r1 = c(r1, str_count(df[i,rEf],i0))}
			r2=c(); for(i0 in eX){r2 = c(r2, str_count(df[i,rEf],i0))}
			if(sum(r1)<=sum(r2)){rC = rC[-which(rC==i)]}
		}}; df[unique(rC),nAm] = 1
	return(df)}

#rEc0 = rEc
##### break species collections into separate columns (highest consistent resolution Genus) #####
## category -> col
x = c("pseudofischeri", "citrobacter", "microccus", "klebsiella", "haemophilus", "elizabethkingia", "moraxella", "delftia", "proteus", "brevundimonas", "sedosporium", "rahnella", "burkholderia", "pneumococcus", "lelliottia", "coliform", "paecilomyces", "rasamsonia", "cupriavidus", "capnocytophaga", "stenotrophomonas", "variovorax", "rothia", "pasturella", "neisseria", "bordetella", "lomentospora", "marcescenes", "pnemocystis", "trichosporon", "nocardia", "lichtheimia", "aeromonas")

x = c("yersinia", "acinomycens", "pneumococcus", "pneumocystis", "cryptococcus")
for(i in x){rEc = nwCol(i)};rm(i,x)

## text equivalents in a single category
#uQ = unique(rEc[,"s05culturespeciesresistotherspec"]); uQ = uQ[order(uQ)]
#uQ[grep("ye", uQ)]
#uQ = gsub("ye","(Z)",uQ) ##!!! enterobacter pittii
#x9 = c();for(i in x9){print(uQ[grep(i, uQ)])};rm(i)
#x9 = c();for(i in x9){uQ = gsub(i,"(Z)",uQ)};rm(i)
rEc = nwCol("", c(""))

xOTH = c("ruh", "xy[lc]o", "insu", "r do[lc]", "r den", " dele", "wof", " piec", "inso", "marp", "muci")
rEc = nwCol("AchromobacterMucicolens", xOTH[11])
rEc = nwCol("AchromobacterMarplatensis", xOTH[10])
rEc = nwCol("AchromobacterInsolitus", xOTH[9])
rEc = nwCol("AchromobacterPiechaudii", xOTH[8])
rEc = nwCol("AchromobacterLwoffii", xOTH[7], c("acin"))
rEc = nwCol("AchromobacterDeleyi", xOTH[6])
rEc = nwCol("AchromobacterDenitrificans", xOTH[5])
rEc = nwCol("AchromobacterDolens", xOTH[4])
rEc = nwCol("AchromobacteriInsuavis", xOTH[3])
rEc = nwCol("AchromobacteriXylosoxidans", xOTH[2])
rEc = nwCol("AchromobacteriRuhlandii", xOTH[1])
rEc = nwCol("achromobacter", c("achr", "achromobacter sp"), c("baum", xOTH))
rEc = nwCol("chimera", "", c("mycobacterium")) # virus

xOTH = c("muil","muc[eio]",)
rEc = nwCol("RothiaMucilaginosa", xOTH, c("rhodo", "coryne", "roseomonas", "r mucilaginosa", "mucoid", "mycobac", "mucor[- ]", "glutinis[.] mucilag", "achro"))
rEc = nwCol("rothia", c("roth"), xOTH)

cOTH = c("gani", "morgani", "morgagnii", "morgan ll")
rEc = nwCol("Morganella morganii", cOTH)
rEc = nwCol("morganella", c("morg"), cOTH)

xMCH = c("cat","osl")
rEc = nwCol("MoraxellaOsloensis", xMCH[2])
rEc = nwCol("MoraxellaCatarrhalis", xMCH[1], c("[io]cat"))
rEc = nwCol("moraxella", c("mora"), xMCH)
rEc = nwCol("MannheimiaHaemolytica", c("mannhaeimia haemolytica"))

xHPH = c("lus parah", "h[.]parah"); xHHL = c("lus h", "haem h", "haemoph[.] h"); xHPI = c("haem.parai", "lus parai", "haem paaraaflu", "haemph[.] parainfluenzae", "lusparai", "lus paran", "lus paraainf", "lus para in", "lus para- in", "lua parai", "lis parai", "haemoph[.] parai", "haemop[.] parai", "haemoph pari", "yticus parai", "haem. parai", "haem para", "h. parai"); xHIF = c("lus inf", "lis inf", "haemi inf")
rEc = nwCol("HaemophilusParahaemolyticus", xHPH)
rEc = nwCol("HaemophilusHaemolyticus", xHHL)
rEc = nwCol("HaemophilusParainfluenzae", xHPI)
rEc = nwCol("HaemophilusInfluenzae", xHIF)
rEc = nwCol("haemophilus", c("haemop", "hamop"), c("beta ", xHPH, xHHL, xHPI, xHIF))

xOTH = c("braa", "kos", "freu", "frendii")
rEc = nwCol("CitrobacterBraakii", xOTH[1])
rEc = nwCol("CitrobacterKoseri", xOTH[2])
rEc = nwCol("CitrobacterFreundii", xOTH[3:4])
rEc = nwCol("citrobacter", c("citr"), xOTH)
rEc = nwCol("roseobacter", c("rcb")) # Roseobacter clade bacteria

xSLF = c("liq", "serratia ligue"); xSMC = c("marc", "marsc", "marse", "mares", "marac", "merc", "maras", "serratia m", "mace", "masc", "serratia m"), c("12th")
rEc = nwCol("SerratiaLiquefaciens", xSLF)
rEc = nwCol("SerratiaMarcescens", xSMC, c("12th"))
rEc = nwCol("serratia", c("seratia", "serr", "sertia"), c(xSLF, xSMC))

xEMS = c("meningos", "meningo s", "maningos"); xEMC = c("miric", "minic", "minc", "miraco")
rEc = nwCol("ElizabethkingiaMeningoseptica", xEMS)
rEc = nwCol("ElizabethkingiaMiricola", xEMC)
rEc = nwCol("elizabethkingia", c("eliz"), c(xEMS, xEMC))

xOTH = c("aerog", "asbur", "kobei")xECC = c("cloa", "clocae", "colacae", "cloccae", "clocacae", "clobcae", "chloache")
rEc = nwCol("EnterobacterAerogenes", xOTH[1])
rEc = nwCol("EnterobacterAsburiae", xOTH[2])
rEc = nwCol("EnterobacterKobei", xOTH[3])
rEc = nwCol("EnterobacterCloacae", xECC)
rEc = nwCol("enterobacter", c("enterobact", "e spp", "esbl"), c(xOTH, xECC))
rEc = nwCol("KlebsiellaOxytoca", c("oxy"), c("carboxy", "oxyspurum"))

xOTH = c("metapneumo", "pneumococcus", "pneumocystis", "pneumotropica", "sppneumonia", "mycoplasma")
rEc = nwCol("KlebsiellaPneumoniae", c("pneu", "kleb pn"), c("strep", xOTH))
rEc = nwCol("klebsiella", c("kleb"), c("pneu", "oxy", "ticked above", "kleb pn", "previously called"))
rEc = nwCol("StreptococcusPyogenes", c("pyog"))
rEc = nwCol("StreptococcusPneumoniae", c("pneu", "strep pn", "piemonaie", "pnemonae", "pnuemonaie"), c("klebsiella", xOTH))
rEc = nwCol("streptococcus", c("step", "strep", "bhs"), c("pneu", "pyog", "strep pn", "piemonaie", "pnemonae", "pnuemonaie"))

xOTH = c("juni", "johnson", "wof", "pitii", "pitti", "ursing", "baum", "acinetobacter b")
rEc = nwCol("AcinetobacterJunii", xOTH[1])
rEc = nwCol("AcinetobacterJohnsonii", xOTH[2])
rEc = nwCol("AcinetobacterLwoffii", xOTH[3], "achr")
rEc = nwCol("AcinetobacterPittii", xOTH[4:5])
rEc = nwCol("AcinetobacterUrsingii", xOTH[6])
rEc = nwCol("AcinetobacterBaumanii", xOTH[7:8])
rEc = nwCol("acinetobacter", c("acin"), xOTH)
rEc = nwCol("mrsa", c("meticillin resistant s. aureus"))
rEc = nwCol("StaphylococcusAureus", c("aureus", "mssa", "tococcus a"), c("meticillin resistant", "agal", "angi"))
rEc = nwCol("staphylococcus", c("stap"), c("aureus", "indicate"))


rEc = nwCol("exophiala", c("exo"))

rEc = nwCol("CandidaAlbicans", c("alb"), c("ochrobactrum"))
rEc = nwCol("candida", c("cand"), c("albicans", "aspergillus candidus"))

rEc = nwCol("escherichia", c("escheric", "e.hermannii"), c("pseudallescheria", "coli"))
rEc = nwCol("ecoli", c("e coli,", "-coli","escherichia coli", "escherechia coli", "e[.]coli", "e[.] coli", "escherichiac coil"), c("positive coliform")) # avoid overkill with coliform & entries with both

rEc = nwCol("coliform", c("coli-f", "colifor", "colyform"))
rEc = nwCol("bacillus", c("baci", " rod", "zieh"))

rEc = nwCol("unidentified", c("carbapenem", "commensal", "flora", "org2", "mixed anerobes", "s mero", "urtf", "lac[.] ferm", "mucor[- ]", "r mucilaginosa"), c("(unidentified)", "unidentified fungus", "rhizomucor"))
# urtf = upper respiratory tract flora
rEc = nwCol("penicillium", c("pen"), c("penumoniae", "carbapenem", "-pencillin")) # worst: peniallium, penn spp, pen spp
rEc = nwCol("fungi", c("fung", "mould", "fugal"), c("not seen"))
rEc = nwCol("yeast", c("ye"), "yersinia") # common name, worst = yes
#########################################################################################
rEc = nwCol("pandoraea", c("pandoreae"))
rEc = nwCol("saccharomyces", c("saccharamyces"))
rEc = nwCol("serratia", c("s marc", "seratia", "serrratia"))
rEc = nwCol("exophiala", c("exophialia", "exo. dermatitidis", "exo.dermatitidis", "exo dermatitidis", "exophilia", "exophalia", "exo. dermatitidas"))
rEc = nwCol("C.Difficile", c("c difficile", "c. diff"))
rEc = nwCol("enterococci", c("entero facalis", "enterococcus", "enterococcous"))
rEc = nwCol("enterobacter", c("enterbacter"))
rEc = nwCol("mycobacterium", c("chimaera", "mac", "mycobacterius", "m. abcessus", "m. abcesses", "m.abcessus", "chelonae"))
rEc = nwCol("staphylococcus", c("staph"))
rEc = nwCol("streptococcus", c("strep", "bhs", "streptoccocus"))
rEc = nwCol("inquilinus", c("inquinolus", "inqulinis", "inquilininus"))
rEc = nwCol("CandidaAlbicans", c("candida sp/albicans", "candida albicans", "Candida albicans"))
rEc = nwCol("chryseobacterium", c("ch.indologenes", "chrtseobacterium indologenes"))
rEc = nwCol("penicillium", c("penicillin", "pencillium", "penici;;ium", "pencillin"))
rEc = nwCol("virus", c("parainfluenzae", "influenza ", "flu a", "h1n1", "chimera", "equivocal cmv"))
rEc = nwCol("EColi", c("e coli", "escherichia coli", "e-coli"))
rEc = nwCol("MSSA", c("mssa"))
rEc = nwCol("bacillus", c("bacilli", "ziehl-neilsen"))
rEc = nwCol("scedosporium", c("scedosporidium", "scedosporidium"))
rEc = nwCol("raoultella", c("raoutella"))
rEc = nwCol("geotrichum", c("geotrichum candidum"))
rEc = nwCol("morganella", c("morg. morganii", "morgarella"))
rEc = nwCol("Roseobacter", c("rcb"))
rEc = nwCol("rhodotorula", c("rhodoturala"))
rEc = nwCol("rhizopus", c("rhizpus"))
rEc = nwCol("pandorea", c("pandora"))
rEc = nwCol("COVID", c("covid-19"))
rEc = nwCol("ochrobactrum", c("ochrobactrium"))
rEc = nwCol("aspergillus", c("asp.fumigatus", "asp flavus"))
rEc = nwCol("achromobacter", c("achromobartar"))
rEc = nwCol("hafnia", c("havnia"))
rEc = nwCol("otherPseudomonas", c("pseudomonas", "psedomonas", "ps fluorescens"))
rEc = nwCol("otherCandida", c("candida", "Candida", "candidida"), c("albicans", "(not")) # every C but CA
rEc = nwCol("otherFungi", c("fungus", "mould", "fungal", "fungi"), "not seen")
rEc = nwCol("unidentified", c("e spp", "pen spp", "commensals", "sengnilparus rugosis", "segnilioparis rugalus", "segniliporus rugosus", "urtf", "flora", "contaminants", "gram positive cocci", "esbl"), "pandora")
rEc = nwCol("acinetobacter", c("ac.baumanii"))

head(unique(rEc[,"s05culturespeciesresistotherspec"]))
rEc[,"s05culturespeciesresistotherspec"] = NULL
