#!/bin/env python3
# author: ph-u
# script: medNameExtraction.py
# desc: collect standardized info
# in: python3 medNameExtraction.py
# out: drug-wiki.csv
# arg: 0
# date: 20220727

##### env #####
nAm = "drug";
import bs4, requests, re, os, copy, sys, timeit, string;
from collections import Counter;
os.chdir("../data"); sEp=",";
u = ["www.drugs.com/","www.medicines.org.uk/emc/product/","go.drugbank.com/drugs/","www.ndrugs.com/?s=","www.sdrugs.com/?c=drug&s=","pillintrip.com/medicine/","medlineplus.gov/"]

##### func #####
def limText(sTr): ## 20220729
	"""limit text length to max 6 words"""
	return sTr # no text limitation
	s1 = sTr.split()
	if len(s1)>6: s0 = [s1[0],s1[1],s1[2],s1[len(s1)-3],s1[len(s1)-2],s1[len(s1)-1]];
	else: s0 = s1;
	return ' '.join(s0)

def drugMTM(tExt): ## 20220728
	"""drugs.com mtm/other page info extraction"""
	w0 = ["" for w in range(4)]
	w1 = 1
	for w in range(len(tExt)):
		if isinstance(re.search("Drug class",tExt[w]),re.Match): w0[0] = ';'.join([w0[0],limText(tExt[w].split(":")[1].strip().lower())])
		if isinstance(re.search("Dosage form",tExt[w]),re.Match): w0[3] = ';'.join([w0[3],limText(tExt[w].split(":")[1].strip().lower())])
		for q in [" given by "]: # ingestion form
			if isinstance(re.search(q,tExt[w]),re.Match): w0[1] = ';'.join([w0[1],limText(tExt[w].split(q)[1].strip().lower())])
## drug function skip warnings/side effects
		if isinstance(re.search("How should",tExt[w]),re.Match) or isinstance(re.search("Dosing",tExt[w]),re.Match) or isinstance(re.search("More about ",tExt[w]),re.Match): w1 = 1
		if isinstance(re.search("Warnings",tExt[w]),re.Match) or isinstance(re.search(" side effects",tExt[w]),re.Match): w1 = 0
		if (w1==1) and isinstance(re.search("ask ",tExt[w].lower()),type(None)):
			for q in [" produce "," can "," helps "," help "]: # drug function
				if isinstance(re.search(q,tExt[w]),re.Match): w0[1] = ';'.join([w0[1],limText(tExt[w].split(q)[1].strip().lower())])
			for q in [" treat "," lead to ", " used in "," result from "," caused by "]: # disease
				if isinstance(re.search(q,tExt[w]),re.Match): w0[2] = ';'.join([w0[2],limText(tExt[w].split(q)[1].strip().lower())])
		if isinstance(re.search("Therapeutic Area:",tExt[w]),re.Match): w0[3] = ';'.join([w0[3],limText(tExt[w].split(":")[1].strip().lower())]) # diseasee
	return w0

def drugPRO(tExt): ## 20200728
	"""drugs.com pro page info extraction"""
	w0 = ["" for w in range(4)]
	for w in range(len(tExt)):
		if isinstance(re.search("Drug class",tExt[w]),re.Match): w0[0] = ";".join([w0[0],limText(tExt[w].split(":")[1].strip().lower())])
		if isinstance(re.search("Dosage form",tExt[w]),re.Match): w0[3] = ";".join([w0[3],limText(tExt[w].split(":")[1].strip().lower())])
		if isinstance(re.search("Indications ",tExt[w]),re.Match) and isinstance(re.search(" for ",tExt[w]),re.Match): w0[1] = w0[2] = ";".join([w0[1],limText(tExt[w+1].strip().lower())])
	return w0

def drugIGD(tExt): ## 20220730
	"""drugs.com ingredient page info extraction"""
	w0 = ["" for w in range(4)]
	tExt = "@".join(tExt)
	for i in list(string.ascii_uppercase):
		tExt = re.sub(i,"@"+i,tExt)
	tExt = [x for x in tExt.split("@") if len(x)>=3]
	w1 = 0
	for w in range(len(tExt)):
		if isinstance(re.search("rug class",tExt[w]),re.Match): w0[0] = ";".join([w0[0],limText(tExt[w].split(":")[1].strip().lower())]) # drug class
## disease bulk text
		if isinstance(re.search("Brand name",tExt[w]),re.Match) or isinstance(re.search("See also",tExt[w]),re.Match): w1 = 0
		if w1==1: w0[2] = ";".join([w0[2],limText(tExt[w].strip().lower())])
		if isinstance(re.search("treatment of",tExt[w]),re.Match): w1 = 1
	return w0

def drugINT(tExt): ## 20220731
	"""drugs.com international page info extraction"""
	w0 = ["" for w in range(4)]
	for w in range(len(tExt)):
		if isinstance(re.search("Therapeutic Category",tExt[w]),re.Match): w0[0] = ";".join([w0[0],limText(re.sub("Therapeutic Category","@",re.sub("Chemical Names","@",re.sub(":",";",tExt[w]))).split("@")[1].strip().lower())])
	return w0

def drugCONS(tExt): ## 20220731
	"""drugs.com cons page info extraction"""
	w0 = ["" for w in range(4)]
	w1 = 0; w2 = 0; w3 = 0
	for w in range(len(tExt)):
## drug class bulk
		if isinstance(re.search("Description",tExt[w]),re.Match): w1 = 0
		if w1 == 1: w0[0] = w0[0] = ";".join([w0[0],limText(re.sub(",",";",tExt[w])).strip().lower()])
		if isinstance(re.search("Category",tExt[w]),re.Match): w1 = 1
## disease
		for q in [" for "," lead to "," condition called "]:
			if isinstance(re.search(q,tExt[w]),re.Match): w0[2] = ";".join([w0[2],limText(tExt[w].split(q)[1]).strip().lower()])
		if isinstance(re.search(",",tExt[w]),re.Match): w2 = 0
		if w2 == 1: w0[2] = ";".join([w0[2],limText(tExt[w]).strip().lower()])
		if isinstance(re.search("Some conditions may",tExt[w]),re.Match) and isinstance(re.search("These include",tExt[w]),re.Match): w2 = 1
## ingestion form
		if isinstance(re.search("Importance of",tExt[w]),re.Match): w3 = 0
		if w3 == 1 and isinstance(re.search("\)",tExt[w].lower()),type(None)): w0[3] = ";".join([w0[3],limText(tExt[w]).strip().lower()])
		if isinstance(re.search("dosage forms:",tExt[w]),re.Match): w3 = 1
	return w0

def drugCOND(tExt): ## 20220731
	"""drugs.com disease condition page info extraction"""
	tExt = [x for x in list(filter(None,re.sub("\.","\n",bs4.BeautifulSoup(requests.get(uRl+"?page_all=1").text, "html.parser").get_text()).split("\n"))) if len(x)>=7]; # re-capture full webpage listing all documented drugs
	w0 = ["" for w in range(4)]
	for w in range(len(tExt)):
		if isinstance(re.search("Drug class",tExt[w]),re.Match): w0[0] = ";".join([w0[0],limText(re.sub("\t","",re.sub(", ",";",tExt[w+1])).strip().lower())])
		for q in [" treatments for "," treatment for "]: # disease
			if isinstance(re.search(q,tExt[w]),re.Match): w0[2] = ';'.join([w0[2],limText(tExt[w].split(q)[1].strip().lower())])
		if isinstance(re.search("Symptoms and treatments",tExt[w]),re.Match): w0[2] = ";".join([w0[2],limText(tExt[w+1].strip().lower())])
		if isinstance(re.search("Other names:",tExt[w]),re.Match): w0[2] = ";".join([w0[2],limText(tExt[w].split(":")[1].strip().lower())])
	return w0

def drugNPC(tExt): ## 20220731
	"""drugs.com npc page info extraction"""
	w0 = ["" for w in range(4)]
	w1 = 0; w2 = 0
	for w in range(len(tExt)):
		if isinstance(re.search("Drug class:",tExt[w]),re.Match): w0[0] = ";".join([w0[0],limText(tExt[w]).split(":")[1].strip().lower()])
## disease collection
		if isinstance(re.search(" treat ",tExt[w]),re.Match): w0[2] = ";".join([w0[2],limText(tExt[w]).split(" treat ")[1].strip().lower()])
		if isinstance(re.search("Further info",tExt[w]),re.Match): w1 = 0
		if w1 == 1: w0[2] = ";".join([w0[2],limText(tExt[w]).strip().lower()])
		if isinstance(re.search("Related treatment guides",tExt[w]),re.Match): w1 = 1
## info bulk collection
		if isinstance(re.search("What are the precautions",tExt[w]),re.Match): w2 = 0
		if w2 == 1:
			for q in [" help "]: # drug function collection
				if isinstance(re.search(q,tExt[w]),re.Match): w0[1] = ';'.join([w0[1],limText(tExt[w].split(q)[1].strip().lower())])
			for q in [" take "," put"]: # ingestion form
				if isinstance(re.search(q,tExt[w]),re.Match): w0[3] = ';'.join([w0[3],limText(tExt[w].split(q)[1].strip().lower())])
		if isinstance(re.search("What is this product used for",tExt[w]),re.Match): w2 = 1
	return w0

def mediORG(tExt): ## 20220729
	"""medicines.org.uk info extraction"""
	w0 = ["" for w in range(4)]
	w1 = 0
	for w in range(len(tExt)):
		if isinstance(re.search("Pharmacotherapeutic group",tExt[w]),re.Match): w0[0] = ";".join([w0[0],limText(tExt[w].split(":")[1].strip().lower())]) # drug class
		for q in [" substance"," leads to "," causes "]: # drug function
			if isinstance(re.search(q,tExt[w]),re.Match): w0[1] = ";".join([w0[1],tExt[w].split(q)[1].strip().lower()])
		for q in [" treatment of ","Drugs for "]: # disease
			if isinstance(re.search(q,tExt[w]),re.Match): w0[2] = ";".join([w0[2],limText(tExt[w].split(q)[1].strip().lower())])
## ingestion form bulk data collection
		if isinstance(re.search("Clinical particulars",tExt[w]),re.Match): w1 = 0
		if w1==1: w0[3] = ";".join([w0[3],limText(tExt[w].strip().lower())])
		elif w1==2: w1 = 1
		if isinstance(re.search("quantitative composition",tExt[w]),re.Match) and isinstance(re.search(" Name of ",tExt[w]),type(None)): w1 = 1
		if isinstance(re.search("Pharmaceutical form",tExt[w]),re.Match): w1 = 2
	return w0

def drugBank(tExt): ## 20220729
	"""drugbank info extraction"""
	w0 = ["" for w in range(4)]
	w1 = 0
	for w in range(len(tExt)):
		for q in [" is a "," as a "]: # drug function
			if isinstance(re.search(q,tExt[w]),re.Match): w0[0] = ";".join([w0[0],limText(tExt[w].split(q)[1].strip().lower())]) # drug class
		if isinstance(re.search(" treat",tExt[w]),re.Match): w0[2] = ";".join([w0[2],limText(tExt[w].split(" treat")[1].strip().lower())]) # disease
		if isinstance(re.search(" via ",tExt[w]),re.Match): w0[3] = ";".join([w0[3],limText(tExt[w].split(" via ")[1].strip().lower())]) # ingestion form
## disease bulk
		if isinstance(re.search("Contraindications",tExt[w]),re.Match): w1 = 0
		if w1==1: w0[2] = ";".join([w0[2],tExt[w].strip().lower()])
		if isinstance(re.search("Associated Conditions",tExt[w]),re.Match): w1 = 1
		for q in ["enhance","assist"]: # drug function
			if isinstance(re.search(q,tExt[w]),re.Match): w0[1] = ";".join([w0[1],limText(tExt[w].split(q)[1].strip().lower())])
		if isinstance(re.search("effect",tExt[w]),re.Match): w0[1] = ";".join([w0[1],limText(tExt[w].split("effect")[0].strip().lower())]) # drug function
	return w0

def nDrugs(tExt): ## 20220730
	"""ndrugs info extraction"""
	w0 = ["" for w in range(4)]
	for w in range(len(tExt)):
		for q in [" is a "]: # drug class
			if isinstance(re.search(q,tExt[w]),re.Match): w0[0] = ";".join([w0[0],limText(tExt[w].split(q)[1].strip().lower())])
		for q in [" effect "]: # drug class
			if isinstance(re.search(q,tExt[w]),re.Match): w0[0] = ";".join([w0[0],limText(tExt[w].split(q)[0].strip().lower())])
		for q in [" help "," works by "," induce "," indicated to "]: # function
			if isinstance(re.search(q,tExt[w]),re.Match): w0[1] = ";".join([w0[1],limText(tExt[w].split(q)[1].strip().lower())])
		for q in [" treat "," treatment of "]: # disease
			if isinstance(re.search(q,tExt[w]),re.Match): w0[2] = ";".join([w0[1],limText(tExt[w].split(q)[1].strip().lower())])
		for q in ["Dosage form"," dosage "," doses of "," dose of "]: # ingestion form
			if isinstance(re.search(q,tExt[w]),re.Match): w0[3] = ";".join([w0[3],limText(tExt[w].split(q)[1].strip().lower())])
	return w0

def sDrugs(tExt): ## 20220731
	"""sdrugs info extraction"""
	w0 = ["" for w in range(4)]
	w1 = 0
	for w in range(len(tExt)):
## drug class
		if isinstance(re.search("pharmaceutical companies",tExt[w]),re.Match): w1 = 0
		if (w1==1) and isinstance(re.search(" what ",tExt[w]),type(None)) and isinstance(re.search("defined as",tExt[w]),type(None)) and isinstance(re.search(" code",tExt[w]),type(None)) and isinstance(re.search("example",tExt[w]),type(None)): w0[0] = ";".join([w0[0],limText(tExt[w].strip().lower())])
		if isinstance(re.search("destination \| category:",tExt[w]),re.Match): w1 = 1
		for q in [" help "," helps "," works by "," induce "," indicated to "," essential for "," application of "]: # function
			if isinstance(re.search(q,tExt[w]),re.Match): w0[1] = ";".join([w0[1],limText(tExt[w].split(q)[1].strip().lower())])
		for q in [" treat "," treatment of "," caused by "," used in "]: # disease
			if isinstance(re.search(q,tExt[w]),re.Match): w0[2] = ";".join([w0[1],limText(tExt[w].split(q)[1].strip().lower())])
		for q in ["Dosage form"," dosage "," doses of "," dose of "," administered ","DOSAGE AND ADMINISTRATION"]: # ingestion form
			if isinstance(re.search(q,tExt[w]),re.Match): w0[3] = ";".join([w0[3],limText(tExt[w].split(q)[1].strip().lower())])
	return w0

def pIt(tExt): ## 20220728
	"""pillintrip info extraction"""
	w0 = ["" for w in range(4)]
	for w in range(len(tExt)):
		if isinstance(re.search("Therapeutic indications",tExt[w]),re.Match): # disease
			w1 = tExt[w+1].lower().split(".")
			for q in range(len(w1)):
				for q0 in [" treatment of "," treat "," helps "," caused by "]:
					if isinstance(re.search(q0,w1[q]),re.Match): w0[2] = ";".join([limText(w1[q].split(q0)[1]),w0[2]])
				if isinstance(re.search(" treatment ",w1[q]),re.Match): w0[2] = ";".join([limText(w1[q].split(" treatment ")[0]),w0[2]])
		if isinstance(re.search(" dose:",tExt[w]),re.Match): # ingestion form
			w1 = tExt[w].lower().split(".")
			for q in range(len(w1)):
				if isinstance(re.search(" dose:",w1[q]),re.Match): w0[3] = ";".join([limText(w1[q].split(":")[1].strip().lower()),w0[3]])
	return w0

def mLp(tExt): ## 20220728
	"""medlineplus info extraction"""
	w0 = ["" for w in range(4)]
	w2 = 0
	for w in range(len(tExt)):
		if isinstance(re.search("for more information",tExt[w]),re.Match): w2 = 0
		if w2 == 1: w0[1] = ";".join([w0[1],tExt[w].strip().lower()]) # drug mechanism
		if isinstance(re.search(" treat ",tExt[w]),re.Match):
			w1 = tExt[w].split(" treat ")
			w0[0] = ";".join([w0[0],limText(w1[0].strip().lower())]) # drug class
			w0[2] = ";".join([w0[2],limText(w1[1].strip().lower())]) # disease
			w2 = 1
		if isinstance(re.search(" treatment of ",tExt[w]),re.Match): w0[2] = ";".join([w0[2],tExt[w].split(" treatment of ")[1].strip().lower()]) # disease
		for q in [" available as "," comes as "]: # ingestion form
			if isinstance(re.search(q,tExt[w]),re.Match): w0[3] = ";".join([w0[3],tExt[w].split(q)[1].strip().lower()]) # ingestion form
	return w0

def drugINFO(dataBase,lIne,tExt):
	"""API for multiple webpage format info grabber"""
	if (dataBase==0):
		if isinstance(re.search("/",lIne[int(dataBase) + 1]),re.Match):
			if (lIne[int(dataBase) + 1].split("/")[0]=="pro"): t0 = drugPRO(tExt) ## drugs.com-PRO
			elif (lIne[int(dataBase) + 1].split("/")[0]=="ingredient"): t0 = drugIGD(tExt) ## drugs.com-ingredient
			elif (lIne[int(dataBase) + 1].split("/")[0]=="condition"): t0 = drugCOND(tExt) ## drugs.com-condition
			elif (lIne[int(dataBase) + 1].split("/")[0]=="cons"): t0 = drugCONS(tExt) ## drugs.com-cons
			elif (lIne[int(dataBase) + 1].split("/")[0]=="international"): t0 = drugINT(tExt) ## drugs.com-international
			elif (lIne[int(dataBase) + 1].split("/")[0]=="npc"): t0 = drugNPC(tExt) ## drugs.com-NPC
			else: t0 = drugMTM(tExt) ## drugs.com-MTM
		else: t0 = drugMTM(tExt) ## general filter on drugs.com pages
	elif (dataBase==1): t0 = mediORG(tExt) ## medicines.org.uk
	elif (dataBase==2): t0 = drugBank(tExt) ## drugbank
	elif (dataBase==3): t0 = nDrugs(tExt) ## ndrugs
	elif (dataBase==4): t0 = sDrugs(tExt) ## sdrugs
	elif (dataBase==5): t0 = pIt(tExt) ## pillintrip
	elif (dataBase==6): t0 = mLp(tExt) ## medlineplus
	for i1 in range(4): lIne[len(lIne)-4+i1] = ";".join([re.sub(",","!",t0[i1]),lIne[len(lIne)-4+i1]])
	return lIne

##### import #####
f = open(nAm + "-url.csv", 'r');
fL = f.readlines();f.close();
fW = open(nAm + "-wiki.csv", 'w');
sT = timeit.default_timer();

##### fill in details from webpages #####
for i in range(len(fL)):
	if (i>0):
		LN = fL[i].split(",");
		for i0 in range(len(u)): # collect website data
			if (LN[int(i0) + 1]!=""):
				if isinstance(re.search(";",LN[int(i0) + 1]),re.Match):
					uList = LN[int(i0) + 1].split(";")
				else: uList = [LN[int(i0) + 1]]
				for i1 in uList:
					uRl = "https://" + u[i0] + i1; # set target webpage
					try: rEq = requests.get(uRl); # double trial max collect webpage info
					except requests.exceptions.ConnectionError: rEq = requests.get(uRl);
					tX = [x for x in list(filter(None,re.sub("\.","\n",bs4.BeautifulSoup(rEq.text, "html.parser").get_text()).split("\n"))) if len(x)>=7]; # filter for webpage long text lines
					LN = drugINFO(i0,LN,tX); # gather webpage drug details

		for i0 in range(4): # filter for keywords
			a0 = Counter([x for x in re.sub("label","",re.sub(LN[0],"",re.sub(r"[;()]"," ",LN[len(LN)-4+i0]))).split() if len(x)>=5]).most_common(5);
			a1 = "";
			for i1,i2 in a0: a1 = ' '.join([a1,i1])
			LN[len(LN)-4+i0] = re.sub(" ",";",a1.strip());

		fW.write(str(','.join(LN) + "\n"));
	else:
		fW.write(fL[i]);

	t1 = timeit.default_timer()-sT;
	print('\rProgress: ' + str(round((i+1)/len(fL)*100,1)) + " %; ETA: " + str(round((t1/(i+1)*len(fL)-t1)/(60**2),1)) + " hours", end="");

print('\rProgress: 100');print("Web info collection total runtime: " + str(round(t1/(60**2),1)) + " hours");fW.close();
