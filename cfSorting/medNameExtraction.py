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
import bs4, requests, re, os, copy, sys, timeit;
from collections import Counter;
os.chdir("../data"); sEp=",";
u = ["www.drugs.com/","www.medicines.org.uk/emc/product/","go.drugbank.com/drugs/","www.ndrugs.com/?s=","pillintrip.com/medicine/","medlineplus.gov/druginfo/meds/"]

##### func #####
def limText(sTr): ## 20220729
	"""limit text length to max 6 words"""
	s1 = sTr.split()
	if len(s1)>6: s0 = [s1[0],s1[1],s1[2],s1[len(s1)-3],s1[len(s1)-2],s1[len(s1)-1]];
	else: s0 = s1;
	return ' '.join(s0)

def drugMTM(tExt): ## 20220728
	"""drugs.com mtm/other page info extraction"""
	w0 = ["" for w in range(4)]
	for w in range(len(tExt)):
		if isinstance(re.search("Drug class",tExt[w]),re.Match): w0[0] = ';'.join([w0[0],limText(tExt[w].split(":")[1].strip().lower())])
		if isinstance(re.search("Dosage form",tExt[w]),re.Match): w0[3] = ';'.join([w0[3],limText(tExt[w].split(":")[1].strip().lower())])
		if isinstance(re.search("ask ",tExt[w].lower()),type(None)) and isinstance(re.search(" conditions ",tExt[w].lower()),type(None)):
			for q in [" produce "," can "," helps "," help "]:
				if isinstance(re.search(q,tExt[w]),re.Match): w0[1] = ';'.join([w0[1],limText(tExt[w].split(q)[1].strip().lower())])
			for q in [" treat "," lead to "]:
				if isinstance(re.search(q,tExt[w]),re.Match): w0[2] = ';'.join([w0[2],limText(tExt[w].split(q)[1].strip().lower())])
	return w0

def drugPRO(tExt): ## 20200728
	"""drugs.com pro page info extraction"""
	w0 = ["" for w in range(4)]
	for w in range(len(tExt)):
		if isinstance(re.search("Drug class",tExt[w]),re.Match): w0[0] = ";".join([w0[0],limText(tExt[w].split(":")[1].strip().lower())])
		if isinstance(re.search("Dosage form",tExt[w]),re.Match): w0[3] = ";".join([w0[3],limText(tExt[w].split(":")[1].strip().lower())])
		if isinstance(re.search("Indications ",tExt[w]),re.Match) and isinstance(re.search(" for ",tExt[w]),re.Match): w0[1] = w0[2] = ";".join([w0[1],limText(tExt[w+1].strip().lower())])
	return w0

def drugIGD(tExt):
	"""drugs.com ingredient page info extraction"""
	w0 = ["" for w in range(4)]
#	for w in range(len(tExt)):
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

def nDrugs(tExt): ## 20220729
	"""ndrugs info extraction"""
	w0 = ["" for w in range(4)]
#	for w in range(len(tExt)):
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
		if isinstance(re.search(" comes as ",tExt[w]),re.Match): w0[3] = ";".join([w0[3],tExt[w].split(" comes as ")[1].strip().lower()]) # ingestion form
	return w0

def drugINFO(dataBase,lIne,tExt):
	"""API for multiple webpage format info grabber"""
	if (dataBase==0) and (lIne[int(dataBase) + 1].split("/")[0]=="pro"): t0 = drugPRO(tX) ## drugs.com-PRO
	elif (dataBase==0) and (lIne[int(dataBase) + 1].split("/")[0]=="ingredient"): t0 = drugIGD(tX) ## drugs.com-ingredient
	elif (dataBase==0): t0 = drugMTM(tX) ## drugs.com-MTM
	elif (dataBase==1): t0 = mediORG(tX) ## medicines.org.uk
	elif (dataBase==2): t0 = drugBank(tX) ## drugbank
	elif (dataBase==3): t0 = nDrugs(tX) ## ndrugs
	elif (dataBase==4): t0 = pIt(tX) ## pillintrip
	elif (dataBase==5): t0 = mLp(tX) ## medlineplus
	for i1 in range(4): lIne[i1+9] = ";".join([re.sub(",","!",t0[i1]),lIne[i1+9]])
	return lIne

##### import #####
f = open(nAm + "-url.csv", 'r')
fL = f.readlines();f.close();
fW = open(nAm + "-wiki.csv", 'w');
sT = timeit.default_timer()

##### fill in details from webpages #####
for i in range(len(fL)):
	if (i>0):
		LN = fL[i].split(",");
		for i0 in range(len(u)): # collect website data
			if (LN[int(i0) + 1]!=""):
				uRl = "https://" + u[i0] + LN[int(i0) + 1]; # set target webpage
				try: rEq = requests.get(uRl); # double trial max collect webpage info
				except requests.exceptions.ConnectionError: rEq = requests.get(uRl);
				tX = [x for x in list(filter(None,re.sub("\.","\n",bs4.BeautifulSoup(rEq.text, "html.parser").get_text()).split("\n"))) if len(x)>=7]; # filter for webpage long text lines
				LN = drugINFO(i0,LN,tX); # gather webpage drug details

		for i0 in range(4): # filter for keywords
			a0 = Counter([x for x in re.sub(r"[;()]"," ",LN[i0+9]).split() if len(x)>=5]).most_common(5)
			a1 = "";
			for i1,i2 in a0: a1 = ' '.join([a1,i1])
			LN[i0+9] = re.sub(" ",";",a1.strip())

		fW.write(','.join(LN));
	else:
		fW.write(fL);
	if (i%(len(fL)/100)): print('\rProgress: '+str(round(i/len(fL)*100))+" %", end="");

print('\rProgress: 100');print("Web info collection total runtime: ", timeit.default_timer()-sT);fW.close();
