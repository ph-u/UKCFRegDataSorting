#!/bin/env python3
# author: ph-u
# script: searchGoogle.py
# desc: search text in Google for reproducible text correction
# in: python3 searchGoogle.py
# out: otherSp-wiki.csv
# arg: 0
# date: 20220216

##### env #####
import bs4, requests, re, os, copy, sys;
os.chdir("../data"); sEp=",";tXb = "Taxonomy browser";
wEb = " Wikipedia '" + tXb + "' OR microbewiki OR PubMed OR NCBI OR StatPearls OR DSMZ OR GBIF";

##### Google search #####
fW = open("otherSp-wiki.csv", 'w');
with open(sys.argv[1], 'r') as f:
	fL = f.readlines();
	for i in range(len(fL)):
		qR = re.sub("\n","",fL[i]); #"geeksforgeeks";
		if (i>0):
			if (len(qY)<8):
				qY = qR + " microorganism";
			else:
				qY = copy.deepcopy(qR);
			uRl = "https://google.com/search?q="+ re.sub(r" ", "+", qY + wEb + "&lr=lang_en");
			rEq = requests.get(uRl);
			tX = bs4.BeautifulSoup(rEq.text, "html.parser").find_all('h3');
			if (len(tX)>0): hD = re.sub(sEp,"!",tX[0].getText());
			else: hD = "";
		else:
			qY = copy.deepcopy(qR);
			hD = "GoogleFirstHit" + sEp + "sourceWeb";
		if re.search(tXb,hD): hD = re.sub("\)", "", re.sub(tXb+" \(", "", re.sub(tXb+" linkout page \(","",hD)));
		fW.write(qR + sEp + re.sub(r" - ",sEp,hD,1)+'\n'); # rev: [::-1]
		if (i%(len(fL)/100)): print('\rProgress: '+str(round(i/len(fL)*100))+" %", end="");

print('\rProgress: 100%');
fW.close();
