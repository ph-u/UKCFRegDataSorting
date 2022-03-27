#!/bin/env python3
# author: ph-u
# script: searchAuto.py
# desc: search text in Wikipedia for mostly reproducible text correction
# in: python3 searchAuto.py [source_file] [basename]
# out: otherSp-wiki.csv
# arg: 0
# date: 20220216

tK = "" #"wiki"
##### env #####
if (tK=="wiki"):
	import wikipedia, re, os, sys;
	nT = type(None);
else:
	import bs4, requests, re, os, copy, sys, time;
	y = sys.argv[2];
	if (y.find("otherSp")==0):
		tXb = "Taxonomy browser";
		wEb = " Wikipedia '" + tXb + "' OR microbewiki OR PubMed OR NCBI OR StatPearls OR DSMZ OR GBIF";
	else: wEb = " Wikipedia medicine";
os.chdir("../data"); sEp=",";

##### rate-limited download ##### 20220324
#def gURL(uRl):
#	"""get url response at rate-limited manner"""
#	rQ = requests.request('GET', uRl, stream=True)
#	for data in rQ.iter_content(chunk_size=3096):
#		time.sleep(.001)
#
#	return rQ

##### automated search ##### 20220317 (necessary slow network)
fW = open(sys.argv[2] + "-wiki.csv", 'w');
with open(sys.argv[1], 'r') as f:
	fL = f.readlines();
	for i in range(len(fL)):
		qR = re.sub("\n","",fL[i]); #"geeksforgeeks";
		if (i>0):
			if (tK=="wiki"):
				rEs = wikipedia.suggest(qR);
				if (isinstance(rEs,nT)): rEs = "";
				hD = rEs + sEp + "Wikipedia";
#			rE = wikipedia.search(qR, results = 1);
#			if (len(rE)==0):
#				try: # https://stackoverflow.com/questions/64423312/disambiguationerror-and-guessedatparserwarning-in-wikipedia-api-in-python
#					rE = wikipedia.summary(qR).split()[0:2];
#					rEs = rE[0] + " " + rEs[1]
#				except wikipedia.exceptions.DisambiguationError as e: rEs = e.options[0];
#				except wikipedia.exceptions.PageError: rEs = "";
#			else: rEs = rE[0];
			else:
				if (len(qY)<8):
					qY = qR + " microorganism";
				else:
					qY = copy.deepcopy(qR);
				uRl = "https://google.com/search?q="+ re.sub(r" ", "+", qY + wEb + "&lr=lang_en");
#				uRl = "https://uk.search.yahoo.com/search?p="+ re.sub(r" ", "+", qY + wEb);
# limit bandwidth: https://stackoverflow.com/questions/17691231/how-to-limit-download-rate-of-http-requests-in-requests-python-library
#				try: rEq = gURL(uRl);
				try: rEq = requests.get(uRl);
				except requests.exceptions.ConnectionError: rEq = requests.get(uRl); # second try
				tX = bs4.BeautifulSoup(rEq.text, "html.parser").find_all('h3');
				if (len(tX)>0): hD = re.sub(sEp,"!",tX[0].getText());
				else: hD = "";
		else:
			if (tK=="wiki"): hD = "WikiFirstHit" + sEp + "sourceWeb";
			else:
				qY = copy.deepcopy(qR);
				hD = "GoogleFirstHit" + sEp + "sourceWeb";
		if (tK!="wiki") and (sys.argv[2]=="otherSp") and re.search(tXb,hD): hD = re.sub("\)", "", re.sub(tXb+" \(", "", re.sub(tXb+" linkout page \(","",hD)));
		fW.write(qR + sEp + re.sub(r" - ",sEp,hD,1)+'\n'); # rev: [::-1]
		if (i%(len(fL)/100)): print('\rProgress: '+str(round(i/len(fL)*100))+" %", end="");
		#if (tK!="wiki") and (i>0): time.sleep(.7) # protect from IP block (useless protection, require a slow network instead)

print('\rProgress: 100');
fW.close();
