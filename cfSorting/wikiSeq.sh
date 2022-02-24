#!/bin/bash
# author: ph-u
# script: sortSeq.sh
# desc: rearrange other species entries in single column
# in: bash sortSeq.sh
# out: ../data/otherSp-wiki.csv
# arg: 0
# date: 20220215

cd ../data
[[ -f otherSp-wiki.csv ]] && rm otherSp-wiki.csv

##### naive wiki search #####
while read -r i;do
	if [[ ${i} -eq "input" ]];then
		i2="WikipediaTopHit"
	else
		#echo -e ${i}
##### set wiki search text string #####
		if [[ `echo -e "${i}" | grep -o " " | wc -l` -lt 1 ]] && [[ `echo -e "${i}" | wc -c` -lt 7 ]];then
			i1=" microorganism"
		else
			i1=""
		fi
		i0=`echo -e "${i}${i1}" | sed -e "s/ /+/g"`
##### wiki search #####
		i2=`curl "https://en.wikipedia.org/w/index.php?search=${i0}&title=Special%3ASearch&go=Go&ns0=1" | grep -e "The page" | sed -e "s/class/\n/g" | grep -e "search-result-heading" | sed -e "s/title/\ntitle/g" | grep -e "title" | head -n 1 | cut -f 2 -d "\""`
	fi
	echo -e "${i},${i2}" >> otherSp-wiki.csv
done < otherSp-mod.csv

exit
