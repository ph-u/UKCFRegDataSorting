#!/bin/bash
# author: ph-u
# script: searchAuto.sh
# desc: reorganise mannual-sorted species list and Wikipedia search
# in: bash searchAuto.sh [input_csv_basename]
# out: ../data/tmp [intermediate]
# arg: 0
# date: 20220219

[[ -z $1 ]] && grep -e "^\# desc:\|^\# in:" $0 | cut -f 2 -d ":" | sed -e "s/^ //" && exit
a=$1
head -n 1 ../data/${a}.csv > ../data/tmp
tail -n +2 ../data/${a}.csv | sort | uniq >> ../data/tmp

##### filter only lines with >1 letter #####
#echo -e "`date`: sort unique searchable names"
#while read -r L;do
#	if [[ `echo -e ${L} | wc -c` -gt 3 ]];then echo -e ${L} >> ../data/tmp0;fi
#done < ../data/tmp
#mv ../data/tmp0 ../data/tmp

b=`echo -e ${a} | cut -f 1 -d "-"`
echo -e "`date`: auto-searches"
python3 searchAuto.py tmp ${b}
rm ../data/tmp
echo -e "`date`: auto-searches done"
exit
