#!/bin/bash
# author: ph-u
# script: sortSeq.sh
# desc: rearrange other species entries in single column
# in: bash sortSeq.sh
# out: ../data/otherSp-raw.csv
# arg: 0
# date: 20220215

cd ../data
[[ -f otherSp-raw.csv ]] && rm otherSp-raw.csv

cut -f 2- -d "," otherSp.csv | tail -n +2 |\
sed -e "s/\&/,/g" |\
sed -e "s/+ /,/g" |\
sed -e "s/ x /,/g" |\
sed -e "s/ and /,/g" |\
sed -e "s/\/\//,/g" |\
sed -e "s/_x00d_//g" | sed -e "s/_x000d_//g" | sed -e "s/10/0/g" |\
sed -e "s/,/\n/g" | sed -e "s/^ //" | sed -e "s/ $//" | sed -e'/^\r$/d' |\
sort | uniq | grep -vwE '\w{1,1}'> otherSp-raw.csv

exit
