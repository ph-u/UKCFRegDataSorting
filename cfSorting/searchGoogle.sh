#!/bin/bash
# author: ph-u
# script: searchGoogle.sh
# desc: reorganise mannual-sorted species list and Google search
# in: bash searchGoogle.sh
# out: ../data/tmp [intermediate]
# arg: 0
# date: 20220219

head -n 1 ../data/otherSp-mod.csv > ../data/tmp
tail -n +2 ../data/otherSp-mod.csv | sort | uniq >> ../data/tmp
python3 searchGoogle.py tmp
rm ../data/tmp
exit
