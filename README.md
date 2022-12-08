# UK Cystic Fibrosis Registry data (ref \#: 425) handling script sequence

- Raw data files (xlsx, doc, docx, pdf) in CF Registry data directory
- I/O data files (rda, csv) in `data` directory
- Figures (pdf) in thesis figure directory

script | in | out | desc
--- | --- | --- | ---
rearrangeCF425.r | Data (xlsx) | cf425.rda | convert raw data to RData (0.5hr runtime)
describeCF425.r | cf425.rda | tsPatients.pdf, attributesTotal.csv | data time-series descriptions (time-series plot, overall column data structure)
(manual) | Data guides (doc,docx,pdf) | cf425\_extraDataGuide.csv | reference R-readable column names source \& meaning from given documents
colRefCF425.r | Data dictionary (xlsx), attributesTotal.csv, cf425\_extraDataGuide.csv | colDict-ExcSep.csv | generate column names equivalence table for UKCF-Reg data
dataMapCF425.r | cf425.rda | cf425FULL.rda, genCols.csv | remap full data into one dataframe by regid and year (tabulate medication record)
otherSp.r | cf425FULL.rda | otherSp.csv | extract other species record
medName.r | cf425FULL.rda | medName.csv | extract medicine name
(manual) | otherSp.csv | otherSp-mod.csv | filter searchable names: include antimicrobial resistence
searchAuto.sh | {otherSp-mod,genCols,medName}.csv | tmp | generate unique text list temporary for Wikipedia search (optimize reproducibility, traceability and efficiency trade-off)
searchAuto.py | tmp | {otherSp,genCols,medName}-wiki.csv | text correction using mostly reproducible Wikipedia search top hit
(manual) | {otherSp,genCols,medName}-wiki.csv | {otherSp,genCols,medName}-qcList.csv | quality check automated search result; manual evaluation used [duckduckgo](https://duckduckgo.com/) search engine
otherSpPlots.r | otherSp-qcList.csv | otherSp-qcREF.csv, GoogleSearchEfficiency.pdf, DataIrregularity.pdf | (construct a reference frame for other species), summarise efficiency and effort for the manual text correction process
medNameExtraction.r | medName-qcList.csv | drug.csv | get standardized medication and active ingredient names
(manual) | drug.csv | drug-url.csv | standardize medical classes and details
medNameExtraction.py | drug-url.csv | drug-wiki.csv | collect standardized info
reArrange.r | cf425FULL.rda, {otherSp,medName,genCols}-qcList.csv | cf425MedMic.rda, {otherSp-qcREF,cf425Medic,cf425Micro}.csv | rearrange columns to medical,microbe dataframes (data sort log: Rscript reArrange.r >> ../data/reArrangeRec.txt; 0.5hr runtime)
selCheck.r | cf425MedMic.rda | NA | select 30 manual check data rows in cf425
mutSort.r | [1st run] 425A\_DataRequest\_final.xlsx, mutation\_lookup\_table\_April2022.xlsx [others] mutLib\_r.csv, mutLib.csv | [1st run] mutLib.{csv,txt} [others] mutPWCF.csv | sort CF mutation of pwCF
micro2Genus.r | cf425Micro.csv / cf425MedMic.rda | genusCF425\_gLV.csv, cf425Genus.pdf | sort mIcro data in genus time-series
genusTSCut.r | genusTimeSeries\_gLV.csv | gTS\_{startYr}{endYr}-gLV.csv |  cut genus time-series data in multiple csv
cftrM\_sep.r | cf425MedMic.rda | ../cftrM\_raw/cftrM\_{abs,mod,int}\_\[start\]\[end\]\_gLV.csv | categorize people with CF (pwCF) into CFTR-modifiers (cftrM) and interacting drugs categories + separate into different time-series
plotGenusTS.r | gTS_{0820,0811,0813,1015,1113,1315,1619}\_gLV-{log,sam}.csv | gTS\_overall.pdf | plot time-series genus data with fitted simulations
clusterHeatmap.r | cftrM\_data/\*.csv | graph/cftrM/heatmap/heatMap_\*.pdf | heatmap clustering grphs for CFTRm segregated time-series data
F508delTSPlots.r | F508del\_data/\*.csv | graph/F508del/\*.pdf | rolling mean plots of 3 consecutive years for F508del mutation CFTRm segregated time-series
F508delStats.r | [selected -eco.csv] | [none] | statistical tests for F508del group
F508delCorrelation.r | cf425MedMic.rda, cftrm\_interaction.csv, antiInfectives.csv, mutPWCF.csv, F508del\_data/\*-eco.csv | F508\_Correlation.csv | Correlation between CFTRm yearly diff with gLV model ability

Notes for Computing

- 20 mins for each prior search (4.5 CPU hours for 293 time-series; 20 Earth mins)
- 7-11 hours for each BI-MCMC run (19,000 CPU hours for 293 time-series with 7 replicates each, 2051 runs in total partitioned into max limit of 448 simultaneous jobs in CSD3 cluster; 47 Earth hours)
- a few seconds for each ecology analysis run (0.0 CPU hours recorded; a few Earth seconds)
