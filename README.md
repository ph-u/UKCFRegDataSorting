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
(manual) | drug.csv | drug-qcList.csv | standardize medical classes and details
reArrange.r | cf425FULL.rda, {otherSp,medNamegenCols}-qcList.csv | cf425MedMic.rda, {otherSp-qcREF,cf425Medic,cf425Micro}.csv | rearrange columns to medical,microbe dataframes (data sort log: Rscript reArrange.r >> ../data/reArrangeRec.txt; 0.5hr runtime)
selCheck.r | cf425MedMic.rda | NA | select 30 manual check data rows in cf425
micro2Genus.r | cf425Micro.csv / cf425MedMic.rda | genusCF425\_gLV.csv, cf425Genus.pdf | sort mIcro data in genus time-series
genusTSCut.r | genusTimeSeries\_gLV.csv | gTS\_{startYr}{endYr}-gLV.csv |  cut genus time-series data in multiple csv
plotGenusTS.r | gTS_{0820,0811,0813,1015,1113,1315,1619}\_gLV-{log,sam}.csv | gTS\_overall.pdf | plot time-series genus data with fitted simulations
