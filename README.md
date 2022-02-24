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
medReconstruct.r (dropped) | cf425.rda, colDict-ExcSep.csv | remapCF425.rda, cf425Cols.txt | tabulate medication records (efficiency issue, replaced by dataMapCF425.r)
dataMapCF425.r | cf425.rda, colDict-ExcSep.csv | cf425FULL.rda, otherSp.csv | remap full data into one dataframe by regid and year (tabulate medication record), extract other species record
sortSeq.sh | otherSp.csv | otherSp-raw.csv | rearrange and sort out unique species names
(manual) | otherSp-raw.csv | otherSp-mod.csv | filter searchable names
searchGoogle.sh | otherSp-mod.csv | tmp | generate unique text list temporary for reproducible Google search (optimize reproducibility, traceability and efficiency trade-off)
searchGoogle.py | tmp | otherSp-wiki.csv | text correction using reproducible Google search (English language) focusing top hit from a selected list of databases
wikiSeq.sh (dropped) | otherSp-mod.csv | otherSp-wiki.csv | wikipedia search
(manual) | otherSp-wiki.csv | otherSp-qc.csv | quality check automated search result; manual evaluation used [duckduckgo](https://duckduckgo.com/) search engine
dataStd.r | cf425FULL.rda, (ongoing) | (ongoing) | tabulate presence/absence full species data
