# codes


## ClinVar Search
ClinVar Search helps extracts data from the ClinVar page using EUtils package. To get all the short variants i.e. Variants < 50kb (SNV and Short INDELS). 
Make sure to pass a csv file with column name 'coordinates' and the genomic coordinates in GRCh38:chromosome:start:end format.
Packages Required:
* pandas
* ratelimit

`python clinvarSearch.py --inputFile <path-to-csv-file> --word <any-word-to-search> --outPrefix <prefix-for-output>`