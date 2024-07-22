# codes


## ClinVar Search
ClinVar Search helps extracts data from the ClinVar page using EUtils package. To get all the short variants i.e. Variants < 50kb (SNV and Short INDELS). 
Make sure to pass a csv file with column name 'coordinates' and the genomic coordinates in GRCh38:chromosome:start:end format.
Packages Required:
* pandas
* ratelimit

```python clinvarSearch.py --inputFile <path-to-csv-file> --word <any-word-to-search> --outPrefix <prefix-for-output>```

This code uses `concurrent.futures` for simulateneously get the ClinVar ID for multiple Genomic Coordinates at once, and since EUtils has api call limit that is 3 query per second so `ratelimit` package is used to control that flow.