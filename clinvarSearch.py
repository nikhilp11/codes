'''
This Script takes a csv file with genomic coordinate.
Based on the genomic coordinates it retrives Clinvar Ids.
Using those clinvar ids, other informations are extracted from the records.
Along with this a list of words is provided to search for those words in those clinvar ids.
It generates 3 csv file:-
-First csv File with genomic coordinates with clinvar ids and the total count of words found in them
-Second csv file with sentences containing those words in each clinvar ids
-Third csv file with Clinvar information

python clinvarSearch.py --inputFile <path-to-csv-file> --word <any-word-to-search> --outPrefix <prefix-for-output>
'''

#Importing all the required librarires
import os
import sys
import re
import requests
import pandas as pd
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from ratelimit import limits, sleep_and_retry



#Creating a argeparse list
import argparse
parser=argparse.ArgumentParser(description='Getting Clinvar Information like Clinvar ID, Clinical Significance, etc.\nAnd search for all the hits mentioned in the clinvar Id.\nThere will be Three output CSV Files.')
parser.add_argument('--inputFile',help='Path of the file with genomic cooridirnates with column "coordinates", with format GRCh:chr#:start:end')
parser.add_argument('--word',help='Mention the word to search in resulted clinvar Id, comma seperated.')
parser.add_argument('--left_flanking',help="Flanking to left i.e. number of bp before start position coordinate",type=int,default=0)
parser.add_argument('--right_flanking',help="Flanking to right i.e. number of bp after end position coordinate",type=int,default=0)
parser.add_argument('--outPrefix',help='Prefix for the output Files')
parser.add_argument('--test',help='Running the script in test mode',default='False')
args=parser.parse_args()

if args.test == 'True':
    valueList = ['retin', 'eye', 'ocular', 'optom', 'cancer', 'tumor', 'cytokine', 'signaling']
    Path(f'{os.getcwd()}/clinvarSearch_Test/').mkdir(parents=True, exist_ok=True)
    outputPath = f'{os.getcwd()}/clinvarSearch_Test/Test_'
    inputFile = ['GRCh38:chr10:50555728:50555751', 'GRCh38:chr1:94027324:94027347', 'GRCh38:chr11:67459747:67459770']
    #Reading the csv file with Genomic coordinates
    targetDf = pd.DataFrame(inputFile,columns=['coordinates'])
    leftFlank = 200
    rightFlank = 200
else:
    #Reading the csv file with Genomic coordinates
    inputFile = args.inputFile
    targetDf = pd.read_csv(inputFile,usecols=['coordinates'])
    leftFlank = args.left_flanking
    rightFlank =args.right_flanking
    #Default key words to search in clinvar id
    valueList=[]
    if args.word == None:
        print('Please enter some words to search. Use --word <words-to-search>')
        sys.exit()
    elif ',' in args.word:
        valueList = args.word.split(",")
    else:
        valueList.append(args.word)

    #A prefix for the output or takes the name from input file
    if args.outPrefix:
        outputPath=f'{Path(args.inputFile).parent}/{args.outPrefix}_'
    else:
        outputPath=f'{Path(args.inputFile).parent}/{Path(args.inputFile).stem}_'


print('####################################################################\n\n')
print(f'Input File: - {inputFile}')
print(f'Words To Search: - {valueList}')
print(f'Flanking to the left (before start coordinate): - {leftFlank}')
print(f'Flanking to the right (after end coordinate): - {rightFlank}')
print(f'Output Prefix with Path: - {outputPath}\n\n')
print('####################################################################\n\n')


#Function to get the clinvar VCV id from the ncbi eutils using genomic coordinates
#The query is the search term based on genomic coordinates and filtering results based on clinical significance and short variants <50bp
#We get result in JSON format and it is stored in dataframe in list format
# Define the rate limit: 3 calls per second
CALLS = 3
PERIOD = 1

@sleep_and_retry
@limits(calls=CALLS, period=PERIOD)
def getVcvID(df,leftFlank=leftFlank,rightFlank=rightFlank):
    chrom=df.coordinates.split(":")[1].replace("chr","")
    start38=int(df.coordinates.split(":")[2]) - leftFlank
    end38=int(df.coordinates.split(":")[3]) + rightFlank
    query = f'''({chrom}[CHR] AND 1:{end38}[chrpos38] AND {start38}:2000000000[chrpos38] AND ( ( ("clinsig pathogenic"[Properties] or "clinsig pathogenic low penetrance"[Properties] or "clinsig established risk allele"[Properties]) OR ("clinsig likely pathogenic"[Properties] or "clinsig likely pathogenic low penetrance"[Properties] or "clinsig likely risk allele"[Properties]) OR ("clinsig vus"[Properties] or "clinsig uncertain risk allele"[Properties]) ) AND 0[VARLEN]:49[VARLEN]))'''
    esearchLink =f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term={query}&rettype=uilist&retmode=json'
    decode = (requests.get(esearchLink)).json()
    if 'esearchresult' in decode and 'idlist' in decode['esearchresult']:
        return decode['esearchresult']['idlist']
    return []  # Return an empty list if all retries fail
    

'''
Using those clinvar Id, request is made to get the page as a text by removing all the html tags
A list of sentences containing all the words is created
'''
def getText(clinvarId,valueList):
    headers = {'user-agent': 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/86.0.4240.111 Safari/537.36'}
    url = f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{clinvarId}"
    response = requests.get(url, headers=headers)
    text = re.sub('<[^<]+?>', '', response.text).strip()
    searchRes = [*set(re.sub(r'^[^a-zA-Z0-9]+', '', item) for item in text.split("\n") for word in valueList if word.lower() in item.lower())]
    return searchRes

'''
This function returns the method used for the submission of the clinical interpretation
It scrapes the data from the url'''
def rcvCheck(rcvId):
    url = f"https://www.ncbi.nlm.nih.gov/clinvar/{rcvId}/"
    response = requests.get(url)
    html = response.content
    df = pd.read_html(html)[0]['Method'].value_counts().reset_index()
    method = df.assign(value = df["Method"] + " (" + df["count"].astype(str) + ")")["value"][0]
    return method

'''
Using the clinvar ids we make a request using eutils to get the summary records of that clinvar id
We get the data on JSON format and using that we extract required information
It returns a dictionary
'''
def esumInfo(vcvId):
    esumLink=f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=clinvar&id={vcvId}&retmode=json&rettype=clinvarset'
    esumRes = requests.get(esumLink)
    res=esumRes.json()
    if 'result' in res and str(vcvId) in res['result']:
        decode = res['result'][str(vcvId)]
        rcvList = [rcvCheck(rcvId) for rcvId in decode["supporting_submissions"]["rcv"]]
        esumDict = {'clinsig' : decode['germline_classification']['description'],
                    'review_status' : decode['germline_classification']['review_status'],
                    'gene' : decode['gene_sort'],
                    'trait' : [decode['germline_classification']['trait_set'][x]['trait_name'] for x in range(len(decode['germline_classification']['trait_set']))],
                    'variation_name' : [decode['variation_set'][x]['variation_name'] for x in range(len(decode['variation_set']))],
                    'canonical_spdi' : [decode['variation_set'][x]['canonical_spdi'] for x in range(len(decode['variation_set']))],
                    'protein_change' : decode['protein_change'],
                    'method' : ", ".join(rcvList)}
        return esumDict
    else:
        print(vcvId, "empty")
    

#Apply the getVcvId function to get the clinvar ID
# Helper function to apply getVcvID in parallel
def parallel_get_vcv_id(row):
    vcvIdDf = getVcvID(row, leftFlank, rightFlank)
    return vcvIdDf

# Apply the getVcvID function to get the clinvar ID
with ThreadPoolExecutor() as executor:
    targetDf['VcvIds'] = list(executor.map(parallel_get_vcv_id, [row for _, row in targetDf.iterrows()]))



#Getting all the vcv id in a list format
vcvIdList = [*set(sum(targetDf.VcvIds.to_list(),[]))] # type: ignore
print('Retrieved VCV Clinvar ID from the targets')

#Creating a dictionary of those clinvar id where Id is the key and all the sentences in list format as value
clinvarInfoDict = {clinvarId: getText(clinvarId,valueList) for clinvarId in vcvIdList}


#Creating a dataframe of clinvar id and sentences
infoList = []
for key, values in clinvarInfoDict.items():
    for value in values:
        infoList.append([key, value])
infoDf = pd.DataFrame(infoList, columns=['ID', 'Info'])
print('Info Df generated based on VCV Id and words provided')


# Create the DataFrame from the list of dictionaries containing summary of clinvar id
data_list = []
for key, values in {vcvId:esumInfo(vcvId) for vcvId in vcvIdList}.items():
    data_row = {'Vcv_ID': key}
    data_row.update(values)
    data_list.append(data_row)
esumDf = pd.DataFrame(data_list)
print('Clinvar Info extracted from the VCV ID using Esummary')

#Creating a new column counting total number of hits of each word in the clinvar Id
targetDf['TotalCount'] = targetDf['VcvIds'].apply(lambda x: pd.Series(x,dtype=pd.StringDtype()).map({clinvarId:len(clinvarInfoDict[clinvarId]) for clinvarId in vcvIdList}).tolist())

#Exporting all the csv files
esumDf.to_csv(f"{outputPath}VCV_ID_info.csv",index=False)
targetDf.to_csv(f"{outputPath}Location-VCV_ID_info.csv",index=False)
infoDf.to_csv(f"{outputPath}VCV_ID-WordSearch.csv",index=False)
print('Exported all the dataframes')
