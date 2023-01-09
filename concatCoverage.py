'''The follwoing script is to concat all the total coverage and their data size from a batch directly from the server'''
import pandas as pd
from pathlib import Path
import glob
import os, sys
import json


batch = sys.argv[1]
#Path to all the coverage files and raw data files in seperate variables
statsPath = f"/mnt/NGS1/WES_Analysis/WES_ProcessedData/Batch_{batch}/"
sizePath = f'/mnt/NGS1/WES_Analysis/RawData/Batch{batch}/'


#Check all the subfolders in the above path and get the .csv files
coverageFiles = glob.glob(os.path.join(statsPath , "**/Coverage/*_TotalCov.txt"))
hsMetricsFiles = glob.glob(os.path.join(statsPath, "**/**/*_hsMetrics.txt"))

#empty list to concat all the coverage dfs in the batch
li = []
#looping through all those files and read them into a list
for covFile in coverageFiles:
    df = pd.read_csv(covFile,sep="\t",low_memory=False).set_index("Depth").T
    df.insert(0,'SampleName',(Path(covFile).stem).rsplit("_",1)[0])
    df.insert(1,'Raw Data Size (GiB)',"")
    df.drop("Difference",axis = 1,inplace=True)
    li.append(df)

#creating a list of all the raw data file path in the provided batch
sizesFiles = glob.glob(os.path.join(sizePath, "*.fq.gz"))

#Dictionary comprehension to create a dictionary {SampleName : Data Size}
sampleSize = {Path(file).stem: round(os.path.getsize(file)/1024**3,7) for file in sizesFiles}

#Creating df of the above dictionary so that we can groupby the sample name to sum the data size of forward and reverse files
sizeDf = pd.DataFrame(sampleSize.items(),columns=['Sample','DataSize'])
sizeDf['Sample'] = sizeDf.Sample.str.replace("_[1|2].fq","",regex=True)
sizeDf = (sizeDf.groupby("Sample",as_index=False)["DataSize"].sum()).round(2)


#concating just the last row from all the df saved in the concated df list
statsDf = pd.concat([d.tail(1) for d in li], ignore_index=True)
statsDf.rename(columns = {"#chrom" : "Coverage"},inplace = True)
#mapping the data size based on sample name after grouping the size df
statsDf['Raw Data Size (GiB)'] = statsDf.SampleName.map(dict(zip(sizeDf.Sample, sizeDf['DataSize'])))

#Dictionary Comprehension to load sample name as key and their respective Dataframe as value and then extracting just the Fold80 value
hsMetricsDfDict = {Path(file).stem.rsplit("_",1)[0]: pd.read_csv(file,sep='\t',comment="#",low_memory=False) for file in hsMetricsFiles}
fold80Df = {Path(file).stem: hsMetricsDfDict[file]['FOLD_80_BASE_PENALTY'][0] for file in hsMetricsDfDict.keys()}

#Mapping Fold80 Value to original df
statsDf['Fold_80'] = statsDf.SampleName.map(fold80Df)

preRunStatsDf = pd.read_excel(f"/mnt/NGS1/WES_Analysis/WES_ProcessedData/Batch_{batch}/Batch_{batch}_Stats.xlsx")

#Merging the pre run statistics file with the post stats file
finalDf = pd.merge(preRunStatsDf,statsDf,how="left",on="SampleName")

finalDf = finalDf.sort_values('TotalBases (G)')

finalDf = finalDf.drop_duplicates(subset="SampleName")

l =['Total','',finalDf["TotalBases (G)"].sum().round(3),finalDf["TotalReads (M)"].sum().round(3),'','','','','','','','','','','','','','','','','','','']

filteredCols =['SampleName','Raw Data Size (GiB)','Q30 (%)','TotalBases (G)','TotalReads (M)','DuplicationRates (%)','Read Length','Sex Check','PerBaseQuality_1+2 (Raw)','PerBaseQuality_1+2 (Trimmed)','Fold_80','1X','5X','10X','20X','30X','40X','50X','60X','70X','80X','90X','100X']

finalDf.loc[len(finalDf)] = l
finalDf.to_excel(f"/mnt/NGS1/WES_Analysis/WES_ProcessedData/Batch_{batch}/Batch_{batch}_Stats.xlsx",index=False,columns=filteredCols)
