import pandas
import glob
files=glob.glob("./../Longitudinal_Data/*.csv")

def label(x):
	return x.split("_")[1]

for i in files:
	df=pandas.read_csv(i,index_col=0)
	df["label"]=df.index.map(label)
	df=df[(df.label=="BSL")|(df.label=="EX1")|(df.label=="P1")]
	#Remove patient number starting with 462 since they are not present in all three BSL EX1 and P1.
	df=df.drop(["462_BSL","462_EX1"])
	df=df.drop(["label"],axis=1)
	#Normalization
	df=df.div(df.sum(axis=1), axis=0)
	df=df.multiply(100, fill_value=0)
	#Filtering Non-zero value in atleast 5% samples 
	nsel=df>0	
	nsel=nsel.sum(axis=0)
	y=nsel[nsel<3].index #In 5% of the patients
	df.drop(columns=y,inplace=True)
	#Re Normalization
	df=df.div(df.sum(axis=1), axis=0)
	df=df.multiply(100, fill_value=0)
	df["label"]=df.index.map(label)
	df.to_csv("./Pre_processed_data/"+str(i)[:-4]+"_mod.csv")	
