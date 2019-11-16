from __future__ import division
import pandas

bac=pandas.read_csv("./../../../Data/bacteria.csv",index_col=0)
fun=pandas.read_csv("./../../../Data/fungi.csv",index_col=0)
vir=pandas.read_csv("./../../../Data/virus.csv",index_col=0)
lab=pandas.read_csv("./../../../wSNF_Analysis/Merged_biome_clustering/results/labels.csv",index_col=0)

vir=vir.div(vir.sum(axis=1), axis=0)
vir=vir.multiply(100, fill_value=0) #Normalize

x=bac.join(fun,how="inner")
y=x.join(vir,how="inner")

y=y.div(y.sum(axis=1), axis=0) #Renormalize

#Dropping 0 columns
nsel=y.sum(axis=0)
n=nsel[nsel==0].index
y.drop(columns=n,inplace=True)

#Non-zero value in atleast 10 patients
nsel=y>0.01
nsel=nsel.sum(axis=0)
x=nsel[nsel<10].index #5%
y.drop(columns=x,inplace=True)

y=y.div(y.sum(axis=1), axis=0) #Renormalize

df=y.join(lab,how="inner")

df.rename(columns={"labels":"x"},inplace=True)
df.to_csv("Microbes.csv") 
#Add PatientID as index column name manually
