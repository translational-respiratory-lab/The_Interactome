from __future__ import division
import pandas

bac=pandas.read_csv("./../../../Data/1_Targeted/bacteria.csv",index_col=0)
fun=pandas.read_csv("./../../../Data/1_Targeted/fungi.csv",index_col=0)
vir=pandas.read_csv("./../../../Data/1_Targeted/virus.csv",index_col=0)
lab=pandas.read_csv("./../../../wSNF_Analysis/Merged_biome_clustering/results/labels.csv",index_col=0)

vir=vir.div(vir.sum(axis=1), axis=0)
vir=vir.multiply(100, fill_value=0) #Normalize

x=bac.join(fun,how="inner")
y=x.join(vir,how="inner")

y=y.div(y.sum(axis=1), axis=0) #Renormalize

df=y.join(lab,how="inner")

#Dropping 0 columns
nsel=df.sum(axis=0)
n=nsel[nsel==0].index
df.drop(columns=n,inplace=True)

'''
#Non-zero value in atleast 5 patients
nsel=df>0.01
nsel=nsel.sum(axis=0)
y=nsel[nsel<5].index #5 or greater than 5
df.drop(columns=y,inplace=True)
'''
df.rename(columns={"labels":"x"},inplace=True)
df.to_csv("Microbes.csv")
#Manually add PatientID as the index column name
