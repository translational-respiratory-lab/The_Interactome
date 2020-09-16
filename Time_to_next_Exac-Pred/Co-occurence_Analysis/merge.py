from __future__ import division
import pandas

i="P1"

bac=pandas.read_csv("./../Pre-processing/Pre-processed/bac_"+i+".csv",index_col=0)
fun=pandas.read_csv("./../Pre-processing/Pre-processed/fungi_"+i+".csv",index_col=0)
vir=pandas.read_csv("./../Pre-processing/Pre-processed/vir_"+i+".csv",index_col=0)

lab=vir.loc[:,"Time to next exacerbation "]
bac.drop(["Time to next exacerbation "],axis=1,inplace=True)
fun.drop(["Time to next exacerbation "],axis=1,inplace=True)
vir.drop(["Time to next exacerbation "],axis=1,inplace=True)

vir=vir.div(vir.sum(axis=1), axis=0)
vir=vir.multiply(100, fill_value=0) #Normalize

bac=bac.div(bac.sum(axis=1), axis=0)
bac=bac.multiply(100, fill_value=0) #Normalize

fun=fun.div(fun.sum(axis=1), axis=0)
fun=fun.multiply(100, fill_value=0) #Normalize


x=bac.join(fun,how="inner")
y=x.join(vir,how="inner")

df=y.div(y.sum(axis=1), axis=0) #Renormalize

#Dropping 0 columns
nsel=df.sum(axis=0)
n=nsel[nsel==0].index
df.drop(columns=n,inplace=True)

#Non-zero value in atleast 10 samples
nsel=df>0.01   #atleast 0.01 out of 100
nsel=nsel.sum(axis=0)
y=nsel[nsel<4].index #In 5% of the patients to add in methods
df.drop(columns=y,inplace=True)

df=df.div(df.sum(axis=1), axis=0) #Normalize

lab[lab<12*7]=1 #high-risk
lab[lab>12*7]=2 #low-risk

df=pandas.merge(df,lab,left_index=True,right_index=True)

df1=df[df["Time to next exacerbation "]==1]
df1.drop(["Time to next exacerbation "],axis=1,inplace=True)
df2=df[df["Time to next exacerbation "]==2]
df2.drop(["Time to next exacerbation "],axis=1,inplace=True)

df2.to_csv("Microbes.csv")
