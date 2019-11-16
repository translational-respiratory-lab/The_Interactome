from __future__ import division
import pandas

bac=pandas.read_csv("./../Longitudinal_Data/1_Bacteria.csv",index_col=0)
fun=pandas.read_csv("./../Longitudinal_Data/3_Fungi.csv",index_col=0)
vir=pandas.read_csv("./../Longitudinal_Data/2_Virus.csv",index_col=0)

def label(x):
	return x.split("_")[1]

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
y=nsel[nsel<4].index #In 5% of the patients
df.drop(columns=y,inplace=True)

df=df.div(df.sum(axis=1), axis=0) #Normalize

df["label"]=df.index.map(label)

#Checking the number of patients and its samples
pat_dic={}
def label_pat(x):
	return x.split("_")

pat=df.index.map(label_pat)
for i in pat:
	pat_dic.setdefault(i[0], [])
	pat_dic[i[0]].append(i[1])

print(pat_dic)

#Print the values in various classes
print(df["label"].value_counts())

df.to_csv("Microbes.csv")
