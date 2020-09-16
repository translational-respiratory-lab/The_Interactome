import pandas
import numpy
import matplotlib.pyplot as plt

df=pandas.read_csv("Copy of Time_to_exacerbation.csv",index_col=0)
df=df["Time to next exacerbation "]
'''
df.hist(bins=100)
plt.show()
'''

#Removing 9904 prefix
df.index=[int(str(i).split("9904")[1]) for i in df.index]

bac=pandas.read_csv("2_Virus.csv",index_col=0)

def label(x):
	return x.split("_")[1]

bac["label"]=bac.index.map(label)

for i in ["BSL","EX1","P1"]:
    f_bac=bac[bac.label==i]
    f_bac=f_bac.drop(["label"],axis=1)
    f_bac.index=[int(str(i).split("_")[0]) for i in f_bac.index]
    f_bac=pandas.merge(f_bac,df,how="left",left_index=True,right_index=True)		
    f_bac.to_csv("Pre-processed/"+"vir_"+i+".csv")
    
'''    
bac.index=[int(str(i).split("_")[0]) for i in bac.index]
bac=pandas.merge(bac,df,left_index=True,right_index=True)
bac.to_csv("./Pre-processed/Categorical/"+"vir"+".csv")
'''
