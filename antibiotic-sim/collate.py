from __future__ import division
import pandas

bac=pandas.read_csv("./../Co-occurence_Analysis/Longitudinal/Longitudinal_Data/1_Bacteria.csv",index_col=0)
fun=pandas.read_csv("./../Co-occurence_Analysis/Longitudinal//Longitudinal_Data/3_Fungi.csv",index_col=0)
vir=pandas.read_csv("./../Co-occurence_Analysis/Longitudinal//Longitudinal_Data/2_Virus.csv",index_col=0)

vir=vir.div(vir.sum(axis=1), axis=0)
vir=vir.multiply(100, fill_value=0) #Normalize

bac=bac.div(bac.sum(axis=1), axis=0)
bac=bac.multiply(100, fill_value=0) #Normalize

fun=fun.div(fun.sum(axis=1), axis=0)
fun=fun.multiply(100, fill_value=0) #Normalize

x=bac.join(fun,how="inner")
y=x.join(vir,how="inner")

df=y.div(y.sum(axis=1), axis=0) #Renormalize


#Print the values in various classes
#print(df["label"].value_counts())

# Prediction samples
sam=["14_BSL","28_BSL","125_BSL","169_BSL","205_BSL","217_BSL","276_BSL","334_BSL","394_BSL","425_BSL","515_BSL","540_BSL"]

sam_po=["14_P1","28_P1","125_P1","169_P1","205_P1","217_P1","276_P1","334_P1","394_P1","425_P1","515_P1","540_P1"]

#Microbes
col_red=["Streptococcus","Staphylococcus","Haemophilus","Moraxella","Actinomyces","Arachnia", "Bacteroides", "Bifidobacterium", "Eubacterium","Fusobacterium","Lactobacillus", "Leptotrichia", "Peptococcus", "Peptostreptococcus","Propionibacterium", "Selenomonas", "Treponema", "Veillonella"]

df_copy=df.copy()
#====================For Post-Antibiotic=================================================
df=df.loc[sam_po,]

#Dropping 0 columns
nsel=df.sum(axis=0)
n=nsel[nsel==0].index
df.drop(columns=n,inplace=True)

#Filter
nsel=df>0.01   #atleast 1%
nsel=nsel.sum(axis=0)
y=nsel[nsel<3].index 

y=y.drop("Neisseria") #Due to the OR condition of the filter, need to add neisseria since its selected by pre-antibiotic

df.drop(columns=y,inplace=True)

df=df.div(df.sum(axis=1), axis=0) #Normalize
print(df.shape)
df.to_csv("Microbes_post.csv")
#==================For Pre-Antibiotic================================================
df=df_copy.loc[sam,]

#Dropping 0 columns
nsel=df.sum(axis=0)
n=nsel[nsel==0].index
df.drop(columns=n,inplace=True)

#Filter
nsel=df>0.01   #atleast 1%
nsel=nsel.sum(axis=0)
y=nsel[nsel<3].index 

y=y.drop("Prevotella") #Due to the OR condition of the filter, need to add prevotella since its selected by post-antibiotic

df.drop(columns=y,inplace=True)

df=df.div(df.sum(axis=1), axis=0) #Normalize
print(df.shape)
df.to_csv("Microbes_pre.csv")
#===================For 75% reduction =======================
df=df_copy.loc[sam,]
inter_col=set(df.columns).intersection(col_red)
df.loc[:,inter_col]=0.25*(df.loc[:,inter_col]) 

#Dropping 0 columns
nsel=df.sum(axis=0)
n=nsel[nsel==0].index
df.drop(columns=n,inplace=True)

#Filter
nsel=df>0.01   #atleast 1%
nsel=nsel.sum(axis=0)
y=nsel[nsel<3].index 

y=y.drop("Prevotella") #Due to the OR condition of the filter, need to add prevotella since its selected by post-antibiotic

df.drop(columns=y,inplace=True)

df=df.div(df.sum(axis=1), axis=0) #Normalize
print(df.shape)
df.to_csv("Microbes_75%_reduction.csv")
