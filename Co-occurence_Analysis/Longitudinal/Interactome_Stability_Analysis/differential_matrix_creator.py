#Importing necessary modules
import pandas
import numpy
import seaborn as sns
import matplotlib.pyplot as plt

#Reading the data
bsl=pandas.read_csv("./../Co-occurance_Analysis/Merge_p-val_scores/BSL_Adj_cyto.csv",index_col=0)
ex1=pandas.read_csv("./../Co-occurance_Analysis/Merge_p-val_scores/EX1_Adj_cyto.csv",index_col=0)
p1=pandas.read_csv("./../Co-occurance_Analysis/Merge_p-val_scores/P1_Adj_cyto.csv",index_col=0)

#Checking if the columns are consistent
print(bsl.index==ex1.index)
print(ex1.index==p1.index)
ind=bsl.index.values #Renamed the columns of bsl manually from code names to their biological names hence we observe some false values

#Combining all the values into a 3D matrix, 23 is the number of microbes
#BEP
comb_mat=numpy.empty((23,23,3))
comb_mat[:,:,0]=bsl.values
comb_mat[:,:,1]=ex1.values
comb_mat[:,:,2]=p1.values

#Uncomment the respective parts run BE, BP or EP

'''
#BE
comb_mat=numpy.empty((23,23,2))
comb_mat[:,:,0]=bsl.values
comb_mat[:,:,1]=ex1.values

#BP
comb_mat=numpy.empty((23,23,2))
comb_mat[:,:,0]=bsl.values
comb_mat[:,:,1]=p1.values

#EP
comb_mat=numpy.empty((23,23,2))
comb_mat[:,:,0]=ex1.values
comb_mat[:,:,1]=p1.values
'''

values=numpy.amax(comb_mat,axis=2)-numpy.amin(comb_mat,axis=2) #Differential score calculation max-min 
values_df=pandas.DataFrame(values,index=bsl.index,columns=bsl.columns)

fig, ax = plt.subplots(figsize=(19,15))
sns.heatmap(values_df,cmap="cubehelix",ax=ax)
plt.ylabel("")
plt.ylim(23,0)
plt.xlim(0,23)
locs, labels = plt.xticks()
plt.xticks(locs,labels,rotation=45,rotation_mode='anchor',ha='right')
#plt.show()
plt.savefig("Figure.png",dpi=300)
