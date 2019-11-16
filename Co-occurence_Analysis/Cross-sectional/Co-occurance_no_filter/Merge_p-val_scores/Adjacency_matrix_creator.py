import pandas
import math
for i in ["cluster1","cluster2"]:
	Adj=pandas.read_csv(str(i)+"_scores.csv")
	p=pandas.read_csv(str(i)+"_pvalues.csv")
	Adj.set_index("Unnamed: 0",inplace=True)
	p.set_index("Unnamed: 0",inplace=True)
	p.fillna(1,inplace=True) #To not select them
	p_cut_off=0.001  #Cut-off of p-values
	ind=p.values>p_cut_off
	Adj.values[ind]=0
	Adj.fillna(0,inplace=True)
	Adj.to_csv(str(i)+"_Adj_cyto.csv") #The output adjacency matrix was imported into Cytoscape for further analysis
