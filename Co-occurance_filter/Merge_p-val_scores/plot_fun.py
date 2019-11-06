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
	Adj.to_csv(str(i)+"_Adj_cyto.csv")

#--------Selecting p-values cutoff-------------
def p_select(Adj,p,p_cut_off,adj_cut_off):
    p.fillna(1,inplace=True)
    ind=p>p_cut_off
    Adj[ind]=0
    Adj[abs(Adj)<adj_cut_off]=0
    return(Adj)
#-------------------------------------------
