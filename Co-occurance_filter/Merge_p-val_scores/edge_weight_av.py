import pandas

cl1=pandas.read_csv("cluster1_Adj_cyto.csv",index_col=0)
cl2=pandas.read_csv("cluster2_Adj_cyto.csv",index_col=0)

cl1=cl1.abs()
cl2=cl2.abs()

cl1=cl1.sum(axis=1)
cl2=cl2.sum(axis=1)

print "Cluster1 \n"
print cl1.sort_values()
print '\n'
print "Cluster2 \n"
print cl2.sort_values()
