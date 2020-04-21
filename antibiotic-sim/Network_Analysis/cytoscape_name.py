# To add the name of edges to the diffany output. Execute this script with the network open in Cytoscape and change the data-fram below accrodingly so as to add the right edge names.
import os
import sys
from time import sleep
import pandas 
from py2cytoscape import cyrest
cytoscape=cyrest.cyclient()
cytoscape.version()
#Shows all network names
print(cytoscape.network.get())
#cytoscape.layout.force_directed(isDeterministic=True,network='Ad_cyto_cutoff.csv',defaultEdgeWeight=0.4)

def edge_name(name):
    name="suid:"+str(name)
    print(name)
    sel=cytoscape.network.select(edgeList=str(name),extendEdges=True)
    node_nam=["suid:"+str(i) for i in sel["nodes"]]
    print(node_nam)
    n1=cytoscape.node.get_attribute(columnList="name",nodeList=node_nam[0])[0]
    n2=cytoscape.node.get_attribute(columnList="name",nodeList=node_nam[1])[0]
    edge_nam=str(n1["name"]+"_"+str(n2["name"]))
    print(edge_nam)
    cytoscape.edge.rename(edge=name, newName=edge_nam)
    cytoscape.network.deselect(edgeList=str(name))
    cytoscape.network.deselect(nodeList=node_nam[0])
    cytoscape.network.deselect(nodeList=node_nam[1])

#df=pandas.read_csv("diffany_75%-reduction_edge weights.csv",index_col=0)
df=pandas.read_csv("diffany_post-antibiotic_edge_weights.csv",index_col=0)
d=df.index.to_list()
for i in d:
    edge_name(i)
