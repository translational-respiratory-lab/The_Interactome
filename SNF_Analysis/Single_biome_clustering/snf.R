library("SNFtool")
library("vegan")
library("reticulate")
use_python("/usr/local/bin/python")
source_python("./../sil.py")
source("modified_est-cluster.R")

b_data=read.csv("./../../Data/1_Targeted/bacteria.csv",row.names = 1)
f_data=read.csv("./../../Data/1_Targeted/fungi.csv",row.names = 1)
v_data=read.csv("./../../Data/1_Targeted/virus.csv",row.names = 1)

#Filtering the data
list_sel=list()
count=1
for (d in list(b_data,f_data,v_data)){
  z=colSums(d>0)
  sel_col=row.names(as.data.frame(z[z>=10])) #In 5% patients prevalent
  #print(sel_col)
  list_sel[[count]]<-sel_col
  count=count+1
}

b_data<-b_data[,list_sel[[1]]]
f_data<-f_data[,list_sel[[2]]]
v_data<-v_data[,list_sel[[3]]]
remove(d,list_sel,count,sel_col,z)

b_dsim=vegdist(b_data,method='bray',diag=TRUE,upper=TRUE)
f_dsim=vegdist(f_data,method='bray',diag=TRUE,upper=TRUE)
v_dsim=vegdist(v_data,method='bray',diag=TRUE,upper=TRUE)

v_dsim[is.nan(v_dsim)]<-0 #As disimilarity is zero if both patients dont have any virus

W1=(as.matrix(b_dsim)-1)*-1
W2=(as.matrix(f_dsim)-1)*-1
W3=(as.matrix(v_dsim)-1)*-1


names=c("bacteria","fungi","virus")
count<-1
for (x in list(W1,W2,W3)){
  print(names[count])
  y<-estimateNumberOfClustersGivenGraph(x)
  print("First best")
  labels=spectralClustering(x,y$`Eigen-gap best`)
  print(table(labels))
  print("second best")
  labels=spectralClustering(x,y$`Eigen-gap 2nd best`)
  print(table(labels))
  count<-count+1
  }
#From the above output we identify these as the optimum number of clusters
clusters=c(3,2,4)
count<-1
for (x in list(W1,W2,W3)){
labels=spectralClustering(x,clusters[count])
lab=as.data.frame(labels,row.names = row.names(b_data))
print(silhouette_score(x,labels))
write.csv(lab,paste('./results/',names[count],"_labels",".csv",sep=''))
write.csv(x,paste('./results/',names[count],"_matrix",".csv",sep=''))
count=count+1
}