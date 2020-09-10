library("SNFtool")
library("vegan")
source("function_snf.R")
library("reticulate")
use_python("/usr/local/bin/python")
source_python("sil.py")

b_data=read.csv("./../../Data/1_Targeted/bacteria.csv",row.names = 1)
f_data=read.csv("./../../Data/1_Targeted/fungi.csv",row.names = 1)
v_data=read.csv("./../../Data/1_Targeted/virus.csv",row.names = 1)

#Filterting the dataset
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

#Assigning weight values
weight_b=dim(b_data)[2]
weight_f=dim(f_data)[2]
weight_v=dim(v_data)[2]

weights_snf=c(weight_b,weight_f,weight_v)
sil_values=c()
for (i in 2:217){
  W = SNF_weighted_iter(list(W1,W2,W3),i,20,weight = weights_snf)
  z=estimateNumberOfClustersGivenGraph(W)$`Eigen-gap best`
  labels=spectralClustering(W,z)
  sil_values<-c(sil_values,silhouette_score(W,labels))
}
tuned_k<-which.max(sil_values)+1 #since starts from 2
print(paste(tuned_k,sil_values[tuned_k-1],sep = " "))

W = SNF_weighted_iter(list(W1,W2,W3),tuned_k,20,weight = weights_snf)
z=estimateNumberOfClustersGivenGraph(W)$`Eigen-gap best`
labels=spectralClustering(W,z)
print(table(labels))
lab=as.data.frame(labels,row.names = row.names(b_data))
write.csv(lab,paste("./results/labels",".csv",sep=''))
write.csv(W,paste("./results/matrix",".csv",sep=''))