library("SNFtool")
library("vegan")
library("reticulate")
use_python("/usr/local/bin/python")
source_python("./../sil.py")
source("modified_est-cluster.R")

b_data=read.csv("./../../Data/1_Targeted/bacteria.csv",row.names = 1)
f_data=read.csv("./../../Data/1_Targeted/fungi.csv",row.names = 1)
v_data=read.csv("./../../Data/1_Targeted/virus.csv",row.names = 1)

#Filtering the dataset
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

#Tuning for optimal value of k
names=c("b+f","f+v","b+v")
count<-1
for (x in list(list(W1,W2),list(W2,W3),list(W1,W3))){
  print(names[count])
  for (i in 1:217){
    W = SNF(x,i,20)
    z=estimateNumberOfClustersGivenGraph(W)$`Eigen-gap best`
    print(paste(c(i,z),sep = " "))
    labels=spectralClustering(W,z)
    write.csv(labels,paste(c('Tuning_k/',names[count],"/labels",i,".csv"),collapse = ""))
    write.csv(W,paste(c("Tuning_k/",names[count],'/matrix',i,".csv"),collapse = ""))
  }
  count<-count+1
}

#Calculation of the silhouette values for each k, is done by a python script
#Please run python sil-values_compiler.py from the Tuning_k directory

#We sort the output csv files(of sil-values_compiler.py) to find the optimal k with max silhouette
#(b+f,v+f,b+v)
k<-c(14,5,4)

count<-1
for (x in list(list(W1,W2),list(W2,W3),list(W1,W3))){
  print(names[count])
  W = SNF(x,k[count],20)
  z=estimateNumberOfClustersGivenGraph(W)$`Eigen-gap best`
  labels=spectralClustering(W,z)
  print(table(labels))
  lab=as.data.frame(labels,row.names = row.names(b_data))
  write.csv(lab,paste('./results/',names[count],"_labels",".csv",sep=''))
  write.csv(W,paste('./results/',names[count],"_matrix",".csv",sep=''))
  count=count+1
}