library("SNFtool")
library("vegan")

b_data=read.csv("./../../Data/1_Targeted/bacteria.csv",row.names = 1)
f_data=read.csv("./../../Data/1_Targeted/fungi.csv",row.names = 1)
v_data=read.csv("./../../Data/1_Targeted/virus.csv",row.names = 1)

#Filter
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

#217 is the number of patients
  for (i in 2:217){
    W = SNF(list(W1,W2,W3),i,20)
    z=estimateNumberOfClustersGivenGraph(W)$`Eigen-gap best`
    print(paste(c(i,z),sep = " "))
    labels=spectralClustering(W,z)
    write.csv(labels,paste(c("Tuning_k/labels",i,".csv"),collapse = ""))
    write.csv(W,paste(c('Tuning_k/matrix',i,".csv"),collapse = ""))
  }

#Run sil_compile.py from the Tuning_k, to compile the silhouette values
#Sort the k_sil.csv to find the optimum k (i.e max silhouette)
#Best k-was chosed to be 9
  W = SNF(list(W1,W2,W3),9,20)
  z=estimateNumberOfClustersGivenGraph(W)$`Eigen-gap best`
  labels=spectralClustering(W,z)
  print(table(labels))
  lab=as.data.frame(labels,row.names = row.names(b_data))
  write.csv(lab,paste("./results/labels",".csv",sep=''))
  write.csv(W,paste("./results/matrix",".csv",sep=''))