library("SNFtool")
library("vegan")
source("function_snf.R")
library("reticulate")
use_python("/usr/local/bin/python3")
source_python("sil.py")

b_data=read.csv("./Pre_processed_data/1_Bacteria_mod.csv",row.names = 1)
f_data=read.csv("./Pre_processed_data/3_Fungi_mod.csv",row.names = 1)
v_data=read.csv("./Pre_processed_data/2_Virus_mod.csv",row.names = 1)
pat_label=read.csv("./Pre_processed_data/patient_labels_color.csv")

groups_=b_data$label
groups_p=as.factor(pat_label$groups)
b_data$label<-NULL
f_data$label<-NULL
v_data$label<-NULL

b_dsim=vegdist(b_data,method='bray',diag=TRUE,upper=TRUE)
f_dsim=vegdist(f_data,method='bray',diag=TRUE,upper=TRUE)
v_dsim=vegdist(v_data,method='bray',diag=TRUE,upper=TRUE)

v_dsim[is.nan(v_dsim)]<-0 #As disimilarity is zero if both patients dont have any virus

W1=(as.matrix(b_dsim)-1)*-1
W2=(as.matrix(f_dsim)-1)*-1
W3=(as.matrix(v_dsim)-1)*-1

weight_b=dim(b_data)[2]
weight_f=dim(f_data)[2]
weight_v=dim(v_data)[2]

weights_snf=c(weight_b,weight_f,weight_v)
sil_values=c()
for (i in 2:51){
  W = SNF_weighted_iter(list(W1,W2,W3),i,20,weight = weights_snf)
  z=estimateNumberOfClustersGivenGraph(W)$`Eigen-gap best`
  labels=spectralClustering(W,z)
  sil_values<-c(sil_values,silhouette_score(W,labels))
}
tuned_k<-which.max(sil_values)+1 #since starts from 2
print(paste(tuned_k,sil_values[tuned_k-1],sep = " "))

W = SNF_weighted_iter(list(W1,W2,W3),tuned_k,20,weight = weights_snf)

#Converting Similarity to Distance
W=(0.5)-W

#MDS using vegan
ord<-metaMDS(W,k=2,trymax = 500)
ordiplot(ord)
#BSL,EX1,P1
ordihull(ord,groups = groups_,draw="polygon",col=c("blue","red","green"))