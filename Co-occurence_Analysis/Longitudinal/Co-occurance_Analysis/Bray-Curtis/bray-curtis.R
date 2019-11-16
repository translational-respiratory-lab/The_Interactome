#Reading the dataset
data=read.csv("./../Microbes.csv",row.names = 1)

c1_data=subset(data,label=="BSL")
c2_data=subset(data,label=="EX1")
c3_data=subset(data,label=="P1")
c1_data$label<-NULL
c2_data$label<-NULL
c3_data$label<-NULL
data$label<-NULL
#Checking for relative abundance structure
r_s=c()
for (i in 1:18){
  r_s=c(r_s,sum(data[i,]))
}
print(var(r_s))
remove(r_s,i,data)

#Bray-Curtis
#transpose data
dd=c3_data #Change the cluster variable and re-run to create output files for different groups
data=t(dd)
library(vegan)
dis=vegdist(data,method='bray',diag=TRUE,upper=TRUE)
sim=(as.matrix(dis)-1)*-1  #Is the adjacency matrix
write.csv(sim,"Adj.csv")

#Bootstrap
max_iter=100
cols=colnames(data)
boot_arr<-array(dim=c(dim(sim)[1],dim(sim)[1],max_iter))
for (i in 1:max_iter){
  t_x=sample(cols,replace = TRUE)
  dis=vegdist(data[,t_x],method='bray',diag=TRUE,upper=TRUE)
  sim=(as.matrix(dis)-1)*-1
  boot_arr[,,i]<-sim
}

#perm and renorm
max_iter=100
perm_arr<-array(dim=c(dim(sim)[1],dim(sim)[1],max_iter))
for (i in 1:max_iter){
  x<-t(apply(data,1,FUN =sample)) #permutation
  x<-t(t(x)/colSums(x)) #renormalization
  dis=vegdist(x,method='bray',diag=TRUE,upper=TRUE)
  sim=(as.matrix(dis)-1)*-1
  perm_arr[,,i]<-sim
}

#Converting NA to 0
boot_arr[is.na(boot_arr)]<-0
perm_arr[is.na(perm_arr)]<-0

p_val=array(dim=c(dim(sim)[1],dim(sim)[1]))
for (i in 1:dim(sim)[1]){
  for (j in 1:dim(sim)[1]){
    t=wilcox.test(perm_arr[i,j,],boot_arr[i,j,],alternative = "two.sided",paired = FALSE,exact = TRUE)
    p_val[i,j]<-t$p.value
  }
}

p<-p.adjust(p_val,method = "fdr")
dim(p)<-dim(p_val)
write.csv(p,"p_val.csv")
#Please rename the output files by suffixing "BSL_", "EX1_",.. as appropriate based on the cluster variable