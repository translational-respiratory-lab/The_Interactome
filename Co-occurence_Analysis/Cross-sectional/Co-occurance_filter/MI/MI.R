#Reading the dataset
data=read.csv("./../Microbes.csv",row.names = 1)

c1_data=subset(data,x==1)
c2_data=subset(data,x==2)
c1_data$x<-NULL; c2_data$x<-NULL

remove(data)

count=1
files=c("cluster1","cluster2")
for (dd in list(c1_data,c2_data)){
data=dd
#Mutual Information
library("minet")
m<-build.mim(data,estimator = "mi.mm",disc = "equalwidth")
net<-aracne(m)
write.csv(net,paste(files[count],"Adj.csv",sep=''))

#Bootstrap
max_iter=100
rows=rownames(data)
boot_arr<-array(dim=c(dim(net)[1],dim(net)[1],max_iter))
for (i in 1:max_iter){
  t_x=sample(rows,replace = TRUE)
  m=build.mim(data[t_x,],estimator = "mi.mm",disc = "equalwidth")
  net<-aracne(m)
  boot_arr[,,i]<-net
}

#perm and renorm
max_iter=100
perm_arr<-array(dim=c(dim(net)[1],dim(net)[1],max_iter))
for (i in 1:max_iter){
  x<-apply(data,2,FUN =sample) #permutation
  x<-x/rowSums(x) #renormalization
  m=build.mim(x,estimator = "mi.mm",disc = "equalwidth")
  net<-aracne(m)
  print(sum(is.na(net)))
  perm_arr[,,i]<-net
}

#Converting NA to 0
boot_arr[is.na(boot_arr)]<-0
perm_arr[is.na(perm_arr)]<-0

p_val=array(dim=c(dim(net)[1],dim(net)[1]))
for (i in 1:dim(net)[1]){
  for (j in 1:dim(net)[1]){
    t=wilcox.test(perm_arr[i,j,],boot_arr[i,j,],alternative = "two.sided",paired = FALSE,exact = TRUE)
    p_val[i,j]<-t$p.value
  }
}

p<-p.adjust(p_val,method = "fdr")
dim(p)<-dim(p_val)
write.csv(p,paste(files[count],"p_val.csv",sep=''))
count=count+1
}

#All Nan instances and warnings are due to same values ie, all elements are 0 in both
#boot and permutation array
