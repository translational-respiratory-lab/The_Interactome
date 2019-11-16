scoring<-function(arr){
  m<-max(abs(arr),na.rm = TRUE)
  arr<-(arr/m)*100
  return(arr)
}

m_score=0
p_m=array(0,c(501,501,5)) #Change based on the number of microbes present

cluster="cluster2"

count=1
methods=c("Bray-Curtis","GBLM","MI","Pearsons","Spearman")
w=c(1/6,1/2,1/6,1/12,1/12) #weights
for (i in methods){
  data<-read.csv(paste("./../",i,"/",cluster,"Adj.csv",sep=""))
  p<-read.csv(paste("./../",i,"/",cluster,"p_val.csv",sep=""))
  row.names(data)<-data$X
  row.names(p)<-data$X
  p$X<-NULL
  colnames(p)<-data$X
  data$X<-NULL
  p_m[,,count]<-as.matrix(p)
  m_score<-m_score+w[count]*abs(scoring(data))
  count=count+1
}
# Calculating the sign of the score
m_score_sign=0
methods=c("GBLM","Pearsons","Spearman")
count=1
w=c(1/2,1/4,1/4)
for (i in methods){
  data<-read.csv(paste("./../",i,"/",cluster,"Adj.csv",sep=""))
  row.names(data)<-data$X
  data$X<-NULL
  m_score_sign<-m_score_sign+w[count]*scoring(data)
  count=count+1
}

m_score=m_score*sign(m_score_sign)

write.csv(m_score,"Merged_scores.csv")

#------------Merging p-values------------------
row.names(p_m)<-row.names(p)
colnames(p_m)<-colnames(p)

sime<-function(p){
  library(gMCP)
  x<-simes.test(p,weights=c(1/6,1/2,1/6,1/12,1/12))
  return(x)
}
p_f<-apply(p_m,c(1,2),sime)
p_f[p_f==0]<-NA
#p-values with NA become 0 because if NA is present it is assigned 0

write.csv(p_f,"merged_pvalues.csv")