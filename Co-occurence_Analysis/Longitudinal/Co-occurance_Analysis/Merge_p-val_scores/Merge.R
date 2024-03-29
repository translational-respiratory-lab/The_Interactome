scoring<-function(arr){
  m<-max(abs(arr),na.rm = TRUE)
  arr<-(arr/m)*100
  return(arr)
}

m_score=0
p_m=array(0,c(23,23,5)) #Change based on the number of microbes present

cluster="P1" #Change the cluster variable and re-run to create output files for different groups
methods=c("Bray-Curtis","GBLM","MI","Pearsons","Spearman")
count=1
w=c(1/6,1/2,1/6,1/12,1/12)
for (i in methods){
  data<-read.csv(paste("./../",i,"/",cluster,"_Adj.csv",sep=""))
  p<-read.csv(paste("./../",i,"/",cluster,"_p_val.csv",sep=""))
  row.names(data)<-data$X
  row.names(p)<-data$X
  p$X<-NULL
  colnames(p)<-data$X
  data$X<-NULL
  p_m[,,count]<-as.matrix(p)
  #ind<-p>0.01   #Unnecessary if combined p-values is used
  #data[ind]<-0
  m_score<-m_score+w[count]*abs(scoring(data))
  count=count+1
}
# Calculating the sign of the score
m_score_sign=0
methods=c("GBLM","Pearsons","Spearman")
count=1
w=c(1/2,1/4,1/4)
for (i in methods){
  data<-read.csv(paste("./../",i,"/",cluster,"_Adj.csv",sep=""))
  row.names(data)<-data$X
  data$X<-NULL
  m_score_sign<-m_score_sign+w[count]*scoring(data)
  count=count+1
}
m_score=m_score*sign(m_score_sign)

write.csv(m_score,paste(cluster,"merged_scores.csv",sep='-'))

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

#library(mppa)
#p_f<-apply(p_m,c(1,2),simes.test)

write.csv(p_f,paste(cluster,"merged_pvalues.csv",sep="-"))
#Please rename the output files by suffixing "BSL_", "EX1_",.. as appropriate based on the cluster variable