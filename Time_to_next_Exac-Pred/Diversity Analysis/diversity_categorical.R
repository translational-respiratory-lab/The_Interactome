library("vegan")
library(splitstackshape)
library(ggplot2)

for (i in c("BSL","EX1","P1")){
  
  bac<-read.csv(paste0("./../Pre-processing/Pre-processed/bac_",i,".csv"),row.names = 1)
  fun<-read.csv(paste0("./../Pre-processing/Pre-processed/fungi_",i,".csv"),row.names = 1)
  vir<-read.csv(paste0("./../Pre-processing/Pre-processed/vir_",i,".csv"),row.names=1)
  
  lab<-bac$Time.to.next.exacerbation
  bac$Time.to.next.exacerbation<-NULL
  fun$Time.to.next.exacerbation<-NULL
  vir$Time.to.next.exacerbation<-NULL
  
  bac<-(bac/rowSums(bac))*100
  fun<-(fun/rowSums(fun))*100
  vir<-(vir/rowSums(vir))*100; vir[is.na(vir)]<-0
  
  vir$Time.to.next.exacerbation<-lab
  
  #Merged data
  tmp_data<-merge(bac,fun,by="row.names",all=TRUE)
  row.names(tmp_data)<-tmp_data$Row.names;tmp_data$Row.names<-NULL
  data<-merge(tmp_data,vir,by="row.names",all=TRUE)
  row.names(data)<-data$Row.names;data$Row.names<-NULL
  
  X<-data;X$Time.to.next.exacerbation<-NULL
  X<-(X/rowSums(X))*100
  
  H<-diversity(X,"shannon") #Shannon
  # H<-diversity(as.matrix(X),"berger-parker",method = "bray-curtis")$berger.parker.D
  

plt_data=as.data.frame(cbind(H,data$Time.to.next.exacerbation))
plt_data$V2<-ifelse(plt_data$V2>12*7,">12","<12")
plt_data$V2=as.factor(plt_data$V2)

png(paste0(i,".png"),height =3.5 ,width = 3,units = "in",res = 300)
p <- ggplot(plt_data,aes(x=V2, y=H,fill=V2)) + 
  scale_fill_manual(values=c("#F8766D","#619CFF"))+
  geom_boxplot(outlier.shape = NA) +geom_jitter(width = 0.03, height = 0)+
  xlab("Time to next exacerbation (weeks)")+
  ylab("Shannon Diversity index")+ggtitle(i)+ylim(1,2.5)+
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"),legend.title = element_blank(),legend.position = "none")
print(p)
dev.off()
print(i)
print(wilcox.test(H~V2,plt_data))
}

#========Not different=============