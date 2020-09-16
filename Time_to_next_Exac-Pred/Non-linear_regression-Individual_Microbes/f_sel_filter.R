library(compositions)

f_selected=list()


for(i in c("BSL","EX1","P1")){
  
  X_m<-NA
  for(j in c("bac","fungi","vir")){
    nam=paste0(j,"_",i)
    print(nam)
    data<-read.csv(paste0("./../Pre-processing/Pre-processed/",nam,".csv"),row.names = 1)
    t_ex=data$Time.to.next.exacerbation
    y<-ifelse(data$Time.to.next.exacerbation<7*12,"1","2")
    X<-data[,1:(dim(data)[2]-1)]
    X=as.data.frame(clr(X))
    X_m=cbind(X_m,X)
  }

  X_m$X_m<-as.factor(y)

  #Feature_selection

  data<-read.csv(paste0("./../Co-occurence_Analysis/",i,"_C1/Microbes.csv"),row.names = 1)
  sel<-colnames(data)
  f_sel=X_m[c("X_m",sel)]
  f_selected[[i]]<-f_sel  
}


#====Logistic regressoion===========

#Co-relation plot
#===========Correlation Plot================
library("Hmisc")
library(corrplot)

data<-f_selected$BSL

data$X_m<-t_ex
colnames(data)[1]<-"Time to next exacerbation"


res<-rcorr(as.matrix(data),type="spearman")
col1 <- colorRampPalette(c("#3B0EEB","white","#FF3815"))
# png("1A_1.png",width=18, height=18,units="in",res=300)
# corrplot(res$r, method="color",type="upper", title="", p.mat = res$P, sig.level = 0.05,
#          insig = "blank",hclust.method = "ward.D2",tl.cex=1,tl.col="black",tl.srt = 90,col=col1(200),addgrid.col = "black",
#          order="original",cl.cex = 1,na.label="square",na.label.col = "white")

ind<-(res$P[1,]<0.05)
ind[1]<-TRUE #first
y<-as.data.frame(res$r[1,])
y[(ind == FALSE),]<-0

write.csv(y,"BSL.csv")

library(earth)

model<-earth(X_m ~ .,data = f_selected$P1,degree = 10)
sink("P1.txt")
summary(model)
sink()

data=evimp(model)

data = as.data.frame(unclass(data[,c(3,4,6)]))
data$var<-row.names(data)

png("P1_fimp.png",width=6, height=3,units="in",res=300)
library(ggplot2)
# Basic barplot
p<-ggplot(data=data, mapping = aes(x = reorder(var, gcv), gcv)) +
  geom_bar(stat="identity",width = 0.3)+ylab("Importance(gcv)")+xlab("")
# Horizontal bar plot
p<-p + coord_flip()+scale_color_grey() +theme(axis.text.y=element_text(face="italic"),panel.background = element_blank())
print(p)
dev.off()