library(lionessR)
library(SummarizedExperiment)
data1=read.csv("./../Co-occurence_Analysis/P1_C1/Microbes.csv",row.names = 1)
data2=read.csv("./../Co-occurence_Analysis/P1_C2/Microbes.csv",row.names = 1)
data=rbind(data1,data2)

source("GBLM.R")


data<-as.matrix(data)
data<-t(data)

set.seed(123)
res<-lioness(data,GBLM)

print(dim(assay(res)))

res<-as.data.frame(assay(res))
#Can convert to edge-list format, after splitting on _

res<-t(res)

#filtering edges with no-variance
ind<-apply(res,2,var)>0
res<-res[,ind]

write.csv(res,"edge_weights-across_patients.csv")


y<-read.csv("./../Pre-processing/Pre-processed/bac_BSL.csv",row.names = 1)
y=as.data.frame(y$Time.to.next.exacerbation,row.names = row.names(y))
res<-merge(res,y,by="row.names",all=TRUE);row.names(res)<-res$Row.name;res$Row.names<-NULL


#==========Correlation plot=============
library("Hmisc")
library(corrplot)

data=res

res<-rcorr(as.matrix(data),type="spearman")
# col1 <- colorRampPalette(c("#3B0EEB","white","#FF3815"))
# diag(res$r)<-NA
# png("P1_corr.png",width=18, height=18,units="in",res=300)
# corrplot(res$r, method="color",type="upper", title="", p.mat = res$P, sig.level = 0.05,
#          insig = "blank",hclust.method = "ward.D2",tl.cex=1,tl.col="black",tl.srt = 90,col=col1(200),addgrid.col = "black",
#          order="original",cl.cex = 1,na.label="square",na.label.col = "white")
# dev.off()

ind<-(res$P[1,]<0.05)
ind[1]<-TRUE #first
y<-as.data.frame(res$r[1,])
y[(ind == FALSE),]<-0

write.csv(y,"P1.csv")

library(earth)

model<-earth(`y$Time.to.next.exacerbation` ~ .,data = data,degree = 10)
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