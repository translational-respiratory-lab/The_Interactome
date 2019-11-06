library("SNFtool")
library("vegan")
source("modified_functions.R")

b_data=read.csv("./../Data/bacteria.csv",row.names = 1)
f_data=read.csv("./../Data/fungi.csv",row.names = 1)
v_data=read.csv("./../Data/virus.csv",row.names = 1)

#Filter and finding the weights
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

library(dunn.test)
X_data=read.csv("./../Data/clinical_attr.csv",row.names = 1)
count<-1
for (nam in names){
  labels=read.csv(paste("./results/",nam,"_labels.csv",sep=""),row.names = 1)
  
  Y_data=merge(X_data,labels,by="row.names",all.y = TRUE)
  row.names(Y_data)<-Y_data$Row.names;Y_data$Row.names<-NULL
  Y_data$x<-factor(Y_data$labels)
  Y_data$labels<-NULL
  
  Y_data$Matching<-NULL
  Y_data$Status<-NULL
  for (i in colnames(Y_data)){print(i);print(class(Y_data[[i]]));print(summary(Y_data[[i]]))}
  Y_data$Gender..Male.0..Female.1.<-factor(Y_data$Gender..Male.0..Female.1.)
  Y_data$Smoking.status..Past.0..Current.1..Never.2.<-factor(Y_data$Smoking.status..Past.0..Current.1..Never.2.)
  Y_data$Ethinicity..CHinese.0..malay.1..indian.2..european.3..vietnamese.4..burmese.5..pakistani.6..<- NULL
  Y_data$Broncoodialator..Yes.0..No.1..<-factor(Y_data$Broncoodialator..Yes.0..No.1..)
  Y_data$Inhaled.corticosteroids..Yes.0..No.1.<-factor(Y_data$Inhaled.corticosteroids..Yes.0..No.1.)
  Y_data$Mucolytic..Yes.0..No.1.<-factor(Y_data$Mucolytic..Yes.0..No.1.)
  Y_data$Long.term.antibiotics..Yes.0..No.1.<-factor(Y_data$Long.term.antibiotics..Yes.0..No.1.)
  Y_data$Antifungal..Yes.0..No.1.<-factor(Y_data$Antifungal..Yes.0..No.1.)
  Y_data$Oral.steroids..Yes.0..No.1.<-factor(Y_data$Oral.steroids..Yes.0..No.1.)
  Y_data$Pneumococcal.vaccine..Yes.0..No.1.<-factor(Y_data$Pneumococcal.vaccine..Yes.0..No.1.)
  Y_data$Influenza.vaccine..Yes.0..No.1.<-factor(Y_data$Influenza.vaccine..Yes.0..No.1.)
  Y_data$Pseudomonas.aeruginosa.status..by.qPCR...positive.0..negative.1.<-factor(Y_data$Pseudomonas.aeruginosa.status..by.qPCR...positive.0..negative.1.)

  dif<-list()
  for (i in 1:80){
    if(class(Y_data[,i])=="integer" | class(Y_data[,i])=="numeric"){
      form=as.formula(paste(names(Y_data)[i],"x",sep="~"))
      print(colnames(Y_data)[i])
      t=wilcox.test(form,data=Y_data)
      #Multiply the p-value value by 2 to account for both side testing
      dif[[colnames(Y_data)[i]]]<-(t$p.value)}
    else{
      # For categorical data comparision
      print(colnames(Y_data)[i])
      f<-fisher.test(Y_data[[i]],Y_data$x,simulate.p.value = TRUE)
      dif[[colnames(Y_data)[i]]]<-(f$p.value)}
  }
  
  #Choosing factors with Alpha value of 0.05
  ch_dif<-dif[dif<0.05]
  sel=names(ch_dif)
  print(sel)
  
  pop1=subset(Y_data,x==1)
  pop2=subset(Y_data,x==2)
  
  sink(paste(nam,"table.csv",sep=''))
  for (i in sel){
    if(class(pop1[[i]])=="integer" | class(pop1[[i]])=="numeric"){
      x=paste(i,median(pop1[[i]],na.rm = TRUE),var(pop1[[i]],na.rm = TRUE),median(pop2[[i]],na.rm = TRUE),var(pop2[[i]],na.rm = TRUE),ch_dif[[i]],sep=",")
      print(x)
    }
    else{
      x=paste(i,toString(table(pop1[[i]])),toString(table(pop2[[i]])),ch_dif[[i]],sep=",")
      print(x)
    }
  }
  sink()
count=count+1
}  
