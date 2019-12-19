library(ggplot2)

r=read.table("orderGlobalEval.csv",header=TRUE,sep=",")

nData=1012 #number of data points
foundAbove50Perc=r[r[,"nFound"]>1012/2.0 & r[,"MCC"]>0.7,]
r50=foundAbove50Perc
ordered=r50[order(-r50$nFound,-r50$MCC),]
plot(ordered[,"similarityThreshold"],ordered[,"nFound"]/nData)
ordered$similarityThreshold <- as.factor(ordered$similarityThreshold)

ggplot(ordered, aes(x=similarityThreshold, y=nFound)) + 
  geom_violin() + geom_boxplot(width=0.1)