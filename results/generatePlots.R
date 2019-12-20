library(ggplot2)

r=read.table("orderGlobalEval.csv",header=TRUE,sep=",")

nData=1012 #number of data points
foundAbove50Perc=r[r[,"nFound"]>1012/2.0 & r[,"MCC"]>0.7,]
r50=foundAbove50Perc
ordered=r50[order(-r50$nFound,-r50$MCC),]
#plot(ordered[,"similarityThreshold"],ordered[,"nFound"]/nData)
ordered$similarityThreshold <- as.factor(ordered$similarityThreshold)


createPlot<-function(data){
  p=ggplot(data, aes(similarityThreshold)) + 
  geom_line(aes(y=MCC,colour="MCC"),size=2) +
  geom_line(aes(y=F1,colour="F1"),size=2) +
  geom_line(aes(y=BA,colour="BA"),size=2) +
  geom_line(aes(y=nFound/nData,colour="Ratio Found"),size=2) + 
  geom_errorbar(aes(ymin=MCC-sdMCC, ymax=MCC+sdMCC), width=.01,color="grey") +
  geom_errorbar(aes(ymin=F1-sdF1, ymax=F1+sdF1), width=.01,color="grey") +
  geom_errorbar(aes(ymin=BA-sdBA, ymax=BA+sdBA), width=.01,color="grey") +
  geom_errorbar(aes(ymin=nFound/nData-sdnFound/nData, ymax=nFound/nData+sdnFound/nData), width=.01,,color="grey") +
  labs(y = "",
                x = "Similarity Threshold",
                colour = "Parameter") + 
               scale_color_brewer(palette="Dark2")+theme_minimal()+
                theme(legend.position="bottom")

  return(p)


}



tanimoto_top_radius_1=r[r[,"FingerprintRadius"]==1 & r[,"SimilarityMetric"]=="TanimotoSimilarity" & r[,"evaluationMode"]=="top",]
p=createPlot(tanimoto_top_radius_1)
# p=ggplot(tanimoto_top_radius_1, aes(similarityThreshold)) + 
#   geom_line(aes(y=MCC,colour="MCC"),color="blue") +
#   geom_line(aes(y=F1,colour="F1"),color="green") +
#   geom_line(aes(y=BA,colour="BA"),color="orange" ) +
#   geom_line(aes(y=nFound/nData,colour="Ratio Found"),color="red" ) + 
#   labs(y = "",
#                 x = "Similarity Threshold",
#                 colour = "Parameter")
print(p)

#p=ggplot(ordered, aes(x=similarityThreshold, y=nFound)) + 
#  geom_violin() + geom_boxplot(width=0.1)
#print(p)