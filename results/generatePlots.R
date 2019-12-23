library(ggplot2)
library(gridExtra)


r=read.table("orderGlobalEval.csv",header=TRUE,sep=",")

nData=1012 #number of data points
foundAbove50Perc=r[r[,"nFound"]>1012/2.0 & r[,"MCC"]>0.7,]
r50=foundAbove50Perc
ordered=r50[order(-r50$nFound,-r50$MCC),]
#plot(ordered[,"similarityThreshold"],ordered[,"nFound"]/nData)
ordered$similarityThreshold <- as.factor(ordered$similarityThreshold)


theme_Publication <- function(base_size=14, base_family="Roboto") {
      library(grid)
      library(ggthemes)
      (theme_foundation(base_size=base_size, base_family=base_family)
       + theme(plot.title = element_text(face = "bold",
                                         size = rel(1.2), hjust = 0.5),
               text = element_text(),
               panel.background = element_rect(colour = NA),
               plot.background = element_rect(colour = NA),
               panel.border = element_rect(colour = NA),
               axis.title = element_text(face = "bold",size = rel(1)),
               axis.title.y = element_text(angle=90,vjust =2),
               axis.title.x = element_text(vjust = -0.2),
               axis.text = element_text(), 
               axis.line = element_line(colour="black"),
               axis.ticks = element_line(),
               panel.grid.major = element_line(colour="#f0f0f0"),
               panel.grid.minor = element_blank(),
               legend.key = element_rect(colour = NA),
               legend.position = "bottom",
               legend.direction = "horizontal",
               legend.key.size= unit(0.2, "cm"),
               legend.margin = unit(0, "cm"),
               legend.title = element_text(face="italic"),
               plot.margin=unit(c(10,5,5,5),"mm"),
               strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
               strip.text = element_text(face="bold")
          ))
      
}

scale_fill_Publication <- function(...){
      library(scales)
      discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)

}

scale_colour_Publication <- function(...){
      library(scales)
      discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)

}


createPlot<-function(data,title){
  p=ggplot(data, aes(similarityThreshold)) + 
  geom_line(aes(y=(MCC*nFound/nData)/(MCC+3*nFound/nData)))+
  geom_line(aes(y=MCC,colour="MCC")) +
  geom_errorbar(aes(ymin=MCC-sdMCC, ymax=MCC+sdMCC,colour="MCC"), width=.01) +
  geom_line(aes(y=F1,colour="F1")) +
  geom_errorbar(aes(ymin=F1-sdF1, ymax=F1+sdF1,colour="F1"), width=.01) +
  geom_line(aes(y=BA,colour="BA")) +
  geom_line(aes(y=nFound/nData,colour="Ratio Found")) + 

  geom_errorbar(aes(ymin=BA-sdBA, ymax=BA+sdBA,colour="BA"), width=.01) +
  geom_errorbar(aes(ymin=nFound/nData-sdnFound/nData, ymax=nFound/nData+sdnFound/nData,colour="Ratio Found"), width=.01) +
  labs(y = "",
                x = "Similarity Threshold",colour="",title=title) + 
                scale_colour_Publication() +
                
                theme_Publication()

  return(p)
  


}



tanimoto_top_radius_1=r[r[,"FingerprintRadius"]==1 & r[,"SimilarityMetric"]=="TanimotoSimilarity" & r[,"evaluationMode"]=="top",]
tanimoto_top_radius_2=r[r[,"FingerprintRadius"]==2 & r[,"SimilarityMetric"]=="TanimotoSimilarity" & r[,"evaluationMode"]=="top",]
tanimoto_top_radius_3=r[r[,"FingerprintRadius"]==3 & r[,"SimilarityMetric"]=="TanimotoSimilarity" & r[,"evaluationMode"]=="top",]

dice_top_radius_1=r[r[,"FingerprintRadius"]==1 & r[,"SimilarityMetric"]=="DiceSimilarity" & r[,"evaluationMode"]=="top",]
dice_top_radius_2=r[r[,"FingerprintRadius"]==2 & r[,"SimilarityMetric"]=="DiceSimilarity" & r[,"evaluationMode"]=="top",]
dice_top_radius_3=r[r[,"FingerprintRadius"]==3 & r[,"SimilarityMetric"]=="DiceSimilarity" & r[,"evaluationMode"]=="top",]

tversky_top_radius_1=r[r[,"FingerprintRadius"]==1 & r[,"SimilarityMetric"]=="TverskySimilarity" & r[,"evaluationMode"]=="top",]
tversky_top_radius_2=r[r[,"FingerprintRadius"]==2 & r[,"SimilarityMetric"]=="TverskySimilarity" & r[,"evaluationMode"]=="top",]
tversky_top_radius_3=r[r[,"FingerprintRadius"]==3 & r[,"SimilarityMetric"]=="TverskySimilarity" & r[,"evaluationMode"]=="top",]

tanimoto_proba_radius_1=r[r[,"FingerprintRadius"]==1 & r[,"SimilarityMetric"]=="TanimotoSimilarity" & r[,"evaluationMode"]=="probability",]
tanimoto_proba_radius_2=r[r[,"FingerprintRadius"]==2 & r[,"SimilarityMetric"]=="TanimotoSimilarity" & r[,"evaluationMode"]=="probability",]
tanimoto_proba_radius_3=r[r[,"FingerprintRadius"]==3 & r[,"SimilarityMetric"]=="TanimotoSimilarity" & r[,"evaluationMode"]=="probability",]

dice_proba_radius_1=r[r[,"FingerprintRadius"]==1 & r[,"SimilarityMetric"]=="DiceSimilarity" & r[,"evaluationMode"]=="probability",]
dice_proba_radius_2=r[r[,"FingerprintRadius"]==2 & r[,"SimilarityMetric"]=="DiceSimilarity" & r[,"evaluationMode"]=="probability",]
dice_proba_radius_3=r[r[,"FingerprintRadius"]==3 & r[,"SimilarityMetric"]=="DiceSimilarity" & r[,"evaluationMode"]=="probability",]

tversky_proba_radius_1=r[r[,"FingerprintRadius"]==1 & r[,"SimilarityMetric"]=="TverskySimilarity" & r[,"evaluationMode"]=="probability",]
tversky_proba_radius_2=r[r[,"FingerprintRadius"]==2 & r[,"SimilarityMetric"]=="TverskySimilarity" & r[,"evaluationMode"]=="probability",]
tversky_proba_radius_3=r[r[,"FingerprintRadius"]==3 & r[,"SimilarityMetric"]=="TverskySimilarity" & r[,"evaluationMode"]=="probability",]

ergfp_top=r[r[,"FingerprintRadius"]==0 & r[,"SimilarityMetric"]=="calc_ergfp" & r[,"evaluationMode"]=="top",]
ergfp_proba=r[r[,"FingerprintRadius"]==0 & r[,"SimilarityMetric"]=="calc_ergfp" & r[,"evaluationMode"]=="probability",]


tt1=createPlot(tanimoto_top_radius_1,"Tanimoto Top MFP1")
tt2=createPlot(tanimoto_top_radius_2,"Tanimoto Top MFP2")
tt3=createPlot(tanimoto_top_radius_3,"Tanimoto Top MFP3")

tp1=createPlot(tanimoto_proba_radius_1,"Tanimoto Prob MFP1")
tp2=createPlot(tanimoto_proba_radius_2,"Tanimoto Prob MFP2")
tp3=createPlot(tanimoto_proba_radius_3,"Tanimoto Prob MFP3")

dt1=createPlot(dice_top_radius_1,"Dice Top MFP1")
dt2=createPlot(dice_top_radius_2,"Dice Top MFP2")
dt3=createPlot(dice_top_radius_3,"Dice Top MFP3")

dp1=createPlot(dice_proba_radius_1,"Dice Prob MFP1")
dp2=createPlot(dice_proba_radius_2,"Dice Prob MFP2")
dp3=createPlot(dice_proba_radius_3,"Dice Prob MFP3")

tvt1=createPlot(tversky_top_radius_1,"Tversky Top MFP1")
tvt2=createPlot(tversky_top_radius_2,"Tversky Top MFP2")
tvt3=createPlot(tversky_top_radius_3,"Tversky Top MFP3")

tvp1=createPlot(tversky_proba_radius_1,"Tversky Prob MFP1")
tvp2=createPlot(tversky_proba_radius_2,"Tversky Prob MFP2")
tvp3=createPlot(tversky_proba_radius_3,"Tversky Prob MFP3")

ergt1=createPlot(ergfp_top,"ERG Top")
ergp1=createPlot(ergfp_proba,"ERG Proba")


topt=grid.arrange(tt1, tt2,tt3,dt1,dt2,dt3,tvt1,tvt2,tvt3 ,nrow = 3)

proba=grid.arrange(tp1,tp2,tp3,dp1,dp2,dp3,tvp1,tvp2,tvp3,nrow=3)

erg=grid.arrange(ergt1,ergp1,nrow=1)

ggsave(file="tanimoto_top_mfp1.svg", plot=topt, width=10, height=10)
ggsave(file="tanimoto_proba_mfp1.svg", plot=proba, width=10, height=10)

ggsave(file="ergfp.svg", plot=erg, width=8, height=4)


# p=ggplot(tanimoto_top_radius_1, aes(similarityThreshold)) + 
#   geom_line(aes(y=MCC,colour="MCC"),color="blue") +
#   geom_line(aes(y=F1,colour="F1"),color="green") +
#   geom_line(aes(y=BA,colour="BA"),color="orange" ) +
#   geom_line(aes(y=nFound/nData,colour="Ratio Found"),color="red" ) + 
#   labs(y = "",
#                 x = "Similarity Threshold",
#                 colour = "Parameter")


#print(p)

#p=ggplot(ordered, aes(x=similarityThreshold, y=nFound)) + 
#  geom_violin() + geom_boxplot(width=0.1)
#print(p)