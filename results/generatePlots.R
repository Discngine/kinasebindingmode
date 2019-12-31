library(ggplot2)
library(gridExtra)
library(reshape2)

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
                                         size = rel(5.2), hjust = 0.5),
               text = element_text(size = rel(5.2)),
               panel.background = element_rect(colour = NA),
               plot.background = element_rect(colour = NA),
               panel.border = element_rect(colour = NA),
               axis.title = element_text(face = "bold",size = rel(1.0)),
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
               legend.key.size= unit(0.5, "cm"),
               legend.margin = unit(0, "cm"),
               legend.title = element_text(face="italic"),
               legend.text = element_text(size = rel(5.2)),
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
                theme_Publication() + 
                ylim(0,1.0)

  
  

  return(p)
  


}

getPerformance<-function(data,threshold){
  data=data[order(-data["MCC"]),]
  idx=which(data[,"nFound"]/nData>=threshold)
  print(idx)
  idx=idx[1]
  return(c(data[idx,"BA"],data[idx,"sdBA"],data[idx,"F1"],data[idx,"sdF1"],data[idx,"MCC"],data[idx,"sdMCC"],data[idx,"similarityThreshold"]))
}


catPerformance<-function(perf,title,file){
  cat(title,"\t",perf[1],"\t",perf[2],"\t",perf[3],"\t",perf[4],"\t",perf[5],"\t",perf[6],"\t",perf[7],"\n",file=file,append=TRUE)

  
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

cat("Combination\tBA\tsdBA\tF1\tsdF1\tMCC\tsdMCC\tsimilarityThreshold",file="top_results_05.csv",sep="\n")
cat("Combination\tBA\tsdBA\tF1\tsdF1\tMCC\tsdMCC\tsimilarityThreshold",file="top_results_08.csv",sep="\n")
cat("Combination\tBA\tsdBA\tF1\tsdF1\tMCC\tsdMCC\tsimilarityThreshold",file="top_results_09.csv",sep="\n")
cat("Combination\tBA\tsdBA\tF1\tsdF1\tMCC\tsdMCC\tsimilarityThreshold",file="top_results_1.csv",sep="\n")
catPerformance(getPerformance(tanimoto_top_radius_1,0.8),"Tanimoto_top_1_0.8%","top_results_08.csv")
catPerformance(getPerformance(tanimoto_top_radius_2,0.8),"Tanimoto_top_2_0.8%","top_results_08.csv")
catPerformance(getPerformance(tanimoto_top_radius_3,0.8),"Tanimoto_top_3_0.8%","top_results_08.csv")
catPerformance(getPerformance(dice_top_radius_1,0.8),"Dice_top_1_0.8%","top_results_08.csv")
catPerformance(getPerformance(dice_top_radius_2,0.8),"Dice_top_2_0.8%","top_results_08.csv")
catPerformance(getPerformance(dice_top_radius_3,0.8),"Dice_top_3_0.8%","top_results_08.csv")
catPerformance(getPerformance(tversky_top_radius_1,0.8),"Tversky_top_1_0.8%","top_results_08.csv")
catPerformance(getPerformance(tversky_top_radius_2,0.8),"Tversky_top_2_0.8%","top_results_08.csv")
catPerformance(getPerformance(tversky_top_radius_3,0.8),"Tversky_top_3_0.8%","top_results_08.csv")

catPerformance(getPerformance(tanimoto_top_radius_1,0.9),"Tanimoto_top_1_0.9%","top_results_09.csv")
catPerformance(getPerformance(tanimoto_top_radius_2,0.9),"Tanimoto_top_2_0.9%","top_results_09.csv")
catPerformance(getPerformance(tanimoto_top_radius_3,0.9),"Tanimoto_top_3_0.9%","top_results_09.csv")
catPerformance(getPerformance(dice_top_radius_1,0.9),"Dice_top_1_0.9%","top_results_09.csv")
catPerformance(getPerformance(dice_top_radius_2,0.9),"Dice_top_2_0.9%","top_results_09.csv")
catPerformance(getPerformance(dice_top_radius_3,0.9),"Dice_top_3_0.9%","top_results_09.csv")
catPerformance(getPerformance(tversky_top_radius_1,0.9),"Tversky_top_1_0.9%","top_results_09.csv")
catPerformance(getPerformance(tversky_top_radius_2,0.9),"Tversky_top_2_0.9%","top_results_09.csv")
catPerformance(getPerformance(tversky_top_radius_3,0.9),"Tversky_top_3_0.9%","top_results_09.csv")

catPerformance(getPerformance(tanimoto_top_radius_1,0.5),"Tanimoto_top_1_0.5%","top_results_05.csv")
catPerformance(getPerformance(tanimoto_top_radius_2,0.5),"Tanimoto_top_2_0.5%","top_results_05.csv")
catPerformance(getPerformance(tanimoto_top_radius_3,0.5),"Tanimoto_top_3_0.5%","top_results_05.csv")
catPerformance(getPerformance(dice_top_radius_1,0.5),"Dice_top_1_0.5%","top_results_05.csv")
catPerformance(getPerformance(dice_top_radius_2,0.5),"Dice_top_2_0.5%","top_results_05.csv")
catPerformance(getPerformance(dice_top_radius_3,0.5),"Dice_top_3_0.5%","top_results_05.csv")
catPerformance(getPerformance(tversky_top_radius_1,0.5),"Tversky_top_1_0.5%","top_results_05.csv")
catPerformance(getPerformance(tversky_top_radius_2,0.5),"Tversky_top_2_0.5%","top_results_05.csv")
catPerformance(getPerformance(tversky_top_radius_3,0.5),"Tversky_top_3_0.5%","top_results_05.csv")

catPerformance(getPerformance(tanimoto_top_radius_1,0.99),"Tanimoto_top_1_0.99%","top_results_099.csv")
catPerformance(getPerformance(tanimoto_top_radius_2,0.99),"Tanimoto_top_2_0.99%","top_results_099.csv")
catPerformance(getPerformance(tanimoto_top_radius_3,0.99),"Tanimoto_top_3_0.99%","top_results_099.csv")
catPerformance(getPerformance(dice_top_radius_1,0.99),"Dice_top_1_0.99%","top_results_099.csv")
catPerformance(getPerformance(dice_top_radius_2,0.99),"Dice_top_2_0.99%","top_results_099.csv")
catPerformance(getPerformance(dice_top_radius_3,0.99),"Dice_top_3_0.99%","top_results_099.csv")
catPerformance(getPerformance(tversky_top_radius_1,0.99),"Tversky_top_1_0.99%","top_results_099.csv")
catPerformance(getPerformance(tversky_top_radius_2,0.99),"Tversky_top_2_0.99%","top_results_099.csv")
catPerformance(getPerformance(tversky_top_radius_3,0.99),"Tversky_top_3_0.99%","top_results_099.csv")

cat("Combination\tBA\tsdBA\tF1\tsdF1\tMCC\tsdMCC\tsimilarityThreshold",file="proba_results_05.csv",sep="\n")
cat("Combination\tBA\tsdBA\tF1\tsdF1\tMCC\tsdMCC\tsimilarityThreshold",file="proba_results_08.csv",sep="\n")
cat("Combination\tBA\tsdBA\tF1\tsdF1\tMCC\tsdMCC\tsimilarityThreshold",file="proba_results_09.csv",sep="\n")
cat("Combination\tBA\tsdBA\tF1\tsdF1\tMCC\tsdMCC\tsimilarityThreshold",file="proba_results_1.csv",sep="\n")

catPerformance(getPerformance(tanimoto_proba_radius_1,0.8),"Tanimoto_proba_1_0.8%","proba_results_08.csv")
catPerformance(getPerformance(tanimoto_proba_radius_2,0.8),"Tanimoto_proba_2_0.8%","proba_results_08.csv")
catPerformance(getPerformance(tanimoto_proba_radius_3,0.8),"Tanimoto_proba_3_0.8%","proba_results_08.csv")
catPerformance(getPerformance(dice_proba_radius_1,0.8),"Dice_proba_1_0.8%","proba_results_08.csv")
catPerformance(getPerformance(dice_proba_radius_2,0.8),"Dice_proba_2_0.8%","proba_results_08.csv")
catPerformance(getPerformance(dice_proba_radius_3,0.8),"Dice_proba_3_0.8%","proba_results_08.csv")
catPerformance(getPerformance(tversky_proba_radius_1,0.8),"Tversky_proba_1_0.8%","proba_results_08.csv")
catPerformance(getPerformance(tversky_proba_radius_2,0.8),"Tversky_proba_2_0.8%","proba_results_08.csv")
catPerformance(getPerformance(tversky_proba_radius_3,0.8),"Tversky_proba_3_0.8%","proba_results_08.csv")

catPerformance(getPerformance(tanimoto_proba_radius_1,0.9),"Tanimoto_proba_1_0.9%","proba_results_09.csv")
catPerformance(getPerformance(tanimoto_proba_radius_2,0.9),"Tanimoto_proba_2_0.9%","proba_results_09.csv")
catPerformance(getPerformance(tanimoto_proba_radius_3,0.9),"Tanimoto_proba_3_0.9%","proba_results_09.csv")
catPerformance(getPerformance(dice_proba_radius_1,0.9),"Dice_proba_1_0.9%","proba_results_09.csv")
catPerformance(getPerformance(dice_proba_radius_2,0.9),"Dice_proba_2_0.9%","proba_results_09.csv")
catPerformance(getPerformance(dice_proba_radius_3,0.9),"Dice_proba_3_0.9%","proba_results_09.csv")
catPerformance(getPerformance(tversky_proba_radius_1,0.9),"Tversky_proba_1_0.9%","proba_results_09.csv")
catPerformance(getPerformance(tversky_proba_radius_2,0.9),"Tversky_proba_2_0.9%","proba_results_09.csv")
catPerformance(getPerformance(tversky_proba_radius_3,0.9),"Tversky_proba_3_0.9%","proba_results_09.csv")

catPerformance(getPerformance(tanimoto_proba_radius_1,0.5),"Tanimoto_proba_1_0.5%","proba_results_05.csv")
catPerformance(getPerformance(tanimoto_proba_radius_2,0.5),"Tanimoto_proba_2_0.5%","proba_results_05.csv")
catPerformance(getPerformance(tanimoto_proba_radius_3,0.5),"Tanimoto_proba_3_0.5%","proba_results_05.csv")
catPerformance(getPerformance(dice_proba_radius_1,0.5),"Dice_proba_1_0.5%","proba_results_05.csv")
catPerformance(getPerformance(dice_proba_radius_2,0.5),"Dice_proba_2_0.5%","proba_results_05.csv")
catPerformance(getPerformance(dice_proba_radius_3,0.5),"Dice_proba_3_0.5%","proba_results_05.csv")
catPerformance(getPerformance(tversky_proba_radius_1,0.5),"Tversky_proba_1_0.5%","proba_results_05.csv")
catPerformance(getPerformance(tversky_proba_radius_2,0.5),"Tversky_proba_2_0.5%","proba_results_05.csv")
catPerformance(getPerformance(tversky_proba_radius_3,0.5),"Tversky_proba_3_0.5%","proba_results_05.csv")

catPerformance(getPerformance(tanimoto_proba_radius_1,0.99),"Tanimoto_proba_1_0.99%","proba_results_099.csv")
catPerformance(getPerformance(tanimoto_proba_radius_2,0.99),"Tanimoto_proba_2_0.99%","proba_results_099.csv")
catPerformance(getPerformance(tanimoto_proba_radius_3,0.99),"Tanimoto_proba_3_0.99%","proba_results_099.csv")
catPerformance(getPerformance(dice_proba_radius_1,0.99),"Dice_proba_1_0.99%","proba_results_099.csv")
catPerformance(getPerformance(dice_proba_radius_2,0.99),"Dice_proba_2_0.99%","proba_results_099.csv")
catPerformance(getPerformance(dice_proba_radius_3,0.99),"Dice_proba_3_0.99%","proba_results_099.csv")
catPerformance(getPerformance(tversky_proba_radius_1,0.99),"Tversky_proba_1_0.99%","proba_results_099.csv")
catPerformance(getPerformance(tversky_proba_radius_2,0.99),"Tversky_proba_2_0.99%","proba_results_099.csv")
catPerformance(getPerformance(tversky_proba_radius_3,0.99),"Tversky_proba_3_0.99%","proba_results_099.csv")


t1_t2_ref <- data.frame(
  name=factor(c("RF","SVM","DNN","Tversky 99%","Tanim. 50%"),levels = c("RF","SVM","DNN","Tversky 99%","Tanim. 50%")),
  BA=c(0.75,0.85,0.88,0.91011,0.96934),
  F1=c(0.66,0.8,0.8,0.97213,0.99194),
  MCC=c(0.68,0.79,0.77,0.76161,0.94208)
)
#  sdBA=c(0.02,0.02,0.02,0.01643,0.0114),
#  sdF1=c(0.04,0.03,0.02,0.00246,0.00252),
#  sdMCC=c(0.03,0.03,0.02,0.02119,0.01739)



t1_t2_ref.m <- melt(t1_t2_ref, id.vars='name')
t1_t2_ref.m=cbind(t1_t2_ref.m,rep(0,length(t1_t2_ref.m[,1])))
t1_t2_ref.m=cbind(t1_t2_ref.m,rep(0,length(t1_t2_ref.m[,1])))
colnames(t1_t2_ref.m)[4]="sd"
colnames(t1_t2_ref.m)[5]="randomsd"
t1_t2_ref.m[4]=c(0.02,0.02,0.02,0.01643,0.0114,0.04,0.03,0.02,0.00246,0.00252,0.03,0.03,0.02,0.02119,0.01739)



p=ggplot(t1_t2_ref.m,aes(name,value,fill = variable)) +
  geom_bar( position = "dodge", stat="identity") + 
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,position=position_dodge(.9)) +
  geom_segment(aes(x = 0.55, xend = 0.85, y = 0.49599, yend = 0.49599),col="black",linetype = "dashed") +
  geom_segment(aes(x = 0.85, xend = 1.15, y = 0.87866, yend = 0.87866),col="black",linetype = "dashed") +
  geom_segment(aes(x = 1.15, xend = 1.45, y = 0.0, yend = 0.0),col="black",linetype = "dashed") +

  geom_segment(aes(x = 1.55, xend = 1.85, y = 0.49599, yend = 0.49599),col="black",linetype = "dashed") +
  geom_segment(aes(x = 1.85, xend = 2.15, y = 0.87866, yend = 0.87866),col="black",linetype = "dashed") +
  geom_segment(aes(x = 2.15, xend = 2.45, y = 0.0, yend = 0.0),col="black",linetype = "dashed") +

  geom_segment(aes(x = 2.55, xend = 2.85, y = 0.49599, yend = 0.49599),col="black",linetype = "dashed") +
  geom_segment(aes(x = 2.85, xend = 3.15, y = 0.87866, yend = 0.87866),col="black",linetype = "dashed") +
  geom_segment(aes(x = 3.15, xend = 3.45, y = 0.0, yend = 0.0),col="black",linetype = "dashed") +

  geom_segment(aes(x = 3.55, xend = 3.85, y = 0.49599, yend = 0.49599),col="black",linetype = "dashed") +
  geom_segment(aes(x = 3.85, xend = 4.15, y = 0.87866, yend = 0.87866),col="black",linetype = "dashed") +
  geom_segment(aes(x = 4.15, xend = 4.45, y = 0.0, yend = 0.0),col="black",linetype = "dashed") +

  geom_segment(aes(x = 4.55, xend = 4.85, y = 0.49599, yend = 0.49599),col="black",linetype = "dashed") +
  geom_segment(aes(x = 4.85, xend = 5.15, y = 0.87866, yend = 0.87866),col="black",linetype = "dashed") +
  geom_segment(aes(x = 5.15, xend = 5.45, y = 0.0, yend = 0.0),col="black",linetype = "dashed") +

  labs(y = "Metric value",
                x = "Prediction Method",colour="",fill="",title="Type I vs type II") + 
                scale_colour_Publication() +
                theme_Publication() + 
                ylim(0,1.0)
#random sampling ba:     0.49599 0.01135 f1:     0.87866 0.00464 mcc     -0.00788        0.02224 found:  1012    notFound        0
print(p)


ggsave(file="comparisonBarplot.svg", dpi=300,plot=p, width=8, height=7)

#Tversky data taken above threshold 0.3
#Tanimoto data taken above threshold 0.4



t1_t2_ref_BA=
t1_t2_ref_F1=
t1_t2_ref_MCC=


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