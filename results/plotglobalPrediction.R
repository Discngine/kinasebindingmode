#library(ROCR)
#library(tidyr)
library(ggplot2)
library(ggrepel)
library(formattable)
library(mccr)

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

scale_colour_Publication2 <- function(...){
      library(scales)
      discrete_scale("colour","Publication",manual_pal(values = c("#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)

}



prepareData<-function(r,returnRandom){
#r[r[r[,"Correct"]=="True"],"True"]=T

    levels(r$Correct) <- c(FALSE,TRUE)

    r$Correct <- as.logical(r$Correct)
    levels(r$randomCorrect) <- c(FALSE,TRUE)
    r$randomCorrect <- as.logical(r$randomCorrect)


    sl=list()
    ll=list()
    accs=c()
    raccs=c() #random accs
    npred=c()
    cutoffs=seq(1.0,0.0,-0.01)
    rfinal_accs=data.frame()#random final accs
    final_accs=data.frame()
    final_npred=data.frame()
    final_type1=data.frame()
    final_type2=data.frame()
    final_type1_2=data.frame()
    for(i in seq(10)){
        accs=c()
        raccs=c()
        ns=c()
        t1=c()
        t2=c()
        t1_2=c()
        print(i)
        for(cutoff in cutoffs){
            sl <- append(sl,list(r$Similarity[r$Cycle==(i-1)]))
            ll<-append(ll,list(r$Correct[r$Cycle==(i-1)]))
        
            stmp=r$Similarity[r$Cycle==(i-1)]
            rstmp=r$randomSimilarity[r$Cycle==(i-1)]
            ctmp=r$Correct[r$Cycle==(i-1)]
            rctmp=r$randomCorrect[r$Cycle==(i-1)]
            cactual=r$ActualClass[r$Cycle==(i-1)]
            cpred=r$Prediction[r$Cycle==(i-1)]
            cpredSubset=cpred[stmp>=cutoff]
            actualSubSet=cactual[stmp>=cutoff]
            actualSubSetCorrect=ctmp[stmp>=cutoff]

            predMask=cpredSubset=="type1"
            trueMask=actualSubSet=="type1"
            type1_mcc=mccr(trueMask,predMask)

            predMask=cpredSubset=="type2"
            trueMask=actualSubSet=="type2"
            type2_mcc=mccr(trueMask,predMask)

            predMask=cpredSubset=="type1_2"
            trueMask=actualSubSet=="type1_2"
            type1_2_mcc=mccr(trueMask,predMask)

            acc=sum(ctmp[stmp>=cutoff]==T)/length(ctmp[stmp>=cutoff])
            racc=sum(rctmp[stmp>=cutoff]==T)/length(rctmp[stmp>=cutoff])
            
            accs=c(accs,acc)
            raccs=c(raccs,racc)
            t1=c(t1,type1_mcc)
            t2=c(t2,type2_mcc)
            t1_2=c(t1_2,type1_2_mcc)
            ns=c(ns,length(ctmp[stmp>=cutoff])/length(ctmp))
            #}
        }
        if(i==1){
        #         plot(cutoffs,accs,ylim=c(0.7,1.0),ty="l")
        #         par(new=T)
                final_accs=as.data.frame(accs)
                final_raccs=as.data.frame(raccs)
                final_npred=as.data.frame(ns)
                final_type1=as.data.frame(t1)
                final_type2=as.data.frame(t2)
                final_type1_2=as.data.frame(t1_2)

        }
        else{
            final_accs=cbind(final_accs,accs)
            final_raccs=cbind(final_raccs,raccs)
            final_npred=cbind(final_npred,ns)
            final_type1=cbind(final_type1,t1)
            final_type2=cbind(final_type2,t2)
            final_type1_2=cbind(final_type1_2,t1_2)
        #    plot(cutoffs,accs,ylim=c(0.7,1.0),ty="l")
        #     par(new=T)
        }
            
        #print(accs)
    }

    #accs_matrix=matrix(accs,nrow=length(cutoffs),byrow=T)
    #print(accs_matrix)

    acc_mean=apply(final_accs,1,mean)
    acc_sd=apply(final_accs,1,sd)
    racc_mean=apply(final_raccs,1,mean)
    racc_sd=apply(final_raccs,1,sd)
    npred_mean=apply(final_npred,1,mean)
    npred_sd=apply(final_npred,1,sd)
    t1_mean=apply(final_type1,1,mean)
    t1_sd=apply(final_type1,1,sd)
    t2_mean=apply(final_type2,1,mean)
    t2_sd=apply(final_type2,1,sd)
    t1_2_mean=apply(final_type1_2,1,mean)
    t1_2_sd=apply(final_type1_2,1,sd)


    data=as.data.frame(cutoffs)
    data=cbind(data,acc_mean)
    data=cbind(data,acc_sd)
    data=cbind(data,npred_mean)
    data=cbind(data,npred_sd)
    data=cbind(data,t1_mean)
    data=cbind(data,t1_sd)
    data=cbind(data,t2_mean)
    data=cbind(data,t2_sd)
    data=cbind(data,t1_2_mean)
    data=cbind(data,t1_2_sd)

    colnames(data)=c("cutoff","mean","sd","n_mean","n_sd","type1_mean","type1_sd","type2_mean","type2_sd","type1_2_mean","type1_2_sd")

    rdata=as.data.frame(cutoffs)
    rdata=cbind(rdata,racc_mean)
    rdata=cbind(rdata,racc_sd)
    colnames(rdata)=c("cutoff","mean","sd")

    if(returnRandom) {
        return(rdata)
    }
    return(data)
}


mfp1=prepareData(read.table("globalPrediction_radius_1_dice.csv",h=T),0)
mfp2=prepareData(read.table("globalPrediction_radius_2_dice.csv",h=T),0)
mfp3=prepareData(read.table("globalPrediction_radius_3_dice.csv",h=T),0)
rfp1=prepareData(read.table("globalPrediction_radius_1_dice.csv",h=T),1)
rfp2=prepareData(read.table("globalPrediction_radius_2_dice.csv",h=T),1)
rfp3=prepareData(read.table("globalPrediction_radius_3_dice.csv",h=T),1)

breaks=c(5,10,20,49,50,74,75,90)
print(round(mfp1[,"n_mean"]*100.0))
p=ggplot(mfp1, aes(x=cutoff, y=mean)) + 
  geom_line(aes(ymin=0.7,ymax=1.0,colour='MFP1'),size=2)+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd,colour='MFP1')) +
  geom_line(data=mfp2,aes(ymin=0.7,ymax=1.0,colour='MFP2'),size=2)+
  geom_errorbar(data=mfp2,aes(ymin=mean-sd, ymax=mean+sd,colour='MFP2')) +
  geom_line(data=mfp3,aes(ymin=0.7,ymax=1.0,colour='MFP3'),size=2)+
  geom_errorbar(data=mfp3,aes(ymin=mean-sd, ymax=mean+sd,colour='MFP3')) +
  geom_line(data=rfp3,aes(ymin=0.7,ymax=1.0,colour='Random'),size=1)+
  geom_errorbar(data=rfp3,aes(ymin=mean-sd, ymax=mean+sd,colour='Random')) +
  geom_hline(  yintercept=0.95, linetype="dashed", "black", size=1) +
  labs(y = "Accuracy",
                x = "Similarity Threshold",colour="",title="Dice Similarity") + 
                scale_colour_Publication() +  
                theme_Publication() + 
                scale_y_continuous(breaks=seq(0.4,1.0,0.1),limits=c(0.4,1.0)) +
                scale_x_continuous(breaks=seq(0.3,1.0,0.1),limits=c(0.3,1.0)) 
                
                
                
print(p)
ggsave(file="generalPredictionDice.svg", plot=p, width=8, height=8)




mfp1=prepareData(read.table("globalPrediction_radius_1_tanimoto.csv",h=T),0)
mfp2=prepareData(read.table("globalPrediction_radius_2_tanimoto.csv",h=T),0)
mfp3=prepareData(read.table("globalPrediction_radius_3_tanimoto.csv",h=T),0)
rfp1=prepareData(read.table("globalPrediction_radius_1_tanimoto.csv",h=T),1)
rfp2=prepareData(read.table("globalPrediction_radius_2_tanimoto.csv",h=T),1)
rfp3=prepareData(read.table("globalPrediction_radius_3_tanimoto.csv",h=T),1)

p=ggplot(mfp1, aes(x=cutoff, y=mean)) + 
  geom_line(aes(ymin=0.7,ymax=1.0,colour='MFP1'),size=2)+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd,colour='MFP1')) +
  geom_line(data=mfp2,aes(ymin=0.7,ymax=1.0,colour='MFP2'),size=2)+
  geom_errorbar(data=mfp2,aes(ymin=mean-sd, ymax=mean+sd,colour='MFP2')) +
  geom_line(data=mfp3,aes(ymin=0.7,ymax=1.0,colour='MFP3'),size=2)+
  geom_errorbar(data=mfp3,aes(ymin=mean-sd, ymax=mean+sd,colour='MFP3')) +
  geom_line(data=rfp3,aes(ymin=0.7,ymax=1.0,colour='Random'),size=1)+
  geom_errorbar(data=rfp3,aes(ymin=mean-sd, ymax=mean+sd,colour='Random')) +
  geom_hline(  yintercept=0.95, linetype="dashed", "black", size=1) +
  labs(y = "Accuracy",
                x = "Similarity Threshold",colour="",title="Tanimoto Similarity") + 
                scale_colour_Publication() +  
                theme_Publication() + 
                scale_y_continuous(breaks=seq(0.4,1.0,0.1),limits=c(0.4,1.0)) +
                scale_x_continuous(breaks=seq(0.3,1.0,0.1),limits=c(0.3,1.0)) +
                geom_text_repel(
                    data = subset(mfp3, (ceiling(n_mean*100) %in% breaks )),
                    aes(label = percent(n_mean)),
                    size = 5,
                    box.padding = unit(1.35, "lines"),
                    point.padding = unit(0.3, "lines")
                )
print(p)
ggsave(file="generalPredictionTanimoto_Accuracy.svg", plot=p, width=8, height=8)


print(mfp3)

p=ggplot(mfp3, aes(x=cutoff, y=type1_mean)) + 
  geom_line(aes(ymin=0.7,ymax=1.0,colour='Type I'),size=2)+
  geom_errorbar(data=mfp3,aes(ymin=type1_mean-type1_sd, ymax=type1_mean+type1_sd,colour='Type I'),position=position_jitter(width=0.001),) +
  geom_line(aes(y=type2_mean,ymin=0.7,ymax=1.0,colour='Type II'),size=2)+
  geom_errorbar(data=mfp3,aes(ymin=type2_mean-type2_sd, ymax=type2_mean+type2_sd,colour='Type II'),position=position_jitter(width=0.001),) +
  geom_line(aes(y=type1_2_mean,ymin=0.7,ymax=1.0,colour='Type I 1/2'),size=2)+
  geom_errorbar(data=mfp3,aes(ymin=type1_2_mean-type1_2_sd, ymax=type1_2_mean+type1_2_sd,colour='Type I 1/2'),position=position_jitter(width=0.001),) +
  labs(y = "MCC",
                x = "Similarity Threshold",colour="",title="Tanimoto Similarity MFP3") + 
                scale_colour_Publication2() +  
                theme_Publication() + 
                scale_y_continuous(breaks=seq(0.65,1.0,0.1),limits=c(0.65,1.0)) +
                scale_x_continuous(breaks=seq(0.3,1.0,0.1),limits=c(0.3,1.0)) +
                geom_text_repel(
                    data = subset(mfp3, (ceiling(n_mean*100) %in% breaks )),
                    aes(label = percent(n_mean)),
                    size = 5,
                    box.padding = unit(4.35, "lines"),
                    point.padding = unit(0.3, "lines")
                )
print(p)
ggsave(file="generalPredictionTanimotoECFP6.svg", plot=p, width=8, height=8)





p=ggplot(mfp2, aes(x=cutoff, y=type1_mean)) + 
  geom_line(aes(ymin=0.7,ymax=1.0,colour='Type I'),size=2)+
  geom_errorbar(data=mfp2,aes(ymin=type1_mean-type1_sd, ymax=type1_mean+type1_sd,colour='Type I'),position=position_jitter(width=0.001),) +
  geom_line(aes(y=type2_mean,ymin=0.7,ymax=1.0,colour='Type II'),size=2)+
  geom_errorbar(data=mfp2,aes(ymin=type2_mean-type2_sd, ymax=type2_mean+type2_sd,colour='Type II'),position=position_jitter(width=0.001),) +
  geom_line(aes(y=type1_2_mean,ymin=0.7,ymax=1.0,colour='Type I 1/2'),size=2)+
  geom_errorbar(data=mfp2,aes(ymin=type1_2_mean-type1_2_sd, ymax=type1_2_mean+type1_2_sd,colour='Type I 1/2'),position=position_jitter(width=0.001),) +
  labs(y = "MCC",
                x = "Similarity Threshold",colour="",title="Tanimoto Similarity") + 
                scale_colour_Publication2() +  
                theme_Publication() + 
                scale_y_continuous(breaks=seq(0.65,1.0,0.1),limits=c(0.65,1.0)) +
                scale_x_continuous(breaks=seq(0.3,1.0,0.1),limits=c(0.3,1.0)) +
                geom_text_repel(
                    data = subset(mfp2, (ceiling(n_mean*100) %in% breaks )),
                    aes(label = percent(n_mean)),
                    size = 5,
                    box.padding = unit(2.35, "lines"),
                    point.padding = unit(0.3, "lines")
                )
print(p)
ggsave(file="generalPredictionTanimotoECFP4.svg", plot=p, width=8, height=8)




p=ggplot(mfp1, aes(x=cutoff, y=type1_mean)) + 
  geom_line(aes(ymin=0.7,ymax=1.0,colour='Type I'),size=2)+
  geom_errorbar(data=mfp1,aes(ymin=type1_mean-type1_sd, ymax=type1_mean+type1_sd,colour='Type I'),position=position_jitter(width=0.001),) +
  geom_line(aes(y=type2_mean,ymin=0.7,ymax=1.0,colour='Type II'),size=2)+
  geom_errorbar(data=mfp1,aes(ymin=type2_mean-type2_sd, ymax=type2_mean+type2_sd,colour='Type II'),position=position_jitter(width=0.001),) +
  geom_line(aes(y=type1_2_mean,ymin=0.7,ymax=1.0,colour='Type I 1/2'),size=2)+
  geom_errorbar(data=mfp1,aes(ymin=type1_2_mean-type1_2_sd, ymax=type1_2_mean+type1_2_sd,colour='Type I 1/2'),position=position_jitter(width=0.001),) +
  labs(y = "MCC",
                x = "Similarity Threshold",colour="",title="Tanimoto Similarity") + 
                scale_colour_Publication2() +  
                theme_Publication() + 
                scale_y_continuous(breaks=seq(0.65,1.0,0.1),limits=c(0.65,1.0)) +
                scale_x_continuous(breaks=seq(0.3,1.0,0.1),limits=c(0.3,1.0)) +
                geom_text_repel(
                    data = subset(mfp1, (ceiling(n_mean*100) %in% breaks )),
                    aes(label = percent(n_mean)),
                    size = 5,
                    box.padding = unit(2.35, "lines"),
                    point.padding = unit(0.3, "lines")
                )
print(p)
ggsave(file="generalPredictionTanimotoECFP2.svg", plot=p, width=8, height=8)




mfp1=prepareData(read.table("globalPrediction_radius_1_tversky.csv",h=T),0)
mfp2=prepareData(read.table("globalPrediction_radius_2_tversky.csv",h=T),0)
mfp3=prepareData(read.table("globalPrediction_radius_3_tversky.csv",h=T),0)
rfp1=prepareData(read.table("globalPrediction_radius_1_tversky.csv",h=T),1)
rfp2=prepareData(read.table("globalPrediction_radius_2_tversky.csv",h=T),1)
rfp3=prepareData(read.table("globalPrediction_radius_3_tversky.csv",h=T),1)

p=ggplot(mfp1, aes(x=cutoff, y=mean)) + 
  geom_line(aes(ymin=0.7,ymax=1.0,colour='MFP1'),size=2)+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd,colour='MFP1')) +
  geom_line(data=mfp2,aes(ymin=0.7,ymax=1.0,colour='MFP2'),size=2)+
  geom_errorbar(data=mfp2,aes(ymin=mean-sd, ymax=mean+sd,colour='MFP2')) +
  geom_line(data=mfp3,aes(ymin=0.7,ymax=1.0,colour='MFP3'),size=2)+
  geom_errorbar(data=mfp3,aes(ymin=mean-sd, ymax=mean+sd,colour='MFP3')) +
  geom_line(data=rfp3,aes(ymin=0.7,ymax=1.0,colour='Random'),size=2)+
  geom_errorbar(data=rfp3,aes(ymin=mean-sd, ymax=mean+sd,colour='Random')) +
  geom_hline(  yintercept=0.95, linetype="dashed", "black", size=1) +
  labs(y = "Accuracy",
                x = "Similarity Threshold",colour="",title="Tversky Similarity") + 
                scale_colour_Publication() +  
                theme_Publication() + 
                scale_y_continuous(breaks=seq(0.4,1.0,0.1),limits=c(0.4,1.0)) +
                scale_x_continuous(breaks=seq(0.3,1.0,0.1),limits=c(0.3,1.0))
print(p)
ggsave(file="generalPredictionTversky.svg", plot=p, width=8, height=8)

