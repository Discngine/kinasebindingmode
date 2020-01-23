#library(ROCR)
#library(tidyr)
library(ggplot2)


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
    cutoffs=seq(1.0,0.0,-0.01)
    rfinal_accs=data.frame()#random final accs
    final_accs=data.frame()
    for(i in seq(10)){
        accs=c()
        raccs=c()
        for(cutoff in cutoffs){
            print(sum(r$Cycle==(i-1)))
            sl <- append(sl,list(r$Similarity[r$Cycle==(i-1)]))
            ll<-append(ll,list(r$Correct[r$Cycle==(i-1)]))
        
            stmp=r$Similarity[r$Cycle==(i-1)]
            rstmp=r$randomSimilarity[r$Cycle==(i-1)]
            ctmp=r$Correct[r$Cycle==(i-1)]
            rctmp=r$randomCorrect[r$Cycle==(i-1)]
            
            #if(cutoff==1){
                #print(ctmp[stmp>=cutoff])
            
            acc=sum(ctmp[stmp>=cutoff]==T)/length(ctmp[stmp>=cutoff])
            racc=sum(rctmp[stmp>=cutoff]==T)/length(rctmp[stmp>=cutoff])
            
            accs=c(accs,acc)
            raccs=c(raccs,racc)

            #}
        }
        if(i==1){
        #         plot(cutoffs,accs,ylim=c(0.7,1.0),ty="l")
        #         par(new=T)
                final_accs=as.data.frame(accs)
                final_raccs=as.data.frame(raccs)

        }
        else{
            final_accs=cbind(final_accs,accs)
            final_raccs=cbind(final_raccs,raccs)
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


    data=as.data.frame(cutoffs)
    data=cbind(data,acc_mean)
    data=cbind(data,acc_sd)
    colnames(data)=c("cutoff","mean","sd")

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
                ylim(0.4,1.0) +
                scale_y_continuous(breaks=seq(0.4,1.0,0.1)) +
                scale_x_continuous(breaks=seq(0.0,1.0,0.1)) 
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
                ylim(0.4,1.0) +
                scale_y_continuous(breaks=seq(0.4,1.0,0.1)) +
                scale_x_continuous(breaks=seq(0.0,1.0,0.1)) 
print(p)
ggsave(file="generalPredictionTanimoto.svg", plot=p, width=8, height=8)





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
                scale_y_continuous(breaks=seq(0.4,1.0,0.1)) +
                scale_x_continuous(breaks=seq(0.0,1.0,0.1)) +
                ylim(0.4,1.0)
print(p)
ggsave(file="generalPredictionTversky.svg", plot=p, width=8, height=8)