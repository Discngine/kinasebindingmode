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



prepareData<-function(inputData,returnRandom,curDateThreshold){

    levels(inputData$Correct) <- c(FALSE,TRUE)
    inputData$Correct <- as.logical(inputData$Correct)
    levels(inputData$randomCorrect) <- c(FALSE,TRUE)
    inputData$randomCorrect <- as.logical(inputData$randomCorrect)

    
    r=inputData[inputData$dateThreshold==curDateThreshold,]
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
    print("here")
    for(i in seq(1)){
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
    data=cbind(data,rep(curDateThreshold,length(t1_2_sd)))

    colnames(data)=c("cutoff","mean","sd","n_mean","n_sd","type1_mean","type1_sd","type2_mean","type2_sd","type1_2_mean","type1_2_sd","date")

    rdata=as.data.frame(cutoffs)
    rdata=cbind(rdata,racc_mean)
    rdata=cbind(rdata,racc_sd)
    colnames(rdata)=c("cutoff","mean","sd")

    if(returnRandom) {
        return(rdata)
    }
    return(data)
}




generatePlot<-function(fname,label){

    r=read.table(fname,h=T)
    uniqueDates=unique(r$dateThreshold)

    thresholds=seq(0.1,1.0,0.05)
    for(th in thresholds){
        threshold=as.character(th)
    
        n=0
        results=c()
        for(currentDate in uniqueDates){
            mfp1=prepareData(read.table(fname,h=T),0,currentDate)
            line=mfp1[mfp1$cutoff==threshold,]
            if(n==0) results=line
            else {
                results=rbind(results,line)
            }
            n=n+1
        }

        outname=sprintf("images/%s Threshold %s.svg",label,threshold)

        results$date=as.Date(results$date)

        mfp3=results
        p=ggplot(results, aes(x=date, y=type1_mean)) + 
        geom_line(aes(colour='Type I'),size=2)+
        geom_line(aes(y=type2_mean,colour='Type II'),size=2)+
        geom_line(aes(y=type1_2_mean,colour='Type I 1/2'),size=2)+
        labs(y = "MCC",
                        x = "Date",colour="",title=label) + 
                        scale_colour_Publication2() +  
                        theme_Publication() + 
                        scale_y_continuous(breaks=seq(0.0,1.0,0.1),limits=c(0.0,1.0)) 
        ggsave(file=outname, plot=p, width=8, height=8)
    }
}


fingerprints=c("Tanimoto","Dice","Tversky")
radii=c(1,2,3)

for(fp in fingerprints){
    for(radius in radii){
        fname=sprintf("predictivity_radius_%d_%s.csv",radius,tolower(fp))
        label=sprintf("%s Similarity MFP%d",fp,radius)
        generatePlot(fname,label)
    }
}

# fname="predictivity_radius_3_tanimoto.csv"
# outname="timelinePredictionTanimotoECFP6.svg"
# label="Tanimoto Similarity MFP3"
# threshold="0.45"
# generatePlot(fname,outname,label,threshold)

# fname="predictivity_radius_2_tanimoto.csv"
# outname="timelinePredictionTanimotoECFP4.svg"
# label="Tanimoto Similarity MFP2"
# threshold="0.45"
# generatePlot(fname,outname,label,threshold)

# fname="predictivity_radius_1_tanimoto.csv"
# outname="timelinePredictionTanimotoECFP2.svg"
# label="Tanimoto Similarity MFP1"
# threshold="0.45"
# generatePlot(fname,outname,label,threshold)

# fname="predictivity_radius_3_dice.csv"
# outname="timelinePredictionDiceECFP6.svg"
# label="Dice Similarity MFP3"
# threshold="0.45"
# generatePlot(fname,outname,label,threshold)
