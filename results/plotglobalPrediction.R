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




r=read.table("globalPrediction_radius_1_tanimoto.csv",h=T)
#r[r[r[,"Correct"]=="True"],"True"]=T
before2=r$Correct
levels(r$Correct) <- c(FALSE,TRUE)
before=r$Correct
r$Correct <- as.logical(r$Correct)
after=r$Correct
sl=list()
ll=list()
accs=c()
cutoffs=seq(1.0,0.0,-0.01)
final_accs=data.frame()
for(i in seq(10)){
    accs=c()
    for(cutoff in cutoffs){
        print(sum(r$Cycle==(i-1)))
        sl <- append(sl,list(r$Similarity[r$Cycle==(i-1)]))
        ll<-append(ll,list(r$Correct[r$Cycle==(i-1)]))
    
        stmp=r$Similarity[r$Cycle==(i-1)]
        ctmp=r$Correct[r$Cycle==(i-1)]
        
        #if(cutoff==1){
            #print(ctmp[stmp>=cutoff])
        
        acc=sum(ctmp[stmp>=cutoff]==T)/length(ctmp[stmp>=cutoff])
        accs=c(accs,acc)
       
        #}
    }
     if(i==1){
    #         plot(cutoffs,accs,ylim=c(0.7,1.0),ty="l")
    #         par(new=T)
             final_accs=as.data.frame(accs)
     }
     else{
         final_accs=cbind(final_accs,accs)
    #    plot(cutoffs,accs,ylim=c(0.7,1.0),ty="l")
    #     par(new=T)
     }
        
    #print(accs)
}

#accs_matrix=matrix(accs,nrow=length(cutoffs),byrow=T)
#print(accs_matrix)

acc_mean=apply(final_accs,1,mean)
acc_sd=apply(final_accs,1,sd)

data=as.data.frame(cutoffs)
data=cbind(data,acc_mean)
data=cbind(data,acc_sd)
colnames(data)=c("cutoff","mean","sd")

#plot(cutoffs,acc_mean,ylim=c(0.7,1.0),ty="l",lwd=3,col="red")
#ggplot(df2, aes(x=dose, y=len, group=supp, color=supp)) + 
#geom_pointrange(aes(ymin=len-sd, ymax=len+sd))
# Use geom_line()+geom_pointrange()
p=ggplot(data, aes(x=cutoff, y=mean)) + 
  geom_line(color="red",aes(ymin=0.7,ymax=1.0))+
  geom_pointrange(aes(ymin=mean-sd, ymax=mean+sd),color="red") +
  scale_colour_Publication() +
                theme_Publication() 
print(p)
#pred=prediction(sl,ll)
#perf <- performance(pred,"tpr","fpr")
#plot(perf,colorize=T)
#plot(perf,avg= "vertical", spread.estimate="boxplot", show.spread.at= seq(0.1, 0.9, by=0.1))