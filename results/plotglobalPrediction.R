library(ROCR)
library(tidyr)

r=read.table("globalPrediction_radius_1_tanimoto.csv",h=T)
#r[r[r[,"Correct"]=="True"],"True"]=T
before2=r$Correct
levels(r$Correct) <- c(FALSE,TRUE)
before=r$Correct
r$Correct <- as.logical(r$Correct)
after=r$Correct
sl=list()
ll=list()
for(i in seq(10)){
    print(sum(r$Cycle==(i-1)))
    sl <- append(sl,list(r$Similarity[r$Cycle==(i-1)]))
    ll<-append(ll,list(r$Correct[r$Cycle==(i-1)]))
    for(cutoff in seq(1.0,0.0,-0.1)){
        stmp=r$Similarity[r$Cycle==(i-1)]
        ctmp=r$Correct[r$Cycle==(i-1)]
        print(cutoff)
        #if(cutoff==1){
            #print(ctmp[stmp>=cutoff])
            acc=sum(ctmp[stmp>=cutoff]==T)/length(ctmp[stmp>=cutoff])
            print(acc)
        #}
    }
}

pred=prediction(sl,ll)
perf <- performance(pred,"tpr","fpr")
plot(perf,colorize=T)
#plot(perf,avg= "vertical", spread.estimate="boxplot", show.spread.at= seq(0.1, 0.9, by=0.1))