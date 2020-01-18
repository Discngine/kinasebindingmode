library(ggplot2)
library(gridExtra)
library(reshape2)

r=read.table("abbvie_stats_I_II_50.csv",header=TRUE,sep="\t")
r2=read.table("abbvie_stats_I_1_2_II_50.csv",h=T,sep="\t")
r3=read.table("abbvie_stats_I_I_1_2_50.csv",h=T,sep="\t")
data=rbind(r,r2)
data=rbind(data,r3)
#meandata=apply(r,2,mean)
#sddata=apply(r,2,sd)
#data=melt(meandata, id.vars='name')

theme_Publication <- function(base_size=14, base_family="Roboto") {
      library(grid)
      library(ggthemes)
      (theme_foundation(base_size=base_size, base_family=base_family)
       + theme(plot.title = element_text(face = "bold",
                                         size = rel(5.2), hjust = 0.5),
              line= element_line(size=rel(1)),
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

#data.m <- melt(data, id.vars='Dataset')
#print(data.m)
# p=ggplot(data.m,aes(Dataset,value,fill = variable)) +
#   geom_boxplot() +
#   labs(y = "Metric value",
#         x = "Prediction Method",colour="",fill="",title="External validation") + 
#         scale_colour_Publication() +
#         theme_Publication() + 
#         ylim(0.0,1.0)



data.m <- melt(data, id.vars='Dataset')
print(data.m)
p=ggplot(data.m,aes(Dataset,value,fill = variable)) +
  geom_boxplot() +
  geom_boxplot(aes(color = variable),
               fatten = NULL, fill = NA, coef = 0, outlier.alpha = 0,
               show.legend = F)+
  labs(y = "Metric value",
        x = "Prediction Method",colour="",fill="",title="External validation") + 
        theme_Publication() + 
        ylim(0.0,1.0)

  #geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,position=position_dodge(.9)) +
  #geom_segment(aes(x = 0.55, xend = 0.85, y = 0.49599, yend = 0.49599),col="black",linetype = "dashed") +

#print(p)

ggsave(file="50_abbvie.svg", plot=p, width=8, height=8)
