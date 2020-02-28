library(ggplot2)
library(dplyr)
library(gridExtra)
library(reshape2)

for(i in 2700:2700) {
    mydataName <- paste0("~/cloud/project/hpc_expTrueRM14A9A/exph_trueRM14A9A_esum_",formatC(i,width=4,flag="0"),".RData")
    load(mydataName)
}

tmpp <- exph_trueRM14A9A_esum_2700

godlist0 <- which(tmpp$ecisL > 0 | tmpp$ecisH <0)
godlist0N <- (1:4710)[-godlist0]

sigdf <- tmpp[godlist0,]
sort_sigdf <- sort(abs(sigdf$log.ecis),decreasing = T)

mydf <- data.frame(prop = 1:length(sort_sigdf) / 3308, absecis = sort_sigdf )

ggplot(mydf) +
  geom_point(aes(x = absecis, y= prop), size= 0.5) + 
  #xlim=rev(range(sort_sigdf)) + 
  labs(title="Cumulative distribution of \n significantly cis-affected genes", y="cumulative proportion of genes", x=expression(paste("|log"["2"], "(","e"["cis"],")|"))) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=0, hjust=1)) +
  scale_x_reverse(breaks= seq(6,0,by=-0.5)) + 
  scale_y_continuous(breaks = seq(0,1,by=0.05)) + 
  theme(
    plot.title = element_text(size = 14, face = "bold",hjust=.5),
    panel.grid = element_line(color="grey85"),  
    axis.text.x = element_text(colour="grey20",size=11,hjust=.5,vjust=.5,face="plain"),
    axis.text.y = element_text(colour="grey20",size=11,hjust=.5,vjust=.5,face="plain"),
    axis.title.x = element_text(colour="grey20",size=13,hjust=.5,vjust=0,face="bold"), 
    axis.title.y = element_text(colour="grey20",size=13,hjust=.5,vjust=.5,face="bold"))

ggsave("~/cloud/project/figures/fig_6B_hist.pdf", device = "pdf", width = 14, height = 12, units = "cm",dpi=200)
     
     