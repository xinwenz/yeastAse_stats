load("~/cloud/project/M6_expdc_ci_hybrid/expdc_hybrid_ci_450.RData")

tmpp <- expdc_hybrid_ci_450
#hist(tmpp$log.ecis, breaks=100,xlim=c(-3,3)) ## histgram 
ggplot(tmpp, aes(x=log.ecis)) + 
    geom_histogram(binwidth = 0.05, color="grey60",fill='grey80') + 
    theme_bw() +
    labs(title="Distribution of cis effect", y="number of genes", x=expression(paste("log"["2"], "(","e"["cis"],")"))) + 
    theme_bw() +
    theme(axis.text.x=element_text(angle=0, hjust=1)) +
    #coord_cartesian(ylim=c(0, 0.6)) +
    scale_x_continuous(breaks=seq(-4,4,by=0.5)) + 
    scale_y_continuous(breaks = seq(0,500,by=100)) + 
    coord_cartesian(xlim=c(-4, 4)) +
    theme(
        panel.grid = element_line(color="grey85"),  
        #legend.position = c(0.875, 0.9),
        #legend.background = element_blank(),
        #legend.title = element_blank(),
        #legend.text = element_text(size=14,face="plain"),
        plot.title = element_text(size = 18, face = "bold",hjust=.5),
        axis.text.x = element_text(angle=40,colour="grey20",size=8,hjust=.5,vjust=.5,face="bold"),
        axis.text.y = element_text(colour="grey20",size=11,hjust=.5,vjust=.5,face="bold"),
        axis.title.x = element_text( colour="grey20",size=15,hjust=.5,vjust=0,face="bold"), 
        axis.title.y = element_text(colour="grey20",size=15,hjust=.5,vjust=.9,face="bold"))

ggsave("~/cloud/project/figures/fig_hist.pdf", device = "pdf", width = 12, height = 12, units = "cm",dpi=200)

sum(tmpp$ecisL > 0 | tmpp$ecisH <0,na.rm = T)
sum(tmpp$ecisL > 0.2 | tmpp$ecisH < -0.2 ,na.rm = T)
sum(tmpp$ecisL > 0.5 | tmpp$ecisH < -0.5,na.rm = T)
sum(tmpp$ecisL > 1 | tmpp$ecisH < -1 ,na.rm = T)

sum(tmpp$ecisL > 0 | tmpp$ecisH <0,na.rm = T) / nrow(tmpp)
sum(tmpp$ecisL > 0.2 | tmpp$ecisH < -0.2 ,na.rm = T) /nrow(tmpp)
sum(tmpp$ecisL > 0.5 | tmpp$ecisH < -0.5,na.rm = T) /nrow(tmpp)
sum(tmpp$ecisL > 1 | tmpp$ecisH < -1,na.rm = T ) /nrow(tmpp)

