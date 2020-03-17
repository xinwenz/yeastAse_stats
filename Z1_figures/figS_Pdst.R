load("~/cloud/project/R1_hpc_null_nostr_ci/rep_null_chs500.RData")
load("~/cloud/project/hpc_pre/Gier2015/hpc_need_gier.RData")

load("~/cloud/project/M6_expdc_ci_hybrid/expdc_hybrid_ci_450.RData")
load("~/cloud/project/hpc_pre/expdc_filter/hpc_need.RData")


lamd12 <- rowSums(expdc[,-1])
allexp <- sum(colSums(expdc[,-1]))/2
lambda12 <- lamd12/allexp

theta1 <-  expdc_hybrid_ci_450$log.rHy
#ggplot() + geom_histogram(aes(x=theta))
o1 <- 2 ^ theta1 * lambda12 
#ggplot() + geom_histogram(aes(x=o1)) + xlim(c(0,2))
p1 <- 2 ^ theta1 * lambda12 / (1 + 2 ^ theta1 * lambda12)
ggplot() + geom_histogram(aes(x=p1)) + 
theme_bw() +
    labs(title="Distribution of gene-specific overdispersion in Data:Xinw2018", y="number of genes", x="p parameter in the Negative-Binomial model") + 
    theme(axis.text.x=element_text(angle=0, hjust=1)) +
    #coord_cartesian(ylim=c(0, 0.6)) +
    scale_x_continuous(breaks=seq(0,1,by=0.1)) + 
    #scale_y_continuous(breaks = seq(0,500,by=100)) + 
    #coord_cartesian(xlim=c(-4, 4)) +
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

#theta2 <- exph_yps_0500$BBf_rHy
#ggplot() + geom_histogram(aes(x=theta2))
#o2 <- 2^ theta2 * lambda12 
#ggplot() + geom_histogram(aes(x=o2)) + xlim(c(0,10))
#p2 <- 2 ^ theta2 * lambda12 / (1 + 2 ^ theta2 * lambda12)   
#ggplot() + geom_histogram(aes(x=p2)) 


lamd12_44 <- rowSums(rep44[,-1])
allexp_44 <- sum(colSums(rep44[,-1]))/2
lambda12_44 <- lamd12_44/allexp_44
theta_44 <-  rep_null_1000$log.rHy
#ggplot() + geom_histogram(aes(x=theta_44))
#o_44 <- 2 ^ theta_44 * lambda12_44 
#ggplot() + geom_histogram(aes(x=o_44)) + xlim(c(0,2))
p_44 <- 2 ^ theta_44 * lambda12_44 / (1 + 2 ^ theta_44 * lambda12_44)
ggplot() + geom_histogram(aes(x=p_44)) +
    theme_bw() +
    labs(title="Distribution of gene-specific overdispersion in Data:Gier2105", y="number of genes", x="p parameter in the Negative-Binomial model") + 
    theme(axis.text.x=element_text(angle=0, hjust=1)) +
    #coord_cartesian(ylim=c(0, 0.6)) +
    scale_x_continuous(breaks=seq(0,1,by=0.1)) + 
    #scale_y_continuous(breaks = seq(0,500,by=100)) + 
    #coord_cartesian(xlim=c(-4, 4)) +
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

