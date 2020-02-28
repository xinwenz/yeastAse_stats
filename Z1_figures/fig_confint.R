load("~/cloud/project/hpc_expTrueRM14A9A/exph_trueRM14A9A_esum_0450.RData")
load("~/cloud/project/hpc_expTrueRM14A9A/exph_trueRM14A9A_esum_0900.RData")
load("~/cloud/project/hpc_expTrueRM14A9A/exph_trueRM14A9A_esum_1800.RData")
load("~/cloud/project/hpc_expTrueRM14A9A/exph_trueRM14A9A_esum_2700.RData")

library(ggplot2)
library(dplyr)
library(gridExtra)
library(reshape2)
library(grid)
library(gridBase)

    siz = 100
    smp_within <- sort(sample(x = 1:250,size = siz,replace = F))
    
    smp_low <- 700:1000
    smp_mid <- 2700:3000
    smp_high <- 4200:4500
    
    
    ## repilicate Number : 3 -- (3 * 150) 
    gfl3 <- exph_trueRM14A9A_esum_0450
    gfl6 <- exph_trueRM14A9A_esum_0900
    gfl12 <- exph_trueRM14A9A_esum_1800
    gfl18 <- exph_trueRM14A9A_esum_2700
    
    g_all <- cbind(gfl3,gfl6,gfl12,gfl18)
    
    colnames(g_all) <- paste0("nm_",rep(1:4,each=9),"_",rep(1:9,4))
    tmp_low_t <- g_all[smp_low,] %>% mutate(avg = (nm_1_1 + nm_2_1 + nm_3_1 + nm_4_1)/4) %>% arrange(avg) %>% filter(
      is.na(nm_1_2) == F & is.na(nm_2_2) ==F & is.na(nm_3_2) ==F & is.na(nm_4_2) ==F  ) 
    tmp_low <-  tmp_low_t %>% slice(sort(sample(x = 1:nrow(tmp_low_t),size = siz,replace = F)))
    
    tmp_mid_t <- g_all[smp_mid,] %>% mutate(avg = (nm_1_1 + nm_2_1 + nm_3_1 + nm_4_1)/4) %>% arrange(avg) %>% filter(
      is.na(nm_1_2) == F & is.na(nm_2_2) ==F & is.na(nm_3_2) ==F & is.na(nm_4_2) ==F   ) 
    tmp_mid <-  tmp_mid_t %>% slice(sort(sample(x = 1:nrow(tmp_mid_t),size = siz,replace = F)))
    
    tmp_high_t <- g_all[smp_high,] %>% mutate(avg = (nm_1_1 + nm_2_1 + nm_3_1 + nm_4_1)/4) %>% arrange(avg) %>% filter(
      is.na(nm_1_2) == F & is.na(nm_2_2) ==F & is.na(nm_3_2) ==F & is.na(nm_4_2) ==F ) 
    tmp_high <-  tmp_high_t %>% slice(sort(sample(x = 1:nrow(tmp_high_t),size = siz,replace = F)))
     
    tmp1 <- cbind(1:siz,tmp_low[,1:3], "Beta-Binom",   'low-expression-genes','replicates:3')
    tmp2 <- cbind(1:siz,tmp_low[,7:9], "Binom",        "low-expression-genes", 'replicates:3')
    tmp3 <- cbind(1:siz,tmp_low[,10:12], "Beta-Binom", 'low-expression-genes','replicates:6')
    tmp4 <- cbind(1:siz,tmp_low[,16:18], "Binom",      "low-expression-genes", 'replicates:6')   
    tmp5 <- cbind(1:siz,tmp_low[,19:21], "Beta-Binom", 'low-expression-genes','replicates:12')
    tmp6 <- cbind(1:siz,tmp_low[,25:27], "Binom",      "low-expression-genes", 'replicates:12')
    tmp7 <- cbind(1:siz,tmp_low[,28:30], "Beta-Binom", 'low-expression-genes','replicates:18')
    tmp8 <- cbind(1:siz,tmp_low[,34:36], "Binom",      "low-expression-genes", 'replicates:18')   
    
    tmp9 <-  cbind(1:siz,tmp_mid[,1:3], "Beta-Binom",   'mid-expression-genes','replicates:3')
    tmp10 <- cbind(1:siz,tmp_mid[,7:9], "Binom",        "mid-expression-genes", 'replicates:3')
    tmp11 <- cbind(1:siz,tmp_mid[,10:12], "Beta-Binom", 'mid-expression-genes','replicates:6')
    tmp12 <- cbind(1:siz,tmp_mid[,16:18], "Binom",      "mid-expression-genes", 'replicates:6')   
    tmp13 <- cbind(1:siz,tmp_mid[,19:21], "Beta-Binom", 'mid-expression-genes','replicates:12')
    tmp14 <- cbind(1:siz,tmp_mid[,25:27], "Binom",      "mid-expression-genes", 'replicates:12')
    tmp15 <- cbind(1:siz,tmp_mid[,28:30], "Beta-Binom", 'mid-expression-genes','replicates:18')
    tmp16 <- cbind(1:siz,tmp_mid[,34:36], "Binom",      "mid-expression-genes", 'replicates:18')   
    
    
    tmp17 <- cbind(1:siz,tmp_high[,1:3], "Beta-Binom",   'high-expression-genes','replicates:3')
    tmp18 <- cbind(1:siz,tmp_high[,7:9], "Binom",        "high-expression-genes", 'replicates:3')
    tmp19 <- cbind(1:siz,tmp_high[,10:12], "Beta-Binom", 'high-expression-genes','replicates:6')
    tmp20 <- cbind(1:siz,tmp_high[,16:18], "Binom",      "high-expression-genes", 'replicates:6')   
    tmp21 <- cbind(1:siz,tmp_high[,19:21], "Beta-Binom", 'high-expression-genes','replicates:12')
    tmp22 <- cbind(1:siz,tmp_high[,25:27], "Binom",      "high-expression-genes", 'replicates:12')
    tmp23 <- cbind(1:siz,tmp_high[,28:30], "Beta-Binom", 'high-expression-genes','replicates:18')
    tmp24 <- cbind(1:siz,tmp_high[,34:36], "Binom",      "high-expression-genes", 'replicates:18')  
    
    bigdataframe <- c()
    for(i in 1:24) {
      ele <- get(paste0("tmp", i))
      colnames(ele) <- c("geneRank","log2ecis","ecisL","ecisH","type","expression","replicates_used")
      bigdataframe <- rbind(bigdataframe, ele)
    }
    
    
#    tx1=textGrob("simulated null datasets: false postive rate with different number of replicates", gp=gpar(fontface="bold", fontsize = 18))
#    tx2=textGrob("false positive rate", gp=gpar(fontface="bold", fontsize=15), rot=90)
#    tx3=textGrob("number of replicates tested: Nr", gp=gpar(fontface="bold", fontsize=15))
  
    ggplot(bigdataframe,aes(x=geneRank)) + 
        geom_linerange(aes(ymin=ecisL,ymax=ecisH, color=type, alpha=type), lwd=1.5)+
        scale_color_manual(values=c("cyan","black")) +
        scale_alpha_manual(values=c(1,0.5)) + 
        geom_abline(slope = 0,intercept = 0) +
        #labs(title = paste(gfln,'replicats ,betabinom vs binom'),x='gene rank') + 
        scale_y_continuous(breaks = seq(-0.75, 0.75, by= 0.25)) +
        coord_cartesian(ylim=c(-0.8, 0.8)) + 
      theme_bw() + 
      labs(title="Xinw2018: confidence interval with \n different levels of replication and expression level ", y=expression(paste("Confidence Interval: ", "log"["2"] ,"(" ,"e" ["cis"] ,")" )), 
           x=expression(paste("Gene Ranked by " ,"log"["2"], "(","e"["cis"],")")) ) +
      #xlab(bquote('   '~CO[2]~ m^-2~' '  )) +
#      labs(x=bquote('   '~CO[2]*m^2~' '  )) +
      theme(
        panel.grid = element_line(color="grey85"),  
        legend.position = c(0.82, 0.95),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=16,face="bold"),
        plot.title = element_text(size = 20, face = "bold",hjust=0.5, vjust=-1),
        #plot.title = element_blank(),
        axis.text.x = element_text(colour="grey20",size=13,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=13,hjust=.5,vjust=.5,face="plain"),
        axis.title.x = element_text(colour="grey20",size=18,hjust=.5,vjust=0,face="bold"),
        #axis.title.x = element_blank(),
        axis.title.y = element_text(colour="grey20",size=18,hjust=.5,vjust=.5,face="bold")) +
        #axis.title.y = element_blank()) + 
       facet_grid(rows = vars(expression),cols=vars(replicates_used), drop=F) +
       theme(strip.text = element_text(size = 15, colour = "midnightblue",face="bold"))
      
ggsave("fig_confint.pdf", device = "pdf", width = 40, height = 25, units = "cm",dpi=300)
  
rowMeans(exph[,-1])[smp_low]
rowMeans(exph[,-1])[smp_mid]
rowMeans(exph[,-1])[smp_high]
