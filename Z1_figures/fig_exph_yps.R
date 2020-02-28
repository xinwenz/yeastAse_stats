library(ggplot2)
library(dplyr)
library(gridExtra)
library(reshape2)

fpo_format <- function(nm,each_nr,nr,nGenes,nbit){
    varName <- nm
    each_nr <- each_nr
    nr <- nr
    nGenes <- nGenes
    
    
    BB_matrix <- matrix(nrow=each_nr,ncol=nr)
    BI_matrix <- matrix(nrow=each_nr,ncol=nr)
    
    for(i in 1:nr) {
        for (j in 1:each_nr) {
            #num <- ifelse( i < nRep, (i-1) * eachRep + j,  (i-1)*eachRep + 1 ) 
            num <- (i-1) * each_nr + j 
            myVarName <- paste0(varName,formatC(num,width=nbit,flag="0",format="fg"))
            tmp <- get(myVarName)
            BB_matrix[j,i] <- sum(tmp$ecisL > 0 | tmp$ecisH < 0,na.rm = T)/nGenes
            BI_matrix[j,i] <- sum(tmp$B.ecisL > 0 | tmp$B.ecisH <0,na.rm = T)/nGenes
        }
    }
    
    rownames(BB_matrix) <- 1:each_nr
    BB_df <- data.frame(BB_matrix)
    #colnames(BB_df) <- paste0("Nr",1:nr)
    BB_df_m <- melt(BB_df)
    colnames(BB_df_m) <- c("Used_replicates","fp_rate")
    
    
    rownames(BI_matrix) <- 1:each_nr
    BI_df <- data.frame(BI_matrix)
    #colnames(BI_df) <- paste0("Nr",1:nr)
    BI_df_m <- melt(BI_df)
    colnames(BI_df_m) <- c("Used_replicates","fp_rate")
    
    dfmm <- rbind(cbind(BB_df_m,type="Beta-Binom"),cbind(BI_df_m,type="Binom")) 
}

for(i in 1:3000) {
    mydataName <- paste0("~/cloud/project/hpc_yps149_nostr/exph_ypsRM14A9A_esum_",formatC(i,width=4,flag="0"),".RData")
    load(mydataName)
}

df_exphyps <- fpo_format(nm="exph_ypsRM14A9A_esum_",each_nr = 150 ,nr = 20 ,nGenes = 4710,nbit = 4)

ggplot(df_exphyps, aes(x=Used_replicates, y=fp_rate, shape=type)) + 
    geom_abline(slope = 0, intercept = 0.05, color='brown') + 
    geom_violin(scale="count",adjust= 4, aes(fill = type, color = type),width=1, position =position_dodge(width = 0.8) , trim=T)  + 
    geom_errorbar(stat="summary", fun.ymin=quantile,fun.ymax= quantile,fun.args = list(0.5), lwd=1,position =position_dodge(0.8))  +
    geom_errorbar(stat="summary", fun.ymin=quantile,fun.ymax= quantile,fun.args = list(0.25), lwd=0.5,linetype="dotted",position =position_dodge(0.8)) +
    geom_errorbar(stat="summary", fun.ymin=quantile,fun.ymax= quantile,fun.args = list(0.75), lwd=0.5,linetype="dotted",position =position_dodge(0.8)) +
    scale_fill_manual(values=c("steelblue2","tan")) +
    scale_color_manual(values=c( "steelblue2","tan")) +
  labs(title="Xinw2018_yps: false posistive rate with different levels of replication", y="false positive rate", x="level of replication (Nr)") +
    theme_bw() +
    theme(axis.text.x=element_text(angle=0, hjust=1)) +
    coord_cartesian(ylim=c(0, 0.6)) +
    scale_x_discrete(labels=1:35) + 
    scale_y_continuous(breaks = seq(0,0.6,by=0.05)) + 
    theme(
        panel.grid = element_line(color="grey85"),  
        legend.position = c(0.875, 0.9),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=14,face="plain"),
        plot.title = element_text(size = 18, face = "bold",hjust=.5),
        axis.text.x = element_text(colour="grey20",size=13,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=13,hjust=.5,vjust=.5,face="plain"),
        axis.title.x = element_text(colour="grey20",size=15,hjust=.5,vjust=0,face="bold"), 
        axis.title.y = element_text(colour="grey20",size=15,hjust=.5,vjust=.5,face="bold"))

ggsave("fig_exph_yps.pdf", device = "pdf", width = 25, height = 15, units = "cm",dpi=300)
