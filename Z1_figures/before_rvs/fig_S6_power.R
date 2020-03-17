library(ggplot2)
library(dplyr)
library(gridExtra)
library(reshape2)

for(i in 1:2700) {
    mydataName <- paste0("~/cloud/project/hpc_expTrueRM14A9A//exph_trueRM14A9A_esum_",formatC(i,width=4,flag="0"),".RData")
    load(mydataName)
}

tmpp <- exph_trueRM14A9A_esum_2700
godlist0 <- which(tmpp$ecisL > 0 | tmpp$ecisH <0)
godlist0N <- (1:4710)[-godlist0]

cal_power <- function(nm,eachRep,nRep,repGenes){
    varName <- nm
    eachRep <- eachRep
    nRep <- nRep
    repGenes <- repGenes
    
    mygodlist <- get(paste0("godlist",0))
    mygodlistN <- get(paste0('godlist',0,'N'))
    BB_matrix <- matrix(nrow=eachRep,ncol=nRep)
    BI_matrix <- matrix(nrow=eachRep,ncol=nRep)
    
    
    for(i in 1:nRep) {
        for (j in 1:eachRep) {
            #num <- ifelse( i < nRep, (i-1) * eachRep + j,  (i-1)*eachRep + 1 ) 
            num <- (i-1) * eachRep + j 
            myVarName <- paste0(varName,formatC(num,width=4,flag="0",format="fg"))
            tmp <- get(myVarName)
            
            
            sigBB_list <- which(tmp$ecisL > 0 | tmp$ecisH < 0)
            sigBI_list <- which(tmp$B.ecisL > 0 | tmp$B.ecisH < 0)
            
            BB_matrix[j,i] <- length(intersect(sigBB_list, mygodlist)) /length(mygodlist)
            BI_matrix[j,i] <- length(intersect(sigBI_list, mygodlist)) /length(mygodlist)
        }
 }
    
    rownames(BB_matrix) <- 1:eachRep
    BB_df <- data.frame(BB_matrix)
    colnames(BB_df) <- paste0("UsRp",1:nRep)
    BB_df_m <- melt(BB_df)
    colnames(BB_df_m) <- c("Used_replicates","fp_rate")
    
    
    rownames(BI_matrix) <- 1:eachRep
    BI_df <- data.frame(BI_matrix)
    colnames(BI_df) <- paste0("UsRp",1:nRep)
    BI_df_m <- melt(BI_df)
    colnames(BI_df_m) <- c("Used_replicates","fp_rate")
    
    dfmm <- rbind(cbind(BB_df_m,type="Beta-Binom"),cbind(BI_df_m,type="Binom")) 
}

ans <- cal_power(nm ="exph_trueRM14A9A_esum_",eachRep = 150,nRep = 18,repGenes = 4710 ) 

ggplot(ans, aes(x=Used_replicates, y=fp_rate, shape=type)) + 
    #geom_abline(slope = 0, intercept = 0.05, color='brown') + 
    geom_violin(scale="count",adjust= 4, aes(fill = type, color = type),width=1, position =position_dodge(width = 0.3) , trim=T)  + 
    geom_errorbar(stat="summary", fun.ymin=quantile,fun.ymax= quantile,fun.args = list(0.5), lwd=0.5,position =position_dodge(0.3))  +
    geom_errorbar(stat="summary", fun.ymin=quantile,fun.ymax= quantile,fun.args = list(0.25), lwd=0.3,linetype="dotted",position =position_dodge(0.3)) +
    geom_errorbar(stat="summary", fun.ymin=quantile,fun.ymax= quantile,fun.args = list(0.75), lwd=0.3,linetype="dotted",position =position_dodge(0.3)) +
    scale_fill_manual(values=c("steelblue2","tan")) +
    scale_color_manual(values=c( "steelblue2","tan")) +
    labs(title="Xinw2018: statistical power with different levels of replication", y="power", x="level of replication (Nr)") +
    theme_bw() +
    theme(axis.text.x=element_text(angle=0, hjust=1)) +
    coord_cartesian(ylim=c(0, 1)) +
    scale_x_discrete(labels=1:18) + 
    scale_y_continuous(breaks = seq(0,1,by=0.05)) + 
    theme(
        panel.grid = element_line(color="grey85"),  
        legend.position = c(0.875, 0.5),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=14,face="plain"),
        plot.title = element_text(size = 18, face = "bold",hjust=.5),
        axis.text.x = element_text(colour="grey20",size=13,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=13,hjust=.5,vjust=.5,face="plain"),
        axis.title.x = element_text(colour="grey20",size=15,hjust=.5,vjust=0,face="bold"), 
        axis.title.y = element_text(colour="grey20",size=15,hjust=.5,vjust=.5,face="bold"))

ggsave("~/cloud/project/figures/fig_S6_power.pdf", device = "pdf", width = 25, height = 15, units = "cm",dpi=300)
