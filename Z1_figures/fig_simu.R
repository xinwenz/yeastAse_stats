library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)
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

for(i in 1:3750) {
    mydataName <- paste0("~/cloud/project/hpc_simu_p1/simu_p1_",formatC(i,width=4,flag="0"),".RData")
    load(mydataName)
}

for(i in 1:3750) {
    mydataName <- paste0("~/cloud/project/hpc_simu_p4/simu_p4_",formatC(i,width=4,flag="0"),".RData")
    load(mydataName)
}

for(i in 1:3750) {
    mydataName <- paste0("~/cloud/project/hpc_simu_p8/simu_p8_",formatC(i,width=4,flag="0"),".RData")
    load(mydataName)
}

for(i in 1:3750) {
    mydataName <- paste0("~/cloud/project/hpc_simu_pc/simu_pc_",formatC(i,width=4,flag="0"),".RData")
    load(mydataName)
}


df_simu1 <- fpo_format(nm="simu_p1_",each_nr = 150 ,nr = 25 ,nGenes = 5000,nbit=4)
df_simu4 <- fpo_format(nm="simu_p4_",each_nr = 150 ,nr = 25 ,nGenes = 5000,nbit=4)
df_simu8 <- fpo_format(nm="simu_p8_",each_nr = 150 ,nr = 25 ,nGenes = 5000,nbit=4)
df_simuc <- fpo_format(nm="simu_pc_",each_nr = 150 ,nr = 25 ,nGenes = 5000,nbit=4)

xfont <- 11 
legfont <- 11
ang <- 40
legp <- 0.91

p1 <- ggplot(df_simu1, aes(x=Used_replicates, y=fp_rate, shape=type)) + 
    geom_abline(slope = 0, intercept = 0.05, color='brown') + 
    geom_violin(scale="count",adjust= 4, aes(fill = type, color = type),width=1, position =position_dodge(width = 0.8) , trim=T)  + 
    geom_errorbar(stat="summary", fun.ymin=quantile,fun.ymax= quantile,fun.args = list(0.5), lwd=0.5,position =position_dodge(0.8))  +
    geom_errorbar(stat="summary", fun.ymin=quantile,fun.ymax= quantile,fun.args = list(0.25), size=0.3,linetype="dotted",position =position_dodge(0.8)) +
    geom_errorbar(stat="summary", fun.ymin=quantile,fun.ymax= quantile,fun.args = list(0.75), size=0.3,linetype="dotted",position =position_dodge(0.8)) +
    scale_fill_manual(values=c("steelblue2","tan")) +
    scale_color_manual(values=c( "steelblue2","tan")) +
    #labs(title="simu_null(p=0.1)", y="false positive rate", x="number of replicates tested:Nr") +
    theme_bw() +
    theme(axis.text.x=element_text(angle=ang, hjust=1)) +
    coord_cartesian(ylim=c(0, 0.5)) +
    scale_x_discrete(labels=1:25) + 
    scale_y_continuous(breaks = seq(0,0.5,by=0.05)) +
   annotate("text", x = 13, y = 0.49, label = "p=0.1", size = 6, fontface="bold", color="cyan4") +
    theme(
        panel.grid = element_line(color="grey85"),  
        legend.position = c(0.875, legp),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=legfont,face="plain"),
        #plot.title = element_text(size = 13, face = "plain",hjust=0.5, vjust=-1),
        plot.title = element_blank(),
        axis.text.x = element_text(colour="grey20",size=xfont,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=13,hjust=.5,vjust=.5,face="plain"),
        #axis.title.x = element_text(colour="grey20",size=15,hjust=.5,vjust=0,face="bold"),
        axis.title.x = element_blank(),
        #axis.title.y = element_text(colour="grey20",size=15,hjust=.5,vjust=.5,face="bold"))
        axis.title.y = element_blank())


p4 <- ggplot(df_simu4, aes(x=Used_replicates, y=fp_rate, shape=type)) + 
    geom_abline(slope = 0, intercept = 0.05, color='brown') + 
    geom_violin(scale="count",adjust= 4, aes(fill = type, color = type),width=1, position =position_dodge(width = 0.8) , trim=T)  + 
    geom_errorbar(stat="summary", fun.ymin=quantile,fun.ymax= quantile,fun.args = list(0.5), lwd=0.5,position =position_dodge(0.8))  +
    geom_errorbar(stat="summary", fun.ymin=quantile,fun.ymax= quantile,fun.args = list(0.25), size=0.3,linetype="dotted",position =position_dodge(0.8)) +
    geom_errorbar(stat="summary", fun.ymin=quantile,fun.ymax= quantile,fun.args = list(0.75), size=0.3,linetype="dotted",position =position_dodge(0.8)) +
    scale_fill_manual(values=c("steelblue2","tan")) +
    scale_color_manual(values=c( "steelblue2","tan")) +
    #labs(title="simu_null(p=0.1)", y="false positive rate", x="number of replicates tested:Nr") +
    theme_bw() +
    theme(axis.text.x=element_text(angle=ang, hjust=1)) +
    coord_cartesian(ylim=c(0, 0.5)) +
    scale_x_discrete(labels=1:25) + 
    scale_y_continuous(breaks = seq(0,0.5,by=0.05)) +
    annotate("text", x = 13, y = 0.49, label = "p=0.4", size = 6, fontface="bold", color="cyan4") +
    theme(
        panel.grid = element_line(color="grey85"),  
        legend.position = c(0.875, legp),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=legfont,face="plain"),
        #plot.title = element_text(size = 13, face = "plain",hjust=0.5, vjust=-1),
        plot.title = element_blank(),
        axis.text.x = element_text(colour="grey20",size=xfont,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=13,hjust=.5,vjust=.5,face="plain"),
        #axis.title.x = element_text(colour="grey20",size=15,hjust=.5,vjust=0,face="bold"),
        axis.title.x = element_blank(),
        #axis.title.y = element_text(colour="grey20",size=15,hjust=.5,vjust=.5,face="bold"))
        axis.title.y = element_blank())

p8 <- ggplot(df_simu8, aes(x=Used_replicates, y=fp_rate, shape=type)) + 
    geom_abline(slope = 0, intercept = 0.05, color='brown') + 
    geom_violin(scale="count",adjust= 4, aes(fill = type, color = type),width=1, position =position_dodge(width = 0.8) , trim=T)  + 
    geom_errorbar(stat="summary", fun.ymin=quantile,fun.ymax= quantile,fun.args = list(0.5), lwd=0.5,position =position_dodge(0.8))  +
    geom_errorbar(stat="summary", fun.ymin=quantile,fun.ymax= quantile,fun.args = list(0.25), size=0.3,linetype="dotted",position =position_dodge(0.8)) +
    geom_errorbar(stat="summary", fun.ymin=quantile,fun.ymax= quantile,fun.args = list(0.75), size=0.3,linetype="dotted",position =position_dodge(0.8)) +
    scale_fill_manual(values=c("steelblue2","tan")) +
    scale_color_manual(values=c( "steelblue2","tan")) +
    #labs(title="simu_null(p=0.1)", y="false positive rate", x="number of replicates tested:Nr") +
    theme_bw() +
    theme(axis.text.x=element_text(angle=ang, hjust=1)) +
    coord_cartesian(ylim=c(0, 0.5)) +
    scale_x_discrete(labels=1:25) + 
    scale_y_continuous(breaks = seq(0,0.5,by=0.05)) +
    annotate("text", x = 13, y = 0.49, label = "p=0.8", size = 6, fontface="bold", color="cyan4") +
    theme(
        panel.grid = element_line(color="grey85"),  
        legend.position = c(0.875, legp),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=legfont,face="plain"),
        #plot.title = element_text(size = 13, face = "plain",hjust=0.5, vjust=-1),
        plot.title = element_blank(),
        axis.text.x = element_text(colour="grey20",size=xfont,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=13,hjust=.5,vjust=.5,face="plain"),
        #axis.title.x = element_text(colour="grey20",size=15,hjust=.5,vjust=0,face="bold"),
        axis.title.x = element_blank(),
        #axis.title.y = element_text(colour="grey20",size=15,hjust=.5,vjust=.5,face="bold"))
        axis.title.y = element_blank())

pc <- ggplot(df_simuc, aes(x=Used_replicates, y=fp_rate, shape=type)) + 
    geom_abline(slope = 0, intercept = 0.05, color='brown') + 
    geom_violin(scale="count",adjust= 4, aes(fill = type, color = type),width=1, position =position_dodge(width = 0.8) , trim=T)  + 
    geom_errorbar(stat="summary", fun.ymin=quantile,fun.ymax= quantile,fun.args = list(0.5), lwd=0.5,position =position_dodge(0.8))  +
    geom_errorbar(stat="summary", fun.ymin=quantile,fun.ymax= quantile,fun.args = list(0.25), size=0.3,linetype="dotted",position =position_dodge(0.8)) +
    geom_errorbar(stat="summary", fun.ymin=quantile,fun.ymax= quantile,fun.args = list(0.75), size=0.3,linetype="dotted",position =position_dodge(0.8)) +
    scale_fill_manual(values=c("steelblue2","tan")) +
    scale_color_manual(values=c( "steelblue2","tan")) +
    #labs(title="simu_null(p=0.1)", y="false positive rate", x="number of replicates tested:Nr") +
    theme_bw() +
    theme(axis.text.x=element_text(angle=ang, hjust=1)) +
    coord_cartesian(ylim=c(0, 0.5)) +
    scale_x_discrete(labels=1:25) + 
    scale_y_continuous(breaks = seq(0,0.5,by=0.05)) +
    annotate("text", x = 13, y = 0.49, label = "p=unif(0,0.8)", size = 6, fontface="bold", color="cyan4") +
    theme(
        panel.grid = element_line(color="grey85"),  
        legend.position = c(0.875, legp),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=legfont,face="plain"),
        #plot.title = element_text(size = 13, face = "plain",hjust=0.5, vjust=-1),
        plot.title = element_blank(),
        axis.text.x = element_text(colour="grey20",size=xfont,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=13,hjust=.5,vjust=.5,face="plain"),
        #axis.title.x = element_text(colour="grey20",size=15,hjust=.5,vjust=0,face="bold"),
        axis.title.x = element_blank(),
        #axis.title.y = element_text(colour="grey20",size=15,hjust=.5,vjust=.5,face="bold"))
        axis.title.y = element_blank())


tx1=textGrob("simulated null datasets: false postive rate with different levels of replication", gp=gpar(fontface="bold", fontsize = 18))
tx2=textGrob("false positive rate", gp=gpar(fontface="bold", fontsize=15), rot=90)
tx3=textGrob("level of replication (Nr)", gp=gpar(fontface="bold", fontsize=15))

pa <- grid.arrange(p1,p4,p8,pc,ncol=2,top=tx1,bottom=tx3, left=tx2)
ggsave(file="fig_simu.pdf",pa, device = "pdf", width = 30, height = 20, units = "cm",dpi=200)
#pdf(pa, file="test.pdf")



