library(ggplot2)
library(dplyr)
library(gridExtra)
library(reshape2)

# function gives out BI_median and BB_median to global environment. 
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
            BB_matrix[j,i] <- sum(tmp$pval1 < 0.05 ,na.rm = T)/nGenes
            BI_matrix[j,i] <- sum(tmp$pval3 < 0.05 ,na.rm = T)/nGenes
        }
    }
    
    rownames(BB_matrix) <- 1:each_nr
    BB_matrix <- cbind(BB_matrix, matrix(NA,nrow=each_nr, ncol=2)) 
    BB_df <- data.frame(BB_matrix)
    BB_median <<- apply(BB_df, 2,median)
    #colnames(BB_df) <- paste0("Nr",1:nr)
    BB_df_m <- melt(BB_df)
    colnames(BB_df_m) <- c("Used_replicates","fp_rate")
    
    
    rownames(BI_matrix) <- 1:each_nr
    BI_matrix <- cbind(BI_matrix, matrix(NA,nrow=each_nr, ncol=2))
    BI_df <- data.frame(BI_matrix)
    BI_median <<- apply(BI_df, 2,median)
    #colnames(BI_df) <- paste0("Nr",1:nr)
    BI_df_m <- melt(BI_df)
    colnames(BI_df_m) <- c("Used_replicates","fp_rate")
    
    
    dfmm <- rbind(cbind(BB_df_m,type="Beta-Binom"),cbind(BI_df_m,type="Binom")) 
}

for(i in 1:360) {
    mydataName <- paste0("~/cloud/project/M4_expdc_liko/expdc_liko_",formatC(i,width=3,flag="0"),".RData")
    load(mydataName)
}

expdc_liko_df <- fpo_format(nm="expdc_liko_",each_nr = 20,nr = 18, nGenes= 4631,nbit=3)


#dt_BI <- cbind(1:19,BI_median[1:19])
#dt_BI <- data.frame(dt_BI)
#colnames(dt_BI) <- c("x","y")
#fit_BI <-  nls(y ~ SSasymp(x, Asym, R0, lrc), data = dt_BI)

#fy1i <- coef(fit_BI)[1]
#fy2i <- coef(fit_BI)[2]
#fy3i <- coef(fit_BI)[3]

#dt_BB <- cbind(1:20,BB_median[1:20])
#dt_BB <- data.frame(dt_BB)
#colnames(dt_BB) <- c("x","y")
#fit_BB <-  nls(y ~ SSasymp(x, Asym, R0, lrc), data = dt_BB)

#fy1b <- coef(fit_BB)[1]
#fy2b <- coef(fit_BB)[2]
#fy3b <- coef(fit_BB)[3]


ggplot(expdc_liko_df, aes(x=Used_replicates, y=fp_rate, shape=type)) + 
    #geom_abline(slope = 0, intercept = fy1i, color='black',linetype="dashed") +
    #geom_abline(slope = 0, intercept = fy1b, color='black',linetype="dashed") +
    #stat_function(fun=function(x) fy1i + (fy2i-fy1i)*exp(-exp(fy3i) * x ) , geom='line', color="orange", linetype="solid" ) + 
    #stat_function(fun=function(x) fy1b + (fy2b-fy1b)*exp(-exp(fy3b) * x ) , geom='line', color="cyan", linetype="solid" ) + 
    geom_violin(scale="count",adjust= 1, aes(fill = type, color = type),width=2, position =position_dodge(width = 0.2) , trim=T)  + 
    geom_errorbar(stat="summary", fun.ymin=quantile,fun.ymax= quantile,fun.args = list(0.5), lwd=0.5,position =position_dodge(0.2), aes(color=type))  +
    geom_errorbar(stat="summary", fun.ymin=quantile,fun.ymax= quantile,fun.args = list(0.25), size=0.25,linetype="dotted",position =position_dodge(0.2), color="brown") +
    geom_errorbar(stat="summary", fun.ymin=quantile,fun.ymax= quantile,fun.args = list(0.75), size=0.25,linetype="dotted",position =position_dodge(0.2), color="brown") +
    scale_fill_manual(values=c("steelblue2","tan")) +
    scale_color_manual(values=c( "steelblue2","tan")) +
    labs(title="M_4:epxdc_liko_0.05", y="rejection rate", x="level of replication (Nr)") + 
    theme_bw() +
    theme(axis.text.x=element_text(angle=0, hjust=1)) +
    coord_cartesian(xlim=c(1,20), ylim=c(0, 0.8)) +
    scale_x_discrete(labels=1:20) +
    #scale_x_log10() + s
    scale_y_continuous(breaks = seq(0,0.8,by=0.05)) + 
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

#ggsave("fig_exph.pdf", device = "pdf", width = 20, height = 15, units = "cm",dpi=300)

#BBBIdt <- rbind(cbind(dt_BI,type="Binom"), cbind(dt_BB,type="Beta-Binom"))
#ggplot(BBBIdt, aes(x=x, y =y, color=type)) +geom_point() + geom_path() + scale_x_log10() + labs(x="log10(x)")

ggplot(expdc_liko_360, aes(x=BBf_ec)) + 
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


#######################################
tmpp <- expdc_liko_360

godlist0 <- which(tmpp$pval1 < 0.05)
#godlist0N <- (1:4710)[-godlist0]

sigdf <- tmpp[godlist0,]
sort_sigdf <- sort(abs(sigdf$BBf_ec),decreasing = T)

mydf <- data.frame(prop = 1:length(sort_sigdf) / 4631, absecis = sort_sigdf )

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




