library(ggplot2)
library(dplyr)
library(gridExtra)
library(reshape2)


setwd(dir = "~/cloud/project/rep_null")
for(i in 1:4400) {
    mydataName <- paste0("rep_null_",formatC(i,width=4,flag="0"),".RData")
    load(mydataName)
}

## 4 plots choase one as a repressnetative for that replicate Number: 2,8,16,22 
plot_rep_null <- function() {
    siz = 150
    smp <- sort(sample(x = 1:nrow(rep_null_0001),size = siz,replace = F))
    
    
    ## repilicate Number : 2 -- (2 * 150) 
    gfl <- rep_null_0300
    gfln <- "2"
    tmp <- gfl[smp,] %>% arrange(`log.ecis`)
    tmp1 <- tmp[,1:3]
    tmp2 <- tmp[,7:9]
    names(tmp2) <- names(tmp1)
    null_two <- rbind(cbind(tmp1,type="Betabinom",od=1:siz),cbind(tmp2,type='binom',od=1:siz)) 
    plot01 <- ggplot(null_two,aes(x=od)) + 
        geom_ribbon(aes(ymin=ecisL,ymax=ecisH,fill=type),alpha=0.3)+
        scale_fill_manual(values=c( "steelblue","yellow"),guide=FALSE) +
        geom_abline(slope = 0,intercept = 0) +
        labs(title = paste(gfln,'replicats ,betabinom vs binom')) + 
        scale_y_continuous(breaks = seq(-0.75, 0.75, by= 0.25), limits=c(-0.8,0.8))
    
    
    ## replicate Numnber 8 
    gfl <- rep_null_0750
    gfln <- "5" 
    tmp <- gfl[smp,] %>% arrange(`log.ecis`)
    tmp1 <- tmp[,1:3]
    tmp2 <- tmp[,7:9]
    names(tmp2) <- names(tmp1)
    null_two <- rbind(cbind(tmp1,type="Betabinom",od=1:siz),cbind(tmp2,type='binom',od=1:siz)) 
    plot02 <- ggplot(null_two,aes(x=od)) + 
        geom_ribbon(aes(ymin=ecisL,ymax=ecisH,fill=type),alpha=0.3)+
        scale_fill_manual(values=c( "steelblue","yellow"),guide=FALSE) +
        geom_abline(slope = 0,intercept = 0) +
        labs(title = paste(gfln,'replicats ,betabinom vs binom')) +
        scale_y_continuous(breaks = seq(-0.75, 0.75, by= 0.25), limits=c(-.8,0.8))
    
    
    gfl <- rep_null_1500
    gfln <- "10" 
    tmp <- gfl[smp,] %>% arrange(`log.ecis`)
    tmp1 <- tmp[,1:3]
    tmp2 <- tmp[,7:9]
    names(tmp2) <- names(tmp1)
    null_two <- rbind(cbind(tmp1,type="Betabinom",od=1:siz),cbind(tmp2,type='binom',od=1:siz)) 
    plot03 <- ggplot(null_two,aes(x=od)) + 
        geom_ribbon(aes(ymin=ecisL,ymax=ecisH,fill=type),alpha=0.3)+
        scale_fill_manual(values=c( "steelblue","yellow"),guide=FALSE) +
        geom_abline(slope = 0,intercept = 0) +
        labs(title = paste(gfln,'replicats ,betabinom vs binom')) +
        scale_y_continuous(breaks = seq(-0.75, 0.75, by= 0.25), limits=c(-.8,0.8))
    
    gfl <- rep_null_3000
    gfln <- "20"
    tmp <- gfl[smp,] %>% arrange(`log.ecis`)
    tmp1 <- tmp[,1:3]
    tmp2 <- tmp[,7:9]
    names(tmp2) <- names(tmp1)
    null_two <- rbind(cbind(tmp1,type="Betabinom",od=1:siz),cbind(tmp2,type='binom',od=1:siz)) 
    plot04 <- ggplot(null_two,aes(x=od)) + 
        geom_ribbon(aes(ymin=ecisL,ymax=ecisH,fill=type),alpha=0.3)+
        scale_fill_manual(values=c( "steelblue","yellow")) +
        geom_abline(slope = 0,intercept = 0) +
        labs(title = paste(gfln,'replicats ,betabinom vs binom')) +
        scale_y_continuous(breaks = seq(-0.75, 0.75, by= 0.25), limits=c(-0.8,0.8))
    
    grid.arrange(plot01,plot02,plot03,plot04,ncol=4)
}
fpo_plot <- function(){
    varName <- "rep_null_"
    eachRep <- 200
    nRep <- 22
    repGenes <- 6023
    
    
    BB_matrix <- matrix(nrow=eachRep,ncol=nRep)
    BI_matrix <- matrix(nrow=eachRep,ncol=nRep)
    for(i in 1:nRep) {
        for (j in 1:eachRep) {
            #num <- ifelse( i < nRep, (i-1) * eachRep + j,  (i-1)*eachRep + 1 ) 
            num <- (i-1) * eachRep + j 
            myVarName <- paste0(varName,formatC(num,width=4,flag="0"))
            tmp <- get(myVarName)
            BB_matrix[j,i] <- sum(tmp$ecisL > 0 | tmp$ecisH < 0,na.rm = T)/repGenes
            BI_matrix[j,i] <- sum(tmp$B.ecisL > 0 | tmp$B.ecisH <0,na.rm = T)/repGenes
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
    
    dfmm <- rbind(cbind(BB_df_m,type="Beta-B"),cbind(BI_df_m,type="Bino")) 
    
    ggplot(dfmm, aes(x=Used_replicates, y=fp_rate,color=type)) + 
        scale_color_manual(values=c( "steelblue","orange")) +
        #geom_jitter(alpha=0.3, width=0.1,size=1.5) +
        geom_point(stat="summary", fun.y="mean",size=2.5) + 
        geom_errorbar(stat="summary", fun.data="mean_se", fun.args = list(mult = 1.96), width=0,size=1) +
        labs(title=paste(varName, "False Posistive Rate vs. Number of replicates used for two models")) 
        #theme_bw() +
        #theme(axis.text.x=element_text(angle=45, hjust=1)) 
}


####### exph_yps #### null ###########  
setwd(dir = "~/cloud/project/exph_yps_esum")
for(i in 1:91) {
  mydataName <- paste0("exph_yps_esum_",formatC(i,width=3,flag="0"),".RData")
  load(mydataName)
}

## 4 plots choase one as a repressnetative for that replicate Number: 2,5,10
plot_exph_yps_esum <- function() {
  siz = 150
  smp <- sort(sample(x = 1:nrow(exph_yps_esum_001),size = siz,replace = F))
  
  ## repilicate Number : 2 -- (2 * 10) 
  gfl <- exph_yps_esum_020
  gfln <- "2"
  tmp <- gfl[smp,] %>% arrange(`log.ecis`)
  tmp1 <- tmp[,1:3]
  tmp2 <- tmp[,7:9]
  names(tmp2) <- names(tmp1)
  null_two <- rbind(cbind(tmp1,type="Betabinom",od=1:siz),cbind(tmp2,type='binom',od=1:siz)) 
  plot01 <- ggplot(null_two,aes(x=od)) + 
    geom_ribbon(aes(ymin=ecisL,ymax=ecisH,fill=type),alpha=0.3)+
    scale_fill_manual(values=c( "steelblue","yellow"),guide=FALSE) +
    geom_abline(slope = 0,intercept = 0) +
    labs(title = paste(gfln,'replicats ,betabinom vs binom')) + 
    scale_y_continuous(breaks = seq(-0.75, 0.75, by= 0.25), limits=c(-0.8,0.8))
  
  
  ## replicate Numnber 5 (5 * 10 ) 
  gfl <- exph_yps_esum_050
  gfln <- "5" 
  tmp <- gfl[smp,] %>% arrange(`log.ecis`)
  tmp1 <- tmp[,1:3]
  tmp2 <- tmp[,7:9]
  names(tmp2) <- names(tmp1)
  null_two <- rbind(cbind(tmp1,type="Betabinom",od=1:siz),cbind(tmp2,type='binom',od=1:siz)) 
  plot02 <- ggplot(null_two,aes(x=od)) + 
    geom_ribbon(aes(ymin=ecisL,ymax=ecisH,fill=type),alpha=0.3)+
    scale_fill_manual(values=c( "steelblue","yellow"),guide=FALSE) +
    geom_abline(slope = 0,intercept = 0) +
    labs(title = paste(gfln,'replicats ,betabinom vs binom')) +
    scale_y_continuous(breaks = seq(-0.75, 0.75, by= 0.25), limits=c(-.8,0.8))
  
  
  gfl <- exph_yps_esum_091
  gfln <- "10" 
  tmp <- gfl[smp,] %>% arrange(`log.ecis`)
  tmp1 <- tmp[,1:3]
  tmp2 <- tmp[,7:9]
  names(tmp2) <- names(tmp1)
  null_two <- rbind(cbind(tmp1,type="Betabinom",od=1:siz),cbind(tmp2,type='binom',od=1:siz)) 
  plot03 <- ggplot(null_two,aes(x=od)) + 
    geom_ribbon(aes(ymin=ecisL,ymax=ecisH,fill=type),alpha=0.3)+
    scale_fill_manual(values=c( "steelblue","yellow")) +
    geom_abline(slope = 0,intercept = 0) +
    labs(title = paste(gfln,'replicats ,betabinom vs binom')) +
    scale_y_continuous(breaks = seq(-0.75, 0.75, by= 0.25), limits=c(-.8,0.8))
  
  
  grid.arrange(plot01,plot02,plot03,ncol=3)
}
fpo_plot <- function(){
  varName <- "exph_yps_esum_"
  eachRep <- 10
  nRep <- 10
  repGenes <- 4710
  
  
  BB_matrix <- matrix(nrow=eachRep,ncol=nRep)
  BI_matrix <- matrix(nrow=eachRep,ncol=nRep)
  for(i in 1:nRep) {
    for (j in 1:eachRep) {
      num <- ifelse( i < nRep, (i-1) * eachRep + j,  (i-1)*eachRep + 1 ) 
      myVarName <- paste0(varName,formatC(num,width=3,flag="0"))
      tmp <- get(myVarName)
      BB_matrix[j,i] <- sum(tmp$ecisL > 0 | tmp$ecisH < 0,na.rm = T)/repGenes
      BI_matrix[j,i] <- sum(tmp$B.ecisL > 0 | tmp$B.ecisH <0,na.rm = T)/repGenes
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
  
  dfmm <- rbind(cbind(BB_df_m,type="Beta-B"),cbind(BI_df_m,type="Bino")) 
  
  ggplot(dfmm, aes(x=Used_replicates, y=fp_rate,color=type)) + 
    scale_color_manual(values=c( "steelblue","orange")) +
    geom_jitter(alpha=0.3, width=0.1,size=1.5) +
    geom_point(stat="summary", fun.y="mean",size=2.5) + 
    geom_errorbar(stat="summary", fun.data="mean_se", fun.args = list(mult = 1.96), width=0,size=1) + 
    labs(title=paste(varName, "False Posistive Rate vs. Number of replicates used for two models")) 
  #theme_bw() +
  #theme(axis.text.x=element_text(angle=45, hjust=1)) 
}




######  plot of exph_yps without 14A #### with eSum ####### 
setwd(dir = "~/cloud/project/exph_ypsRM14A_esum")
for(i in 1:900) {
  mydataName <- paste0("exph_ypsRM14A_esum_",formatC(i,width=3,flag="0"),".RData")
  load(mydataName)
}

plot_exph_ypsRM14A_esum <- function() {
  siz = 150
  smp <- sort(sample(x = 1:nrow(exph_ypsRM14A_esum_001),size = siz,replace = F))
  
  ## repilicate Number : 2 -- (2 * 9) 
  gfl <- exph_ypsRM14A_esum_018
  gfln <- "2"
  tmp <- gfl[smp,] %>% arrange(`log.ecis`)
  tmp1 <- tmp[,1:3]
  tmp2 <- tmp[,7:9]
  names(tmp2) <- names(tmp1)
  null_two <- rbind(cbind(tmp1,type="Betabinom",od=1:siz),cbind(tmp2,type='binom',od=1:siz)) 
  plot01 <- ggplot(null_two,aes(x=od)) + 
    geom_ribbon(aes(ymin=ecisL,ymax=ecisH,fill=type),alpha=0.3)+
    scale_fill_manual(values=c( "steelblue","yellow"),guide=FALSE) +
    geom_abline(slope = 0,intercept = 0) +
    labs(title = paste(gfln,'replicats ,betabinom vs binom')) + 
    scale_y_continuous(breaks = seq(-0.75, 0.75, by= 0.25), limits=c(-0.8,0.8))
  
  
  ## replicate Numnber 5 (5 * 9 ) 
  gfl <- exph_ypsRM14A_esum_045
  gfln <- "5" 
  tmp <- gfl[smp,] %>% arrange(`log.ecis`)
  tmp1 <- tmp[,1:3]
  tmp2 <- tmp[,7:9]
  names(tmp2) <- names(tmp1)
  null_two <- rbind(cbind(tmp1,type="Betabinom",od=1:siz),cbind(tmp2,type='binom',od=1:siz)) 
  plot02 <- ggplot(null_two,aes(x=od)) + 
    geom_ribbon(aes(ymin=ecisL,ymax=ecisH,fill=type),alpha=0.3)+
    scale_fill_manual(values=c( "steelblue","yellow"),guide=FALSE) +
    geom_abline(slope = 0,intercept = 0) +
    labs(title = paste(gfln,'replicats ,betabinom vs binom')) +
    scale_y_continuous(breaks = seq(-0.75, 0.75, by= 0.25), limits=c(-.8,0.8))
  
  
  gfl <- exph_yps_esum_073
  gfln <- "9" 
  tmp <- gfl[smp,] %>% arrange(`log.ecis`)
  tmp1 <- tmp[,1:3]
  tmp2 <- tmp[,7:9]
  names(tmp2) <- names(tmp1)
  null_two <- rbind(cbind(tmp1,type="Betabinom",od=1:siz),cbind(tmp2,type='binom',od=1:siz)) 
  plot03 <- ggplot(null_two,aes(x=od)) + 
    geom_ribbon(aes(ymin=ecisL,ymax=ecisH,fill=type),alpha=0.3)+
    scale_fill_manual(values=c( "steelblue","yellow")) +
    geom_abline(slope = 0,intercept = 0) +
    labs(title = paste(gfln,'replicats ,betabinom vs binom')) +
    scale_y_continuous(breaks = seq(-0.75, 0.75, by= 0.25), limits=c(-.8,0.8))
  
  
  grid.arrange(plot01,plot02,plot03,ncol=3)
}
fpo_plot <- function(){
  varName <- "exph_ypsRM14A_esum_"
  eachRep <- 9
  nRep <- 9
  repGenes <- 4710
  
  
  BB_matrix <- matrix(nrow=eachRep,ncol=nRep)
  BI_matrix <- matrix(nrow=eachRep,ncol=nRep)
  for(i in 1:nRep) {
    for (j in 1:eachRep) {
      num <- ifelse( i < nRep, (i-1) * eachRep + j,  (i-1)*eachRep + 1 ) 
      myVarName <- paste0(varName,formatC(num,width=3,flag="0"))
      tmp <- get(myVarName)
      BB_matrix[j,i] <- sum(tmp$ecisL > 0 | tmp$ecisH < 0,na.rm = T)/repGenes
      BI_matrix[j,i] <- sum(tmp$B.ecisL > 0 | tmp$B.ecisH <0,na.rm = T)/repGenes
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
  
  dfmm <- rbind(cbind(BB_df_m,type="Beta-B"),cbind(BI_df_m,type="Bino")) 
  
  ggplot(dfmm, aes(x=Used_replicates, y=fp_rate,color=type)) + 
    scale_color_manual(values=c( "steelblue","orange")) +
    geom_jitter(alpha=0.3, width=0.1,size=1.5) +
    geom_point(stat="summary", fun.y="mean",size=2.5) + 
    geom_errorbar(stat="summary", fun.data="mean_se", fun.args = list(mult = 1.96), width=0,size=1) + 
    labs(title=paste(varName, "False Posistive Rate vs. Number of replicates used for two models")) 
  #theme_bw() +
  #theme(axis.text.x=element_text(angle=45, hjust=1)) 
}

########## plot exph_yps without 14A and 9A ## eSum ####### 
setwd(dir = "~/cloud/project/exph_ypsRM14A9A_esum")
for(i in 1:900) {
  mydataName <- paste0("exph_ypsRM14A9A_esum_",formatC(i,width=3,flag="0"),".RData")
  load(mydataName)
}

plot_exph_ypsRM14A9A_esum <- function() {
  siz = 100
  smp <- sort(sample(x = 1:nrow(exph_ypsRM14A9A_esum_001),size = siz,replace = F))
  
  ## repilicate Number : 2 -- (2 * 100) 
  gfl <- exph_ypsRM14A9A_esum_200
  gfln <- "2"
  tmp <- gfl[smp,] %>% arrange(`log.ecis`)
  tmp1 <- tmp[,1:3]
  tmp2 <- tmp[,7:9]
  names(tmp2) <- names(tmp1)
  null_two <- rbind(cbind(tmp1,type="Betabinom",od=1:siz),cbind(tmp2,type='binom',od=1:siz)) 
  plot01 <- ggplot(null_two,aes(x=od)) + 
    geom_ribbon(aes(ymin=ecisL,ymax=ecisH,fill=type),alpha=0.3)+
    scale_fill_manual(values=c( "steelblue","yellow"),guide=FALSE) +
    geom_abline(slope = 0,intercept = 0) +
    labs(title = paste(gfln,'replicats ,betabinom vs binom')) + 
    scale_y_continuous(breaks = seq(-0.75, 0.75, by= 0.25), limits=c(-0.8,0.8))
  
  
  ## replicate Numnber 5 (5 * 8 ) 
  gfl <- exph_ypsRM14A9A_esum_400
  gfln <- "4" 
  tmp <- gfl[smp,] %>% arrange(`log.ecis`)
  tmp1 <- tmp[,1:3]
  tmp2 <- tmp[,7:9]
  names(tmp2) <- names(tmp1)
  null_two <- rbind(cbind(tmp1,type="Betabinom",od=1:siz),cbind(tmp2,type='binom',od=1:siz)) 
  plot02 <- ggplot(null_two,aes(x=od)) + 
    geom_ribbon(aes(ymin=ecisL,ymax=ecisH,fill=type),alpha=0.3)+
    scale_fill_manual(values=c( "steelblue","yellow"),guide=FALSE) +
    geom_abline(slope = 0,intercept = 0) +
    labs(title = paste(gfln,'replicats ,betabinom vs binom')) +
    scale_y_continuous(breaks = seq(-0.75, 0.75, by= 0.25), limits=c(-.8,0.8))
  
  
  gfl <- exph_ypsRM14A9A_esum_600
  gfln <- "6" 
  tmp <- gfl[smp,] %>% arrange(`log.ecis`)
  tmp1 <- tmp[,1:3]
  tmp2 <- tmp[,7:9]
  names(tmp2) <- names(tmp1)
  null_two <- rbind(cbind(tmp1,type="Betabinom",od=1:siz),cbind(tmp2,type='binom',od=1:siz)) 
  plot03 <- ggplot(null_two,aes(x=od)) + 
    geom_ribbon(aes(ymin=ecisL,ymax=ecisH,fill=type),alpha=0.3)+
    scale_fill_manual(values=c( "steelblue","yellow"),guide=FALSE) +
    geom_abline(slope = 0,intercept = 0) +
    labs(title = paste(gfln,'replicats ,betabinom vs binom')) +
    scale_y_continuous(breaks = seq(-0.75, 0.75, by= 0.25), limits=c(-.8,0.8))

  gfl <- exph_ypsRM14A9A_esum_900
  gfln <- "9" 
  tmp <- gfl[smp,] %>% arrange(`log.ecis`)
  tmp1 <- tmp[,1:3]
  tmp2 <- tmp[,7:9]
  names(tmp2) <- names(tmp1)
  null_two <- rbind(cbind(tmp1,type="Betabinom",od=1:siz),cbind(tmp2,type='binom',od=1:siz)) 
  plot04 <- ggplot(null_two,aes(x=od)) + 
    geom_ribbon(aes(ymin=ecisL,ymax=ecisH,fill=type),alpha=0.3)+
    scale_fill_manual(values=c( "steelblue","yellow"),guide=FALSE) +
    geom_abline(slope = 0,intercept = 0) +
    labs(title = paste(gfln,'replicats ,betabinom vs binom')) +
    scale_y_continuous(breaks = seq(-0.75, 0.75, by= 0.25), limits=c(-.8,0.8))
  
  grid.arrange(plot01,plot02,plot03,plot04,ncol=4)
}

fpo_plot <- function(){
  varName <- "exph_ypsRM14A9A_esum_"
  eachRep <- 100
  nRep <- 9
  repGenes <- 4710
  
  
  BB_matrix <- matrix(nrow=eachRep,ncol=nRep)
  BI_matrix <- matrix(nrow=eachRep,ncol=nRep)
  for(i in 1:nRep) {
    for (j in 1:eachRep) {
      #num <- ifelse( i < nRep, (i-1) * eachRep + j,  (i-1)*eachRep + 1 ) 
      num <- (i-1) * eachRep + j 
      myVarName <- paste0(varName,formatC(num,width=3,flag="0"))
      tmp <- get(myVarName)
      BB_matrix[j,i] <- sum(tmp$ecisL > 0 | tmp$ecisH < 0,na.rm = T)/repGenes
      BI_matrix[j,i] <- sum(tmp$B.ecisL > 0 | tmp$B.ecisH <0,na.rm = T)/repGenes
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
  
  dfmm <- rbind(cbind(BB_df_m,type="Beta-B"),cbind(BI_df_m,type="Bino")) 
  
  ggplot(dfmm, aes(x=Used_replicates, y=fp_rate,color=type)) + 
    scale_color_manual(values=c( "steelblue","orange")) +
    geom_jitter(alpha=0.3, width=0.1,size=1.5) +
    geom_point(stat="summary", fun.y="mean",size=2.5) + 
    geom_errorbar(stat="summary", fun.data="mean_se", fun.args = list(mult = 1.96), width=0,size=1) + 
    labs(title=paste(varName, "Significant genes fount rate vs. Number of replicates used for two models")) 
  #theme_bw() +
  #theme(axis.text.x=element_text(angle=45, hjust=1)) 
}

###############  exph_true ################
setwd(dir = "~/cloud/project/exph_true")
for(i in 1:381) {
  mydataName <- paste0("exph_true_",formatC(i,width=3,flag="0"),".RData")
  load(mydataName)
}

plot_exph_true <- function() {
  siz = 150
  smp <- sort(sample(x = 1:nrow(exph_ypsRM14A_esum_001),size = siz,replace = F))
  
  ## repilicate Number : 2 -- (2 * 20) 
  gfl <- exph_true_040
  gfln <- "2"
  tmp <- gfl[smp,] %>% arrange(`log.ecis`)
  tmp1 <- tmp[,1:3]
  tmp2 <- tmp[,7:9]
  names(tmp2) <- names(tmp1)
  null_two <- rbind(cbind(tmp1,type="Betabinom",od=1:siz),cbind(tmp2,type='binom',od=1:siz)) 
  plot01 <- ggplot(null_two,aes(x=od)) + 
    geom_ribbon(aes(ymin=ecisL,ymax=ecisH,fill=type),alpha=0.3)+
    scale_fill_manual(values=c( "steelblue","yellow"),guide=FALSE) +
    geom_abline(slope = 0,intercept = 0) +
    labs(title = paste(gfln,'replicats ,betabinom vs binom')) + 
    scale_y_continuous(breaks = seq(-0.75, 0.75, by= 0.25), limits=c(-0.8,0.8))
  
  
  ## replicate Numnber 5 (10 * 20 ) 
  gfl <- exph_true_200
  gfln <- "10" 
  tmp <- gfl[smp,] %>% arrange(`log.ecis`)
  tmp1 <- tmp[,1:3]
  tmp2 <- tmp[,7:9]
  names(tmp2) <- names(tmp1)
  null_two <- rbind(cbind(tmp1,type="Betabinom",od=1:siz),cbind(tmp2,type='binom',od=1:siz)) 
  plot02 <- ggplot(null_two,aes(x=od)) + 
    geom_ribbon(aes(ymin=ecisL,ymax=ecisH,fill=type),alpha=0.3)+
    scale_fill_manual(values=c( "steelblue","yellow"),guide=FALSE) +
    geom_abline(slope = 0,intercept = 0) +
    labs(title = paste(gfln,'replicats ,betabinom vs binom')) +
    scale_y_continuous(breaks = seq(-0.75, 0.75, by= 0.25), limits=c(-.8,0.8))
  
  
  gfl <- exph_true_381
  gfln <- "20" 
  tmp <- gfl[smp,] %>% arrange(`log.ecis`)
  tmp1 <- tmp[,1:3]
  tmp2 <- tmp[,7:9]
  names(tmp2) <- names(tmp1)
  null_two <- rbind(cbind(tmp1,type="Betabinom",od=1:siz),cbind(tmp2,type='binom',od=1:siz)) 
  plot03 <- ggplot(null_two,aes(x=od)) + 
    geom_ribbon(aes(ymin=ecisL,ymax=ecisH,fill=type),alpha=0.3)+
    scale_fill_manual(values=c( "steelblue","yellow")) +
    geom_abline(slope = 0,intercept = 0) +
    labs(title = paste(gfln,'replicats ,betabinom vs binom')) +
    scale_y_continuous(breaks = seq(-0.75, 0.75, by= 0.25), limits=c(-.8,0.8))
  
  
  grid.arrange(plot01,plot02,plot03,ncol=3)
}

fpo_plot <- function(){
  varName <- "exph_true_"
  eachRep <- 20
  nRep <- 20
  repGenes <- 4710
  
  
  BB_matrix <- matrix(nrow=eachRep,ncol=nRep)
  BI_matrix <- matrix(nrow=eachRep,ncol=nRep)
  for(i in 1:nRep) {
    for (j in 1:eachRep) {
      num <- ifelse( i < nRep, (i-1) * eachRep + j,  (i-1)*eachRep + 1 ) 
      myVarName <- paste0(varName,formatC(num,width=3,flag="0"))
      tmp <- get(myVarName)
      BB_matrix[j,i] <- sum(tmp$ecisL > 0 | tmp$ecisH < 0,na.rm = T)/repGenes
      BI_matrix[j,i] <- sum(tmp$B.ecisL > 0 | tmp$B.ecisH <0,na.rm = T)/repGenes
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
  
  dfmm <- rbind(cbind(BB_df_m,type="Beta-B"),cbind(BI_df_m,type="Bino")) 
  
  ggplot(dfmm, aes(x=Used_replicates, y=fp_rate,color=type)) + 
    scale_color_manual(values=c( "steelblue","orange")) +
    geom_jitter(alpha=0.3, width=0.1,size=1.5) +
    geom_point(stat="summary", fun.y="mean",size=2.5) + 
    geom_errorbar(stat="summary", fun.data="mean_se", fun.args = list(mult = 1.96), width=0,size=1) +
    scale_y_continuous(breaks=seq(0.1,0.8,0.05)) + 
    labs(title=paste(varName, "Significant genes found  Rate vs. Number of replicates used for two models")) 
  #theme_bw() +
  #theme(axis.text.x=element_text(angle=45, hjust=1)) 
}


###############  exph_true without 14A and 9A ################
setwd(dir = "~/cloud/project/exph_trueRM14A9A")
for(i in 1:307) {
  mydataName <- paste0("exph_trueRM14A9A_",formatC(i,width=3,flag="0"),".RData")
  load(mydataName)
}

plot_exph_trueRM14A9A <- function() {
  siz = 150
  smp <- sort(sample(x = 1:nrow(exph_ypsRM14A_esum_001),size = siz,replace = F))
  
  ## repilicate Number : 2 -- (2 * 18) 
  gfl <- exph_trueRM14A9A_036
  gfln <- "2"
  tmp <- gfl[smp,] %>% arrange(`log.ecis`)
  tmp1 <- tmp[,1:3]
  tmp2 <- tmp[,7:9]
  names(tmp2) <- names(tmp1)
  null_two <- rbind(cbind(tmp1,type="Betabinom",od=1:siz),cbind(tmp2,type='binom',od=1:siz)) 
  plot01 <- ggplot(null_two,aes(x=od)) + 
    geom_ribbon(aes(ymin=ecisL,ymax=ecisH,fill=type),alpha=0.3)+
    scale_fill_manual(values=c( "steelblue","yellow"),guide=FALSE) +
    geom_abline(slope = 0,intercept = 0) +
    labs(title = paste(gfln,'replicats ,betabinom vs binom')) + 
    scale_y_continuous(breaks = seq(-0.75, 0.75, by= 0.25), limits=c(-0.8,0.8))
  
  
  ## replicate Numnber 5 (10 * 18 ) 
  gfl <- exph_trueRM14A9A_180
  gfln <- "10" 
  tmp <- gfl[smp,] %>% arrange(`log.ecis`)
  tmp1 <- tmp[,1:3]
  tmp2 <- tmp[,7:9]
  names(tmp2) <- names(tmp1)
  null_two <- rbind(cbind(tmp1,type="Betabinom",od=1:siz),cbind(tmp2,type='binom',od=1:siz)) 
  plot02 <- ggplot(null_two,aes(x=od)) + 
    geom_ribbon(aes(ymin=ecisL,ymax=ecisH,fill=type),alpha=0.3)+
    scale_fill_manual(values=c( "steelblue","yellow"),guide=FALSE) +
    geom_abline(slope = 0,intercept = 0) +
    labs(title = paste(gfln,'replicats ,betabinom vs binom')) +
    scale_y_continuous(breaks = seq(-0.75, 0.75, by= 0.25), limits=c(-.8,0.8))
  
  
  gfl <- exph_trueRM14A9A_307
  gfln <- "18" 
  tmp <- gfl[smp,] %>% arrange(`log.ecis`)
  tmp1 <- tmp[,1:3]
  tmp2 <- tmp[,7:9]
  names(tmp2) <- names(tmp1)
  null_two <- rbind(cbind(tmp1,type="Betabinom",od=1:siz),cbind(tmp2,type='binom',od=1:siz)) 
  plot03 <- ggplot(null_two,aes(x=od)) + 
    geom_ribbon(aes(ymin=ecisL,ymax=ecisH,fill=type),alpha=0.3)+
    scale_fill_manual(values=c( "steelblue","yellow")) +
    geom_abline(slope = 0,intercept = 0) +
    labs(title = paste(gfln,'replicats ,betabinom vs binom')) +
    scale_y_continuous(breaks = seq(-0.75, 0.75, by= 0.25), limits=c(-.8,0.8))
  
  
  grid.arrange(plot01,plot02,plot03,ncol=3)
}

fpo_plot <- function(){
  varName <- "exph_trueRM14A9A_"
  eachRep <- 18
  nRep <- 18
  repGenes <- 4710
  
  
  BB_matrix <- matrix(nrow=eachRep,ncol=nRep)
  BI_matrix <- matrix(nrow=eachRep,ncol=nRep)
  for(i in 1:nRep) {
    for (j in 1:eachRep) {
      num <- ifelse( i < nRep, (i-1) * eachRep + j,  (i-1)*eachRep + 1 ) 
      myVarName <- paste0(varName,formatC(num,width=3,flag="0"))
      tmp <- get(myVarName)
      BB_matrix[j,i] <- sum(tmp$ecisL > 0 | tmp$ecisH < 0,na.rm = T)/repGenes
      BI_matrix[j,i] <- sum(tmp$B.ecisL > 0 | tmp$B.ecisH <0,na.rm = T)/repGenes
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
  
  dfmm <- rbind(cbind(BB_df_m,type="Beta-B"),cbind(BI_df_m,type="Bino")) 
  
  ggplot(dfmm, aes(x=Used_replicates, y=fp_rate,color=type)) + 
    scale_color_manual(values=c( "steelblue","orange")) +
    geom_jitter(alpha=0.3, width=0.1,size=1.5) +
    geom_point(stat="summary", fun.y="mean",size=2.5) + 
    geom_errorbar(stat="summary", fun.data="mean_se", fun.args = list(mult = 1.96), width=0,size=1) + 
    scale_y_continuous(breaks=seq(0.1,0.8,0.05)) + 
    labs(title=paste(varName, "Significant genes found Rate vs. Number of replicates used for two models")) 
  #theme_bw() +
  #theme(axis.text.x=element_text(angle=45, hjust=1)) 
}



### e verse. averge read count graph ######### 
tmp <- exph %>% mutate(avg  = rowMeans(.[-1]))
plot(exph_trueRM14A9A_307$log.ecis,tmp$avg)
tmp2 <- cbind(exph_trueRM14A9A_307[,c(1,2,3)], avg = tmp$avg)
tmp2 <- tmp2 %>% mutate(disR = ecisH - ecisL, signi = (ecisL > 0 | ecisH < 0) )  

ggplot(tmp2,aes(x=log.ecis,y=log2(avg),color=disR,shape=signi)) + 
  geom_point(alpha=0.4,stroke=0,size=5) + 
    scale_color_continuous(low="pink",high="black") 

tmp3 <- cbind(tmp2, numofSig = cumsum(ifelse(is.na(tmp2$signi), 0, tmp2$signi)) + tmp2$signi*0)
plot(log2(tmp3$avg),tmp3$numofSig)
plot(1:4710,tmp3$numofSig)
abline(a = 0,b=3308/4710)

