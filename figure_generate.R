library(ggplot2)
library(dplyr)
library(gridExtra)
library(reshape2)

load_rep <- function(x) {
    mydir <- paste0("~/cloud/project/",x)
    setwd(dir = mydir)
}

fpo_plot <- function(nm,eachRep,nRep,repGenes,note){
    varName <- nm
    eachRep <- eachRep
    nRep <- nRep
    repGenes <- repGenes
    
    
    BB_matrix <- matrix(nrow=eachRep,ncol=nRep)
    BI_matrix <- matrix(nrow=eachRep,ncol=nRep)
    for(i in 1:nRep) {
        for (j in 1:eachRep) {
            #num <- ifelse( i < nRep, (i-1) * eachRep + j,  (i-1)*eachRep + 1 ) 
            num <- (i-1) * eachRep + j 
            myVarName <- paste0(varName,formatC(num,width=4,flag="0",format="fg"))
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
    
    ans <- ggplot(dfmm, aes(x=Used_replicates, y=fp_rate,color=type)) + 
        scale_color_manual(values=c( "steelblue","orange")) +
        #geom_jitter(alpha=0.1, width=0.1,size=1.5) +
        #geom_point(stat="summary", fun.y="mean",size=2.5) + 
        #geom_errorbar(stat="summary", fun.data="mean_se", fun.args = list(mult = 1.96), width=0,size=1) +
        geom_boxplot(coef=0,outlier.size = 0) +
        labs(title=paste(varName, "False Posistive Rate vs. Number of replicates used for two models; ",note)) + 
        theme_bw() +
        theme(axis.text.x=element_text(angle=45, hjust=1)) 
    return(ans)
}


#### Fig 1A 
load_rep("hpc_null_nostr")
for(i in 1:17500) {
    mydataName <- paste0("rep_null_",formatC(i,width=5,flag="0"),".RData")
    load(mydataName)
}

fpo_plot(nm="rep_null_",eachRep = 500 ,nRep = 35 ,repGenes = 6023,note='no_structure')

##### Fig 1B 
load_rep("hpc_yps149_nostr")
for(i in 1:3000) {
    mydataName <- paste0("exph_ypsRM14A9A_esum_",formatC(i,width=4,flag="0"),".RData")
    load(mydataName)
}
fpo_plot(nm="exph_ypsRM14A9A_esum_",eachRep = 150,nRep = 20,repGenes = 4710,note="no_structure")

###### Fig 2.A
load_rep("hpc_simu_p1")
for(i in 1:3750) {
    mydataName <- paste0("simu_p1_",formatC(i,width=4,flag="0"),".RData")
    load(mydataName)
}

fpo_plot(nm="simu_p1_",eachRep = 150 ,nRep = 25 ,repGenes = 2000,note='p=0.1')


###### Fig 2.B 
load_rep("hpc_simu_p4")
for(i in 1:3750) {
    mydataName <- paste0("simu_p4_",formatC(i,width=4,flag="0"),".RData")
    load(mydataName)
}

fpo_plot(nm="simu_p4_",eachRep = 150 ,nRep = 25 ,repGenes = 2000,note='p=0.4')

##### Fig 2.C
load_rep("hpc_simu_p8")
for(i in 1:3750) {
    mydataName <- paste0("simu_p8_",formatC(i,width=4,flag="0"),".RData")
    load(mydataName)
}

fpo_plot(nm="simu_p8_",eachRep = 150 ,nRep = 25 ,repGenes = 2000,note='p=0.8')

#### Fig 2.D 
load_rep("hpc_simu_pc")
for(i in 1:3750) {
    mydataName <- paste0("simu_pc_",formatC(i,width=4,flag="0"),".RData")
    load(mydataName)
}

fpo_plot(nm="simu_pc_",eachRep = 150 ,nRep = 25 ,repGenes = 2000,note='p=runif(0 - 0.8)')


##### Fig 3. 
load_rep("hpc_expTrueRM14A9A")
for(i in 1:2700) {
    mydataName <- paste0("exph_trueRM14A9A_esum_",formatC(i,width=4,flag="0"),".RData")
    load(mydataName)
}
fpo_plot(nm="exph_trueRM14A9A_esum_",eachRep = 150,nRep = 18,repGenes = 4710,note="no_structure")


##### Fig 4. 
tmpp <- exph_trueRM14A9A_esum_2700
#hist(tmpp$log.ecis, breaks=100,xlim=c(-3,3)) ## histgram 
ggplot(tmpp, aes(x=log.ecis)) + geom_histogram(binwidth = 0.05, color="black",fill='white') + theme_bw() 

sum(tmpp$ecisL > 0 | tmpp$ecisH <0,na.rm = T)
sum(tmpp$ecisL > 0.2 | tmpp$ecisH < -0.2 ,na.rm = T)
sum(tmpp$ecisL > 0.5 | tmpp$ecisH < -0.5,na.rm = T)
sum(tmpp$ecisL > 1 | tmpp$ecisH < -1 ,na.rm = T)

sum(tmpp$ecisL > 0 | tmpp$ecisH <0,na.rm = T) / nrow(tmpp)
sum(tmpp$ecisL > 0.2 | tmpp$ecisH < -0.2 ,na.rm = T) /nrow(tmpp)
sum(tmpp$ecisL > 0.5 | tmpp$ecisH < -0.5,na.rm = T) /nrow(tmpp)
sum(tmpp$ecisL > 1 | tmpp$ecisH < -1,na.rm = T ) /nrow(tmpp)


############## Fig 5. 
plot_rep_null <- function() {
    siz = 150
    smp <- sort(sample(x = 1:nrow(exph_trueRM14A9A_esum_0001),size = siz,replace = F))
    
    
    ## repilicate Number : 1 -- (1 * 150) 
    gfl <- exph_trueRM14A9A_esum_0150
    gfln <- "1"
    tmp <- gfl[smp,] %>% arrange(`log.ecis`)
    tmp1 <- tmp[,1:3]
    tmp2 <- tmp[,7:9]
    names(tmp2) <- names(tmp1)
    null_two <- rbind(cbind(tmp1,type="Betabinom",od=1:siz),cbind(tmp2,type='binom',od=1:siz)) 
    plot01 <- ggplot(null_two,aes(x=od)) + 
        geom_ribbon(aes(ymin=ecisL,ymax=ecisH,fill=type),alpha=0.3)+
        scale_fill_manual(values=c( "steelblue","yellow"),guide=FALSE) +
        geom_abline(slope = 0,intercept = 0) +
        labs(title = paste(gfln,'replicats ,betabinom vs binom'),x='gene rank') + 
        scale_y_continuous(breaks = seq(-0.75, 0.75, by= 0.25), limits=c(-0.8,0.8))
    
    
    ## replicate Numnber 8 
    gfl <- exph_trueRM14A9A_esum_0900
    gfln <- "6" 
    tmp <- gfl[smp,] %>% arrange(`log.ecis`)
    tmp1 <- tmp[,1:3]
    tmp2 <- tmp[,7:9]
    names(tmp2) <- names(tmp1)
    null_two <- rbind(cbind(tmp1,type="Betabinom",od=1:siz),cbind(tmp2,type='binom',od=1:siz)) 
    plot02 <- ggplot(null_two,aes(x=od)) + 
        geom_ribbon(aes(ymin=ecisL,ymax=ecisH,fill=type),alpha=0.3)+
        scale_fill_manual(values=c( "steelblue","yellow"),guide=FALSE) +
        geom_abline(slope = 0,intercept = 0) +
        labs(title = paste(gfln,'replicats ,betabinom vs binom'),x='gene rank') +
        scale_y_continuous(breaks = seq(-0.75, 0.75, by= 0.25), limits=c(-.8,0.8))
    
    
    gfl <- exph_trueRM14A9A_esum_1800
    gfln <- "12" 
    tmp <- gfl[smp,] %>% arrange(`log.ecis`)
    tmp1 <- tmp[,1:3]
    tmp2 <- tmp[,7:9]
    names(tmp2) <- names(tmp1)
    null_two <- rbind(cbind(tmp1,type="Betabinom",od=1:siz),cbind(tmp2,type='binom',od=1:siz)) 
    plot03 <- ggplot(null_two,aes(x=od)) + 
        geom_ribbon(aes(ymin=ecisL,ymax=ecisH,fill=type),alpha=0.3)+
        scale_fill_manual(values=c( "steelblue","yellow"),guide=FALSE) +
        geom_abline(slope = 0,intercept = 0) +
        labs(title = paste(gfln,'replicats ,betabinom vs binom'),x='gene rank') +
        scale_y_continuous(breaks = seq(-0.75, 0.75, by= 0.25), limits=c(-.8,0.8))
    
    gfl <- exph_trueRM14A9A_esum_2700
    gfln <- "18"
    tmp <- gfl[smp,] %>% arrange(`log.ecis`)
    tmp1 <- tmp[,1:3]
    tmp2 <- tmp[,7:9]
    names(tmp2) <- names(tmp1)
    null_two <- rbind(cbind(tmp1,type="Betabinom",od=1:siz),cbind(tmp2,type='binom',od=1:siz)) 
    plot04 <- ggplot(null_two,aes(x=od)) + 
        geom_ribbon(aes(ymin=ecisL,ymax=ecisH,fill=type),alpha=0.3)+
        scale_fill_manual(values=c( "steelblue","yellow"),guide=FALSE) +
        geom_abline(slope = 0,intercept = 0) +
        labs(title = paste(gfln,'replicats ,betabinom vs binom'),x='gene rank') +
        scale_y_continuous(breaks = seq(-0.75, 0.75, by= 0.25), limits=c(-0.8,0.8))
    
    grid.arrange(plot01,plot02,plot03,plot04,ncol=4)
}

##### Fig 6. #####
tmpp <- exph_trueRM14A9A_esum_2700
godlist0 <- which(tmpp$ecisL > 0 | tmpp$ecisH <0)
godlist0N <- (1:4710)[-godlist0]

godlist5 <- which(tmpp$ecisL > 0.5 | tmpp$ecisH < -0.5)
godlist5N <- (1:4710)[-godlist5]

godlistA <- which(tmpp$ecisL > 1 | tmpp$ecisH < -1 )
godlistAN <- (1:4710)[-godlistA]

transf <- function(clas) {
  if(clas == '0'){
    return(0)
  }else if(clas == '5') {
    return(0.5)
  } else if (clas == 'A') {
    return(1)
  }else {
    return(NA)
  }
}

cal_fp <- function(nm,eachRep,nRep,repGenes,clas){
    varName <- nm
    eachRep <- eachRep
    nRep <- nRep
    repGenes <- repGenes
    
    mygodlist <- get(paste0("godlist",clas))
    mygodlistN <- get(paste0('godlist',clas,'N'))
    BB_matrix <- matrix(nrow=eachRep,ncol=nRep)
    BI_matrix <- matrix(nrow=eachRep,ncol=nRep)
    
    
    for(i in 1:nRep) {
        for (j in 1:eachRep) {
            #num <- ifelse( i < nRep, (i-1) * eachRep + j,  (i-1)*eachRep + 1 ) 
            num <- (i-1) * eachRep + j 
            myVarName <- paste0(varName,formatC(num,width=4,flag="0",format="fg"))
            tmp <- get(myVarName)
            
            
            sigBB_list <- which(tmp$ecisL > transf(clas) | tmp$ecisH < -transf(clas))
            sigBI_list <- which(tmp$B.ecisL > transf(clas) | tmp$B.ecisH < -transf(clas))
            
            BB_matrix[j,i] <- length(intersect(sigBB_list, mygodlistN)) /length(mygodlistN)
            BI_matrix[j,i] <- length(intersect(sigBI_list, mygodlistN)) /length(mygodlistN)
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
    
    ans <- ggplot(dfmm, aes(x=Used_replicates, y=fp_rate,color=type)) + 
        scale_color_manual(values=c( "steelblue","orange")) +
        #geom_jitter(alpha=0.1, width=0.1,size=1.5) +
        #geom_point(stat="summary", fun.y="mean",size=2.5) + 
        #geom_errorbar(stat="summary", fun.data="mean_se", fun.args = list(mult = 1.96), width=0,size=1) +
        geom_boxplot(coef=0,outlier.size = 0) +
        labs(title=paste(varName, "False Posistive Rate vs. Number of replicates used for two models; ",clas)) + 
        theme_bw() +
        theme(axis.text.x=element_text(angle=45, hjust=1)) 
    return(ans)
}

cal_fp('exph_trueRM14A9A_esum_',150,18,4710,'0')
cal_fp('exph_trueRM14A9A_esum_',150,18,4710,'5')
cal_fp('exph_trueRM14A9A_esum_',150,18,4710,'A')



############ Fig 6.2 ##########  
cal_power <- function(nm,eachRep,nRep,repGenes,clas){
  varName <- nm
  eachRep <- eachRep
  nRep <- nRep
  repGenes <- repGenes
  
  mygodlist <- get(paste0("godlist",clas))
  mygodlistN <- get(paste0('godlist',clas,'N'))
  BB_matrix <- matrix(nrow=eachRep,ncol=nRep)
  BI_matrix <- matrix(nrow=eachRep,ncol=nRep)
  
  
  for(i in 1:nRep) {
    for (j in 1:eachRep) {
      #num <- ifelse( i < nRep, (i-1) * eachRep + j,  (i-1)*eachRep + 1 ) 
      num <- (i-1) * eachRep + j 
      myVarName <- paste0(varName,formatC(num,width=4,flag="0",format="fg"))
      tmp <- get(myVarName)
      
      
      sigBB_list <- which(tmp$ecisL > transf(clas) | tmp$ecisH < -transf(clas))
      sigBI_list <- which(tmp$B.ecisL > transf(clas) | tmp$B.ecisH < -transf(clas))
      
      BB_matrix[j,i] <- length(intersect(sigBB_list, mygodlist)) /length(mygodlist)
      BI_matrix[j,i] <- length(intersect(sigBI_list, mygodlist)) /length(mygodlist)
    }
    
  }
  
  rownames(BB_matrix) <- 1:eachRep
  BB_df <- data.frame(BB_matrix)
  colnames(BB_df) <- paste0("UsRp",1:nRep)
  BB_df_m <- melt(BB_df)
  colnames(BB_df_m) <- c("Used_replicates","power")
  
  
  rownames(BI_matrix) <- 1:eachRep
  BI_df <- data.frame(BI_matrix)
  colnames(BI_df) <- paste0("UsRp",1:nRep)
  BI_df_m <- melt(BI_df)
  colnames(BI_df_m) <- c("Used_replicates","power")
  
  dfmm <- rbind(cbind(BB_df_m,type="Beta-B"),cbind(BI_df_m,type="Bino")) 
  
  ans <- ggplot(dfmm, aes(x=Used_replicates, y=power,color=type)) + 
    scale_color_manual(values=c( "steelblue","orange")) +
    #geom_jitter(alpha=0.1, width=0.1,size=1.5) +
    #geom_point(stat="summary", fun.y="mean",size=2.5) + 
    #geom_errorbar(stat="summary", fun.data="mean_se", fun.args = list(mult = 1.96), width=0,size=1) +
    geom_boxplot(coef=0,outlier.size = 0) +
    labs(title=paste(varName, "Power vs. Number of replicates used for two models; ",clas)) + 
    theme_bw() +
    theme(axis.text.x=element_text(angle=45, hjust=1)) + 
    ylim(0,1)
  return(ans)
}

cal_power('exph_trueRM14A9A_esum_',150,18,4710,'0')
cal_power('exph_trueRM14A9A_esum_',150,18,4710,'5')
cal_power('exph_trueRM14A9A_esum_',150,18,4710,'A')


############  Fig 6.3 
cal_fd <- function(nm,eachRep,nRep,repGenes,clas){
  varName <- nm
  eachRep <- eachRep
  nRep <- nRep
  repGenes <- repGenes
  
  mygodlist <- get(paste0("godlist",clas))
  mygodlistN <- get(paste0('godlist',clas,'N'))
  BB_matrix <- matrix(nrow=eachRep,ncol=nRep)
  BI_matrix <- matrix(nrow=eachRep,ncol=nRep)
  
  
  for(i in 1:nRep) {
    for (j in 1:eachRep) {
      #num <- ifelse( i < nRep, (i-1) * eachRep + j,  (i-1)*eachRep + 1 ) 
      num <- (i-1) * eachRep + j 
      myVarName <- paste0(varName,formatC(num,width=4,flag="0",format="fg"))
      tmp <- get(myVarName)
      
      
      sigBB_list <- which(tmp$ecisL > transf(clas) | tmp$ecisH < -transf(clas))
      sigBI_list <- which(tmp$B.ecisL > transf(clas) | tmp$B.ecisH < -transf(clas))
      
      BB_matrix[j,i] <- length(intersect(sigBB_list, mygodlistN)) /length(sigBB_list)
      BI_matrix[j,i] <- length(intersect(sigBI_list, mygodlistN)) /length(sigBI_list)
    }
    
  }
  
  rownames(BB_matrix) <- 1:eachRep
  BB_df <- data.frame(BB_matrix)
  colnames(BB_df) <- paste0("UsRp",1:nRep)
  BB_df_m <- melt(BB_df)
  colnames(BB_df_m) <- c("Used_replicates","fd")
  
  
  rownames(BI_matrix) <- 1:eachRep
  BI_df <- data.frame(BI_matrix)
  colnames(BI_df) <- paste0("UsRp",1:nRep)
  BI_df_m <- melt(BI_df)
  colnames(BI_df_m) <- c("Used_replicates","fd")
  
  dfmm <- rbind(cbind(BB_df_m,type="Beta-B"),cbind(BI_df_m,type="Bino")) 
  
  ans <- ggplot(dfmm, aes(x=Used_replicates, y=fd,color=type)) + 
    scale_color_manual(values=c( "steelblue","orange")) +
    #geom_jitter(alpha=0.1, width=0.1,size=1.5) +
    #geom_point(stat="summary", fun.y="mean",size=2.5) + 
    #geom_errorbar(stat="summary", fun.data="mean_se", fun.args = list(mult = 1.96), width=0,size=1) +
    geom_boxplot(coef=0,outlier.size = 0) +
    labs(title=paste(varName, "False Discovery vs. Number of replicates used for two models; ",clas)) + 
    theme_bw() +
    theme(axis.text.x=element_text(angle=45, hjust=1)) + 
    ylim(0,0.25)
  return(ans)
}
cal_fd('exph_trueRM14A9A_esum_',150,18,4710,'0')
cal_fd('exph_trueRM14A9A_esum_',150,18,4710,'5')
cal_fd('exph_trueRM14A9A_esum_',150,18,4710,'A')


##########  Fig 6.4 correlation test ############  
tmpp <- exph_trueRM14A9A_esum_2700
godtf0 <- tmpp$ecisL > 0 | tmpp$ecisH <0
godtf5 <- tmpp$ecisL > 0.5 | tmpp$ecisH < -0.5
godtfA <- tmpp$ecisL > 1 | tmpp$ecisH < -1 

cal_cor <- function(nm,eachRep,nRep,repGenes,clas){
  varName <- nm
  eachRep <- eachRep
  nRep <- nRep
  repGenes <- repGenes
  
  mytf <- get(paste0("godtf",clas))

  BB_matrix <- matrix(nrow=eachRep,ncol=nRep)
  BI_matrix <- matrix(nrow=eachRep,ncol=nRep)
  
  
for(i in 1:nRep) {
    for (j in 1:eachRep) {
      #num <- ifelse( i < nRep, (i-1) * eachRep + j,  (i-1)*eachRep + 1 ) 
      num <- (i-1) * eachRep + j 
      myVarName <- paste0(varName,formatC(num,width=4,flag="0",format="fg"))
      tmp <- get(myVarName)
      
      
      sigBB_list <- tmp$ecisL > transf(clas) | tmp$ecisH < -transf(clas)
      sigBI_list <- tmp$B.ecisL > transf(clas) | tmp$B.ecisH < -transf(clas)
      
      BB_matrix[j,i] <- cor(mytf,sigBB_list,use="pairwise.complete.obs")
      BI_matrix[j,i] <- cor(mytf,sigBI_list,use="pairwise.complete.obs")
    }
}
  
  rownames(BB_matrix) <- 1:eachRep
  BB_df <- data.frame(BB_matrix)
  colnames(BB_df) <- paste0("UsRp",1:nRep)
  BB_df_m <- melt(BB_df)
  colnames(BB_df_m) <- c("Used_replicates","cor")
  
  
  rownames(BI_matrix) <- 1:eachRep
  BI_df <- data.frame(BI_matrix)
  colnames(BI_df) <- paste0("UsRp",1:nRep)
  BI_df_m <- melt(BI_df)
  colnames(BI_df_m) <- c("Used_replicates","cor")
  
  dfmm <- rbind(cbind(BB_df_m,type="Beta-B"),cbind(BI_df_m,type="Bino")) 
  
  ans <- ggplot(dfmm, aes(x=Used_replicates, y=cor,color=type)) + 
    scale_color_manual(values=c( "steelblue","orange")) +
    #geom_jitter(alpha=0.1, width=0.1,size=1.5) +
    #geom_point(stat="summary", fun.y="mean",size=2.5) + 
    #geom_errorbar(stat="summary", fun.data="mean_se", fun.args = list(mult = 1.96), width=0,size=1) +
    geom_boxplot(coef=0,outlier.size = 0) +
    labs(title=paste(varName, "correlation with gold standard; ",clas)) + 
    theme_bw() +
    theme(axis.text.x=element_text(angle=45, hjust=1)) + 
    ylim(0,1)
  return(ans)
}
cal_cor('exph_trueRM14A9A_esum_',150,18,4710,'0')
cal_cor('exph_trueRM14A9A_esum_',150,18,4710,'5')
cal_cor('exph_trueRM14A9A_esum_',150,18,4710,'A')


#########  the godA_gene_list ##########33  
tmppp <- cbind(tmpp,rowMeans(exph[,-1]))
godA_list <- tmppp[which(godtfA==TRUE),]
a <- rownames(godA_list)
sink("~/cloud/project/godA_list.txt")
for(i in 1:101) { cat(a[i],'\n') }
sink()

sink("~/cloud/project/all_genes.txt")
for(i in 1:4710) { cat(rownames(tmpp)[i],'\n')}
sink()
