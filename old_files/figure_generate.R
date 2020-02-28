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


#### Fig 1A 



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

godA_list_sort <- godA_list %>% mutate(geneN = rownames(godA_list)) %>% arrange(desc(abs(log.ecis)))
