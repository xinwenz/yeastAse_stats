load_rep <- function(x) {
    mydir <- paste0("~/cloud/project/",x)
    setwd(dir = mydir)
}

for(i in 1:3300) {
    mydataName <- paste0("rep_null_",formatC(i,width=4,flag="0"),".RData")
    load(mydataName)
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





################# line extension ###############
fpo_number <- function(nm,eachRep,nRep,repGenes,note){
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
  
  #dfmm <- rbind(cbind(BB_df_m,type="Beta-B"),cbind(BI_df_m,type="Bino")) 
  
  #ans <- ggplot(dfmm, aes(x=Used_replicates, y=fp_rate,color=type)) + 
   # scale_color_manual(values=c( "steelblue","orange")) +
    #geom_jitter(alpha=0.1, width=0.1,size=1.5) +
    #geom_point(stat="summary", fun.y="mean",size=2.5) + 
    #geom_errorbar(stat="summary", fun.data="mean_se", fun.args = list(mult = 1.96), width=0,size=1) +
    #geom_boxplot(coef=0,outlier.size = 0) +
    #labs(title=paste(varName, "False Posistive Rate vs. Number of replicates used for two models; ",note)) + 
    #theme_bw() +
    #theme(axis.text.x=element_text(angle=45, hjust=1)) 
  return(BI_df)
}
ans_BI <- fpo_number(nm="exph_trueRM14A9A_esum_",eachRep = 150,nRep = 18,repGenes = 4710,note="no_structure")
time <- 1:18
fp_rate <- colMeans(ans)
lin.mod2 <- lm(fp_rate ~ I(time^0.5)+time)
pr.lim2 <- predict(lin.mod2)
plot(time,fp_rate)
lines(pr.lim2 ~ time)



