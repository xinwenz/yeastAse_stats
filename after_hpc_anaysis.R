setwd(dir = "/cloud/project/expi_true")
# 50 data points to calculate the range of false postives

for(i in 1:20) {
  mydataName <- paste0("expi_true_",formatC(i,width=2,flag="0"),".RData")
  load(mydataName)
}

# 20 replicates 
library(ggplot2)
## Fig1 and table 1 
ggplot(expi_true_20, aes(x=log.ecis)) + geom_histogram(binwidth = 0.05, color='blue',fill='blue',alpha=0.2)+ xlim(c(-2,2)) + labs(title="Histogram of expression log.ec ", x="log.ec", y="gene number") 


hist(expi_true_20$log.ecis,xlim=c(-1,1),breaks=10)
ans_tmp <- sum(expi_true_20$ecisL > 0 | expi_true_20$ecisH < 0,na.rm = T) 
bb_sig_95 <- c(bb_sig_95,ans_tmp)

sig_genes_0  <- subset(expi_true_20,expi_true_20$`ecisL`> 0 | expi_true_20$`ecisH`< 0 )  
sig_genes_2  <- subset(expi_true_20,expi_true_20$`ecisL`> 0.2 | expi_true_20$`ecisH`< -0.2 )  
sig_genes_5  <- subset(expi_true_20,expi_true_20$`ecisL`> 0.5 | expi_true_20$`ecisH`< -0.5 )  
sig_genes_1  <- subset(expi_true_20,expi_true_20$`ecisL`> 1 | expi_true_20$`ecisH`< -1 )  

c(nrow(sig_genes_0),nrow(sig_genes_2),nrow(sig_genes_5),nrow(sig_genes_1))
c(nrow(sig_genes_0),nrow(sig_genes_2),nrow(sig_genes_5),nrow(sig_genes_1))/4306


#### cali_null analysis #############  
setwd(dir = "/cloud/project/cali")
for(i in 1:20) {
  mydataName <- paste0("cali_null_",formatC(i,width=2,flag="0"),".RData")
  load(mydataName)
}

library(dplyr)
siz = 150
smp <- sample(x = 1:4423,size = siz,replace = F)
tmp1 <- cali_null_06[smp,1:3] %>% arrange(`log.ecis`)
tmp2 <- cali_null_06[smp,7:9] %>% arrange(B.log.ecis)
names(tmp2) <- names(tmp1)
null_two <- rbind(cbind(tmp1,type="bebebinomial",od=1:siz),cbind(tmp2,type='binomial',od=1:siz)) 
ggplot(null_two,aes(x=od)) + 
  geom_ribbon(aes(ymin=ecisL,ymax=ecisH,fill=type),alpha=0.3)+
  scale_fill_manual(values=c( "steelblue","yellow")) +
  geom_abline(slope = 0,intercept = 0)

for(i in 1:20) {
  myVarName <- paste0("cali_null_",formatC(i,width=2,flag="0"))
  print(myVarName)
  tmp <- get(myVarName)
  ans_tmp <- sum(tmp$ecisL > 0 | tmp$ecisH < 0,na.rm = T) 
  print(ans_tmp)
  print(ans_tmp/4423)
}

for(i in 1:20) {
  myVarName <- paste0("cali_null_",formatC(i,width=2,flag="0"))
  print(myVarName)
  tmp <- get(myVarName)
  ans_tmp <- sum(tmp$B.ecisL > 0 | tmp$B.ecisH < 0,na.rm = T) 
  print(ans_tmp)
  print(ans_tmp/4423)
}





sig_genes_0  <- subset(expi_true_20,expi_true_20$`ecisL`> 0 | expi_true_20$`ecisH`< 0 )  ### false positive rate ### 

#### expi_yps_null analysis #############
setwd(dir = "/cloud/project/expi_yps_null")
for(i in 1:10) {
  mydataName <- paste0("expi_yps_null_",formatC(i,width=2,flag="0"),".RData")
  load(mydataName)
}



for(i in 1:10) {
  myVarName <- paste0("expi_yps_null_",formatC(i,width=2,flag="0"))
  print(myVarName)
  tmp <- get(myVarName)
  ans_tmp <- sum(tmp$ecisL > 0 | tmp$ecisH < 0,na.rm = T) 
  print(ans_tmp)
  print(ans_tmp/4306)
}

for(i in 10:10) {
  myVarName <- paste0("expi_yps_null_",formatC(i,width=2,flag="0"))
  print(myVarName)
  tmp <- get(myVarName)
  ans_tmp <- subset(tmp, tmp$B.ecisL > 0 | tmp$B.ecisH < 0) 
  #print(ans_tmp)
  #print(ans_tmp/4306)
}


######### expi_rm_null analysis ############# 
setwd(dir = "/cloud/project/expi_rm_null")
for(i in 1:10) {
  mydataName <- paste0("expi_rm_null_",formatC(i,width=2,flag="0"),".RData")
  load(mydataName)
}

library(dplyr)
siz = 200
smp <- sample(x = 1:4306,size = siz,replace = F)
tmp1 <- expi_rm_null_10[smp,1:3] %>% arrange(`log.ecis`)
tmp2 <- expi_rm_null_02[smp,1:3] %>% arrange(B.log.ecis)
names(tmp2) <- names(tmp1)
null_two <- rbind(cbind(tmp1,type="high_rep_bb",od=1:siz),cbind(tmp2,type='low_rep_bb',od=1:siz)) 
ggplot(null_two,aes(x=od)) + 
  geom_ribbon(aes(ymin=ecisL,ymax=ecisH,fill=type),alpha=0.3)+
  scale_fill_manual(values=c( "steelblue","yellow")) +
  geom_abline(slope = 0,intercept = 0)

for(i in 1:10) {
  myVarName <- paste0("expi_rm_null_",formatC(i,width=2,flag="0"))
  print(myVarName)
  tmp <- get(myVarName)
  ans_tmp <- sum(tmp$ecisL > 0 | tmp$ecisH < 0,na.rm = T) 
  print(ans_tmp)
  print(ans_tmp/4306)
}

for(i in 1:10) {
  myVarName <- paste0("expi_rm_null_",formatC(i,width=2,flag="0"))
  print(myVarName)
  tmp <- get(myVarName)
  ans_tmp <- sum(tmp$B.ecisL > 0 | tmp$B.ecisH < 0,na.rm = T) 
  print(ans_tmp)
  print(ans_tmp/4306)
}
############ expi true significant rate #######################  
setwd(dir = "~/cloud/project/expi_true/")
# 50 data points to calculate the range of false postives
for(i in 1:20) {
  mydataName <- paste0("expi_true_",formatC(i,width=2,flag="0"),".RData")
  load(mydataName)
}
expi_sig_num <- function(){
  ans <- c()
  for(i in 1:20) {
    myVarName <- paste0("expi_true_",formatC(i,width=2,flag="0"))
    tmp <- get(myVarName)
    B_uniq <- sum(tmp$B.ecisL > 0 | tmp$B.ecisH < 0,na.rm = T) /nrow(tmp)
    Bb_uniq <- sum(tmp$ecisL > 0 | tmp$ecisH < 0,na.rm = T) /nrow(tmp)
    
    
    ans <- rbind(ans,c(ans_expi_true_B,ans_expi_true_Bb))
  }
  rownames(ans) <- 1:20
  colnames(ans) <- c("expi_true_binom","expi_true_betaBi")
  return(data.frame(ans))
}
ggplot(expi_sig_num(), aes(1:20)) + 
  geom_line(aes(y = expi_true_binom), linetype="longdash") + 
  geom_line(aes(y = expi_true_betaBi),linetype = "solid")

########check false postive Rate in NULL datasets #######################
get_fal_pos <- function() {
  ans <- c()
  for(i in 1:22) {
  ans_expi_rm_B <- tryCatch({
    myVarName <- paste0("expi_rm_",formatC(i,width=2,flag="0"))
    tmp <- get(myVarName)
    sum(tmp$B.ecisL > 0 | tmp$B.ecisH < 0,na.rm = T) /nrow(tmp)
  },  error = function(e) {
    NA
  }, finally = {
    NA
  })
  
  ans_expi_rm_Bb <- tryCatch({
    myVarName <- paste0("expi_rm_",formatC(i,width=2,flag="0"))
    tmp <- get(myVarName)
    sum(tmp$ecisL > 0 | tmp$ecisH < 0,na.rm = T) /nrow(tmp)
  },  error = function(e) {
    NA
  }, finally = {
    NA
  })
  
  ans_expi_yps_B <- tryCatch({
    myVarName <- paste0("expi_yps_",formatC(i,width=2,flag="0"))
    tmp <- get(myVarName)
    sum(tmp$B.ecisL > 0 | tmp$B.ecisH < 0,na.rm = T) /nrow(tmp)
  },  error = function(e) {
    NA
  }, finally = {
    NA
  })
  
  ans_expi_yps_Bb <- tryCatch({
    myVarName <- paste0("expi_yps_",formatC(i,width=2,flag="0"))
    tmp <- get(myVarName)
    sum(tmp$ecisL > 0 | tmp$ecisH < 0,na.rm = T) /nrow(tmp)
  },  error = function(e) {
    NA
  }, finally = {
    NA
  })
  
  
  myVarName <- paste0("rep_null_",formatC(i,width=2,flag="0"))
  tmp <- get(myVarName)
  ans_rep_null_B <- sum(tmp$B.ecisL > 0 | tmp$B.ecisH < 0,na.rm = T) /nrow(tmp)
  ans_rep_null_Bb <- sum(tmp$ecisL > 0 | tmp$ecisH < 0,na.rm = T) /nrow(tmp)
  

  ans_cali_null_B <- tryCatch({
    myVarName <- paste0("cali_null_",formatC(i,width=2,flag="0"))
    tmp <- get(myVarName)
    sum(tmp$B.ecisL > 0 | tmp$B.ecisH < 0,na.rm = T) /nrow(tmp)
  },  error = function(e) {
    NA
  }, finally = {
    NA
  })
  
  ans_cali_null_Bb <- tryCatch({
    myVarName <- paste0("cali_null_",formatC(i,width=2,flag="0"))
    tmp <- get(myVarName)
    sum(tmp$ecisL > 0 | tmp$ecisH < 0,na.rm = T) /nrow(tmp)
  },  error = function(e) {
    NA
  }, finally = {
    NA
  })
  
  
   ans <- rbind(ans,c(ans_expi_rm_B,ans_expi_rm_Bb, ans_expi_yps_B,ans_expi_yps_Bb,ans_rep_null_B,ans_rep_null_Bb,ans_cali_null_B,ans_cali_null_Bb))
}
  rownames(ans) <- 1:22
  colnames(ans) <- c("rm_null_binom","rm_null_betaBi","yps_null_binom","yps_null_betaBi","rep_null_binom","rep_null_betaBi","cali_null_binom","cali_null_betaBi")
return(ans)
}

ans <- get_fal_pos()
df_ans <- data.frame(ans)
ggplot(df_ans, aes(1:22)) + 
  #geom_line(aes(y = rm_null_binom, colour = "rm_null"), linetype="longdash") + 
  #geom_line(aes(y = rm_null_betaBi, colour = "rm_null"),linetype = "solid") +
  geom_line(aes(y = yps_null_binom, colour = "yps_null"), linetype="longdash") + 
  geom_line(aes(y = yps_null_betaBi, colour = "yps_null"),linetype = "solid") +
  geom_line(aes(y = rep_null_binom, colour = "rep_null"), linetype="longdash") + 
  geom_line(aes(y = rep_null_betaBi, colour = "rep_null"),linetype = "solid") +
  geom_line(aes(y = cali_null_binom, colour = "cali_null"), linetype="longdash") + 
  geom_line(aes(y = cali_null_betaBi, colour = "cali_null"),linetype = "solid") + 
  labs(title='number of significant genes in NULL datasets with increasing bio-replicates')

############ check expi_rm_null  ###########
setwd(dir = "~/cloud/project/expi_rm/")
# 50 data points to calculate the range of false postives
for(i in 1:10) {
  mydataName <- paste0("expi_rm_",formatC(i,width=2,flag="0"),".RData")
  load(mydataName)
}

plot_expi_rm <- function() {
  gfl <- expi_rm_02
  gfln <- "expi_rm_02"
  siz = 150
  smp <- sample(x = 1:nrow(gfl),size = siz,replace = F)
  tmp <- gfl[smp,] %>% arrange(`log.ecis`)
  tmp1 <- tmp[,1:3]
  tmp2 <- tmp[,7:9]
  names(tmp2) <- names(tmp1)
  null_two <- rbind(cbind(tmp1,type="Betabinom",od=1:siz),cbind(tmp2,type='binom',od=1:siz)) 
  plot05 <- ggplot(null_two,aes(x=od)) + 
    geom_ribbon(aes(ymin=ecisL,ymax=ecisH,fill=type),alpha=0.3)+
    scale_fill_manual(values=c( "steelblue","yellow"),guide=FALSE) +
    geom_abline(slope = 0,intercept = 0) +
    labs(title = paste(gfln,'replicats ,betabinom vs binom')) + 
    scale_y_continuous(breaks = seq(-0.75, 0.75, by= 0.25), limits=c(-0.8,0.8))
  
  
  gfl <- expi_rm_05
  gfln <- "expi_rm_05" 
  smp <- sample(x = 1:nrow(gfl),size = siz,replace = F)
  tmp <- gfl[smp,] %>% arrange(`log.ecis`)
  tmp1 <- tmp[,1:3]
  tmp2 <- tmp[,7:9]
  names(tmp2) <- names(tmp1)
  null_two <- rbind(cbind(tmp1,type="Betabinom",od=1:siz),cbind(tmp2,type='binom',od=1:siz)) 
  plot10 <- ggplot(null_two,aes(x=od)) + 
    geom_ribbon(aes(ymin=ecisL,ymax=ecisH,fill=type),alpha=0.3)+
    scale_fill_manual(values=c( "steelblue","yellow"),guide=FALSE) +
    geom_abline(slope = 0,intercept = 0) +
    labs(title = paste(gfln,'replicats ,betabinom vs binom')) +
    scale_y_continuous(breaks = seq(-0.75, 0.75, by= 0.25), limits=c(-.8,0.8))
  
  gfl <- expi_rm_10
  gfln <- "expi_rm_10"
  smp <- sample(x = 1:nrow(gfl),size = siz,replace = F)
  tmp <- gfl[smp,] %>% arrange(`log.ecis`)
  tmp1 <- tmp[,1:3]
  tmp2 <- tmp[,7:9]
  names(tmp2) <- names(tmp1)
  null_two <- rbind(cbind(tmp1,type="Betabinom",od=1:siz),cbind(tmp2,type='binom',od=1:siz)) 
  plot20 <- ggplot(null_two,aes(x=od)) + 
    geom_ribbon(aes(ymin=ecisL,ymax=ecisH,fill=type),alpha=0.3)+
    scale_fill_manual(values=c( "steelblue","yellow")) +
    geom_abline(slope = 0,intercept = 0) +
    labs(title = paste(gfln,'replicats ,betabinom vs binom')) +
    scale_y_continuous(breaks = seq(-0.75, 0.75, by= 0.25), limits=c(-0.8,0.8))
  
  grid.arrange(plot05,plot10,plot20,ncol=3)
}

############ ckeck expi_yps_null ############ 
setwd(dir = "~/cloud/project/expi_yps/")
# 50 data points to calculate the range of false postives
for(i in 1:10) {
  mydataName <- paste0("expi_yps_",formatC(i,width=2,flag="0"),".RData")
  load(mydataName)
}

plot_expi_yps <- function() {
  gfl <- expi_yps_02
  gfln <- "expi_yps_02"
  siz = 150
  smp <- sample(x = 1:nrow(gfl),size = siz,replace = F)
  tmp <- gfl[smp,] %>% arrange(`log.ecis`)
  tmp1 <- tmp[,1:3]
  tmp2 <- tmp[,7:9]
  names(tmp2) <- names(tmp1)
  null_two <- rbind(cbind(tmp1,type="Betabinom",od=1:siz),cbind(tmp2,type='binom',od=1:siz)) 
  plot05 <- ggplot(null_two,aes(x=od)) + 
    geom_ribbon(aes(ymin=ecisL,ymax=ecisH,fill=type),alpha=0.3)+
    scale_fill_manual(values=c( "steelblue","yellow"),guide=FALSE) +
    geom_abline(slope = 0,intercept = 0) +
    labs(title = paste(gfln,'replicats ,betabinom vs binom')) + 
    scale_y_continuous(breaks = seq(-0.75, 0.75, by= 0.25), limits=c(-0.8,0.8))
  
  
  gfl <- expi_yps_05
  gfln <- "expi_yps_05" 
  smp <- sample(x = 1:nrow(gfl),size = siz,replace = F)
  tmp <- gfl[smp,] %>% arrange(`log.ecis`)
  tmp1 <- tmp[,1:3]
  tmp2 <- tmp[,7:9]
  names(tmp2) <- names(tmp1)
  null_two <- rbind(cbind(tmp1,type="Betabinom",od=1:siz),cbind(tmp2,type='binom',od=1:siz)) 
  plot10 <- ggplot(null_two,aes(x=od)) + 
    geom_ribbon(aes(ymin=ecisL,ymax=ecisH,fill=type),alpha=0.3)+
    scale_fill_manual(values=c( "steelblue","yellow"),guide=FALSE) +
    geom_abline(slope = 0,intercept = 0) +
    labs(title = paste(gfln,'replicats ,betabinom vs binom')) +
    scale_y_continuous(breaks = seq(-0.75, 0.75, by= 0.25), limits=c(-.8,0.8))
  
  gfl <- expi_yps_10
  gfln <- "expi_yps_10"
  smp <- sample(x = 1:nrow(gfl),size = siz,replace = F)
  tmp <- gfl[smp,] %>% arrange(`log.ecis`)
  tmp1 <- tmp[,1:3]
  tmp2 <- tmp[,7:9]
  names(tmp2) <- names(tmp1)
  null_two <- rbind(cbind(tmp1,type="Betabinom",od=1:siz),cbind(tmp2,type='binom',od=1:siz)) 
  plot20 <- ggplot(null_two,aes(x=od)) + 
    geom_ribbon(aes(ymin=ecisL,ymax=ecisH,fill=type),alpha=0.3)+
    scale_fill_manual(values=c( "steelblue","yellow")) +
    geom_abline(slope = 0,intercept = 0) +
    labs(title = paste(gfln,'replicats ,betabinom vs binom')) +
    scale_y_continuous(breaks = seq(-0.75, 0.75, by= 0.25), limits=c(-0.8,0.8))
  
  grid.arrange(plot05,plot10,plot20,ncol=3)
}

############ check cali_null data ############# 
setwd(dir = "~/cloud/project/cali_null")
# 50 data points to calculate the range of false postives
for(i in 1:20) {
  mydataName <- paste0("cali_null_",formatC(i,width=2,flag="0"),".RData")
  load(mydataName)
}

plot_cali_null <- function() {
  gfl <- cali_null_05
  gfln <- "cali_null_05"
  siz = 150
  smp <- sample(x = 1:nrow(gfl),size = siz,replace = F)
  tmp <- gfl[smp,] %>% arrange(`log.ecis`)
  tmp1 <- tmp[,1:3]
  tmp2 <- tmp[,7:9]
  names(tmp2) <- names(tmp1)
  null_two <- rbind(cbind(tmp1,type="Betabinom",od=1:siz),cbind(tmp2,type='binom',od=1:siz)) 
  plot05 <- ggplot(null_two,aes(x=od)) + 
    geom_ribbon(aes(ymin=ecisL,ymax=ecisH,fill=type),alpha=0.3)+
    scale_fill_manual(values=c( "steelblue","yellow"),guide=FALSE) +
    geom_abline(slope = 0,intercept = 0) +
    labs(title = paste(gfln,'replicats ,betabinom vs binom')) + 
    scale_y_continuous(breaks = seq(-0.75, 0.75, by= 0.25), limits=c(-0.8,0.8))
  
  
  gfl <- cali_null_10
  gfln <- "cali_null_10" 
  smp <- sample(x = 1:nrow(gfl),size = siz,replace = F)
  tmp <- gfl[smp,] %>% arrange(`log.ecis`)
  tmp1 <- tmp[,1:3]
  tmp2 <- tmp[,7:9]
  names(tmp2) <- names(tmp1)
  null_two <- rbind(cbind(tmp1,type="Betabinom",od=1:siz),cbind(tmp2,type='binom',od=1:siz)) 
  plot10 <- ggplot(null_two,aes(x=od)) + 
    geom_ribbon(aes(ymin=ecisL,ymax=ecisH,fill=type),alpha=0.3)+
    scale_fill_manual(values=c( "steelblue","yellow"),guide=FALSE) +
    geom_abline(slope = 0,intercept = 0) +
    labs(title = paste(gfln,'replicats ,betabinom vs binom')) +
    scale_y_continuous(breaks = seq(-0.75, 0.75, by= 0.25), limits=c(-.8,0.8))
  
  gfl <- cali_null_20
  gfln <- "cali_null_20"
  smp <- sample(x = 1:nrow(gfl),size = siz,replace = F)
  tmp <- gfl[smp,] %>% arrange(`log.ecis`)
  tmp1 <- tmp[,1:3]
  tmp2 <- tmp[,7:9]
  names(tmp2) <- names(tmp1)
  null_two <- rbind(cbind(tmp1,type="Betabinom",od=1:siz),cbind(tmp2,type='binom',od=1:siz)) 
  plot20 <- ggplot(null_two,aes(x=od)) + 
    geom_ribbon(aes(ymin=ecisL,ymax=ecisH,fill=type),alpha=0.3)+
    scale_fill_manual(values=c( "steelblue","yellow")) +
    geom_abline(slope = 0,intercept = 0) +
    labs(title = paste(gfln,'replicats ,betabinom vs binom')) +
    scale_y_continuous(breaks = seq(-0.75, 0.75, by= 0.25), limits=c(-0.8,0.8))
  
  grid.arrange(plot05,plot10,plot20,ncol=3)
}

###################### check rep44 data ############# 
setwd(dir = "~/cloud/project/rep_null")
# 50 data points to calculate the range of false postives
for(i in 1:22) {
  mydataName <- paste0("rep_null_",formatC(i,width=2,flag="0"),".RData")
  load(mydataName)
}

plot_rep_null <- function() {
siz = 150
smp <- sample(x = 1:nrow(rep44),size = siz,replace = F)
  
  
gfl <- rep_null_04
gfln <- "rep_null_04"
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


gfl <- rep_null_08
gfln <- "rep_null_08" 
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


gfl <- rep_null_16
gfln <- "rep_null_16" 
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

gfl <- rep_null_22
gfln <- "rep_null_22"
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
