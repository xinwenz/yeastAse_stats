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

library(dplyr)
siz = 200
smp <- sample(x = 1:4306,size = siz,replace = F)
tmp1 <- expi_yps_null_10[smp,7:9] %>% arrange(`log.ecis`)
tmp2 <- expi_yps_null_02[smp,7:9] %>% arrange(B.log.ecis)
names(tmp2) <- names(tmp1)
null_two <- rbind(cbind(tmp1,type="high_rep_binom",od=1:siz),cbind(tmp2,type='low_rep_binom',od=1:siz)) 
ggplot(null_two,aes(x=od)) + 
  geom_ribbon(aes(ymin=ecisL,ymax=ecisH,fill=type),alpha=0.3)+
  scale_fill_manual(values=c( "steelblue","yellow")) +
  geom_abline(slope = 0,intercept = 0)

for(i in 1:10) {
  myVarName <- paste0("expi_yps_null_",formatC(i,width=2,flag="0"))
  print(myVarName)
  tmp <- get(myVarName)
  ans_tmp <- sum(tmp$ecisL > 0 | tmp$ecisH < 0,na.rm = T) 
  print(ans_tmp)
  print(ans_tmp/4306)
}

for(i in 1:10) {
  myVarName <- paste0("expi_yps_null_",formatC(i,width=2,flag="0"))
  print(myVarName)
  tmp <- get(myVarName)
  ans_tmp <- sum(tmp$B.ecisL > 0 | tmp$B.ecisH < 0,na.rm = T) 
  print(ans_tmp)
  print(ans_tmp/4306)
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
