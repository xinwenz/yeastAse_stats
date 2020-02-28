#!/bin/bash
#$ -t 1-360
#$ -ckpt restart

margs=$(head -n $SGE_TASK_ID expdc_hybridRM14A9A_hpc_com.txt | tail -n 1) 
module load R/3.5.1
Rscript -<<EOF $margs

load("hpc_need.RData") # likelihood functions,paraGetB/BB,expi,cali,expi_Ct,cali_Ct  
require("foreach")
require("iterators")
library("bbmle")

myargs <- commandArgs(trailingOnly = T)
num <- myargs[1]
#num <- "test_allsame"
sp_i <- as.integer(myargs[2:length(myargs)])
spPart1 <- sp_i[1:(length(sp_i)/2)]
spPart2 <- sp_i[(length(sp_i)/2+1) : length(sp_i)]

######################################### define which data frame and coverage vector to use 
orig <- expdc
orig_Ct <- expdc_Ct
f <- "expdc_liko"
#######################################################
x_key <- names(orig_Ct)[spPart1]
y_key <- names(orig_Ct)[spPart2]
#x_key <- names(orig_Ct)[3]
#y_key <- names(orig_Ct)[31]
ans1 <- foreach(d=iter(orig,by="row"),
                .combine=rbind,
                .packages='bbmle') %do%
    paraGetBB_lkr(d,x_key,y_key,orig_Ct)
colnames(ans1) <- c("BBf_ec","BBf_rHy","BBf_lik","BBec0_ec","BBec0_rHy","BBec0_lik")
rownames(ans1) <- unlist(orig[,1])
ans1 <- data.frame(ans1)


ans2 <- foreach(d=iter(orig,by="row"),
                .combine=rbind,
                .packages='bbmle') %do%
    paraGetB_lkr(d,x_key,y_key,orig_Ct)
colnames(ans2) <- c("Bf_ec","Bf_lik","Bec0_lik")
rownames(ans2) <- unlist(orig[,1])
ans2 <- data.frame(ans2)

ans3 <- cbind(ans1,ans2)

pval1 <- pchisq(2 * ( ans3[,"BBf_lik"] - ans3[,"BBec0_lik"] ) ,df =1 ,lower.tail = FALSE) # compare with ec == 0 
pval2 <- pchisq(2 * (ans3[,"BBf_lik"] - ans3[,"Bf_lik"] ) , df=1 , lower.tail = FALSE) # compare with binomial, check if rho is necessage.
pval3 <- pchisq(2 * ( ans3[,"Bf_lik"] - ans3[,"Bec0_lik"] ) , df =1 , lower.tail = FALSE) # in binomial with ec == 0 


ans4 <- cbind(ans3,pval1=pval1, pval2=pval2 , pval3 = pval3)

vn <- paste0(f,'_',num)
assign(x=vn,value=ans4)

save(list = vn,file=paste0(f,"_",num,".RData"))

EOF
