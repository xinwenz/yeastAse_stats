################# expression level v.s. log(ecis) figure  #### 

tmpp <- exph_ypsRM14A9A_esum_2992

count_dat <- exph

avg_expression <- rowMeans(count_dat[,-1])
conf_range_BB <- tmpp$ecisH - tmpp$ecisL
conf_range_BI <- tmpp$B.ecisH - tmpp$B.ecisL
sigi <- (tmpp$ecisH < 0 | tmpp$ecisL > 0) 
tmpp_new <- cbind(tmpp,avg_expression,conf_range_BB,conf_range_BI, sigi)

tmpp_new_sig <- tmpp_new %>% filter(sigi == TRUE)
ggplot(tmpp_new_sig,aes(x=log.ecis,y=log2(avg_expression))) + geom_point(aes(color= conf_range_BB),size=3,alpha=0.4,shape=17) + scale_colour_gradient(low = "orange", high = "black")


tmpp_new_no <- tmpp_new %>% filter(sigi == FALSE | is.na(sigi))
ggplot(tmpp_new_no,aes(x=log.ecis,y=log2(avg_expression))) + geom_point(aes(color= conf_range_BB),size=3,alpha=0.4,shape=19) + scale_colour_gradient(low = "steelblue2", high = "black")



ggplot(tmpp_new,aes(x=log.ecis,y=log2(avg_expression))) + geom_point(aes(color= conf_range_BB, shape=sigi),size=3,alpha=0.4) + scale_colour_gradient(low = "cyan4", high = "black")



########## check the speed of finding significant genes ########## 
#tm <- exph_trueRM14A9A_esum_2700
#tm <- simu_p8_3750

tm <- exph_ypsRM14A9A_esum_2999
sigi <- (tm$ecisH < 0 | tm$ecisL > 0) 
tm_new <- cbind(tm,sigi)
tm_rat <- ( cumsum(ifelse(is.na(tm_new$sigi),0,tm_new$sigi)) + x*0 ) 
tm_final <- tail(tm_rat,1)

plot(1:nrow(tm), tm_rat/tm_final)
abline(a = 0,b= 1/nrow(tm))



######################################## test consequences of using less replicates 
nam <- names(exp_fine4)[2:21]
replicates <- seq(2,20,2)
lev = 0.95
for(i in seq(2,20,2) ){
  
  print(i)
  ale1n <- sample(nam,i,replace = F)
  ale2n <- gsub("^y","r",ale1n)
  exp_tmp <- exp_fine4[,c(ale1n,ale2n)]
  
  y_key <- grep('^y.*[AH]$',names(exp_tmp),value=T)
  r_key <- grep('^r.*[AH]$',names(exp_tmp),value=T)
  totCon <- colSums(exp_tmp[,y_key],na.rm =T) + colSums(exp_tmp[,r_key],na.rm=T)
  print(totCon)
  
  cl <- makeCluster(3)
  registerDoParallel(cl)# initiate cluster # not for windows users
  clusterExport(cl,c("dbetabinom","neglhbetaBinomial_cisonly"))
  ans_tmp <- foreach(d=iter(exp_tmp,by="row"),
                     .combine=rbind,
                     .packages='bbmle') %dopar%
    paraGetBB(d)
  stopCluster(cl)
  
  colnames(ans_tmp) <- c("log.ecis","ecis-","ecis+","log.rHy","rHy-","rhy+")
  
  ans_nam <- paste0("ans_bb_withrep",i)
  assign(ans_nam,cbind(exp_fine4[,1],ans_tmp ))  
  
  cl <- makeCluster(3)
  registerDoParallel(cl)# initiate cluster # not for windows users
  clusterExport(cl,c("paraGetB","neglhBinomial_cisonly"))
  ans_B_tmp <- foreach(d=iter(exp_tmp,by="row"),
  					 .combine=rbind,
  					 .packages='bbmle') %dopar%
  	paraGetB(d)
  stopCluster(cl)
  
  colnames(ans_B_tmp) <- c("log.ecis","ecis-","ecis+")
  ans_nam <- paste0("ans_b_withrep",i)
  assign(ans_nam,cbind(exp_fine4[,1],ans_B_tmp ))  
}



##################  beta binomial model --- result graph #############

cl <- makeCluster(3)
registerDoParallel(cl)# initiate cluster # not for windows users
clusterExport(cl,c("paraGetB","neglhBinomial_cisonly"))

y_key <- grep('^y.*[AH]$',names(exp_final),value=T)
r_key <- grep('^r.*[AH]$',names(exp_final),value=T)

ans_B <- foreach(d=iter(exp_final,by="row"),
                 .combine=rbind,
                 .verbose=F,
                 .packages='bbmle') %dopar%
  paraGetB(d)
stopCluster(cl)

colnames(ans_B) <- c("log.ecis","ecis-","ecis+")
ans_B_geneID <- cbind(exp_final[,c(1,2,3)],ans_B )  
ans_B_geneID_notcover0 <- subset(ans_B_geneID,ans_B_geneID$`ecis-`> 0 | ans_B_geneID$`ecis+`< 0 )   
print(nrow(ans_B_geneID_notcover0))

##########################
rep_B<- seq(2,20,2)
genes_sig_B <- c()

for(i in seq(2,20,2) ){
  
  print(i)
  cl <- makeCluster(3)
  registerDoParallel(cl)# initiate cluster # not for windows users
  clusterExport(cl,c("paraGetB","neglhBinomial_cisonly"))
  exp_tmp <- exp_final[,1:(3+3*i)]
  
  y_key <- grep('^y.*[AH]$',names(exp_tmp),value=T)
  r_key <- grep('^r.*[AH]$',names(exp_tmp),value=T)
  
  ans_B_tmp <- foreach(d=iter(exp_tmp,by="row"),
                       .combine=rbind,
                       .packages='bbmle') %dopar%
    paraGetB(d)
  stopCluster(cl)
  
  colnames(ans_B_tmp) <- c("log.ecis","ecis-","ecis+")
  ans_nam <- paste0("ans_B_tmp",i)
  assign(ans_nam,cbind(exp_final[,c(1,2,3)],ans_B_tmp ))  
  
  aB <- subset(ans_B_tmp,ans_B_tmp[,"ecis-"] > 0 | ans_B_tmp[,"ecis+"] < 0 )   
  
  sig_num <- nrow(aB)
  
  print(sig_num)
  genes_sig_B <- c(genes_sig_B,sig_num)
}

plot(rep_B, genes_sig_B,type='l',col="blue",ylim=c(3200,3700),xaxt="n")
axis(side = 1, at = rep_B,labels = T)
lines(rep_B, genes_sig,col="green")







################## check true postive and false positive and gene missing using ans_B_tmp and ans_tmp ############## 
hard <- 0
ans_bb_20_sig <- subset(ans_bb_withrep20,ans_bb_withrep20$`ecis-`> hard | ans_bb_withrep20$`ecis+`< (-1 * hard) )

gold_standlist <- ans_bb_20_sig$ypsTript
gold_n <- length(gold_standlist)

found_genes <- c()
missing_genes <- c()
false_positive_genes <- c()
inner_genes <- c()

for(i in seq(2,20,2)) {
  tigname <- paste0("ans_bb_withrep",i)
  print(tigname)
  
  tig <- get(tigname)
  abtig <- subset(tig,tig$"ecis-" > hard | tig$"ecis+" < (-1*hard) )
  found <- nrow(abtig)
  print("sig genes found: ")
  print(found)
  found_genes <- c(found_genes,found)
  
  tig_list <- abtig$ypsTript
  inner <- intersect(gold_standlist,tig_list)
  inner_ele <- length(inner)
  inner_genes <- c(inner_genes,inner_ele)
  missing <- gold_n - inner_ele
  print("missing:")
  print(missing)
  missing_genes <- c(missing_genes,missing)
  
  false_positive <- found - length(inner)
  print("false positive, not in gold stand")
  print(false_positive)
  false_positive_genes <-c(false_positive_genes,false_positive)
}


###########3 check binoomial model #########

b_found_genes <- c()
b_missing_genes <- c()
b_false_positive_genes <- c()
b_inner_genes <- c()

for(i in seq(2,20,2)) {
	tigname <- paste0("ans_b_withrep",i)
	print(tigname)
	
	tig <- get(tigname)
	abtig <- subset(tig,tig$"ecis-" > hard| tig$"ecis+" < (-1 * hard) )
	found <- nrow(abtig)
	print("sig genes found: ")
	print(found)
	b_found_genes <- c(b_found_genes,found)
	
	tig_list <- abtig$ypsTript
	inner <- intersect(gold_standlist,tig_list)
	inner_ele <- length(inner)
	b_inner_genes <- c(b_inner_genes,inner_ele)
	missing <- gold_n - inner_ele
	print("missing:")
	print(missing)
	b_missing_genes <- c(b_missing_genes,missing)
	
	false_positive <- found - length(inner)
	print("false positive, not in gold stand")
	print(false_positive)
	b_false_positive_genes <-c(b_false_positive_genes,false_positive)
}



#plot(seq(2,20,2),missing_genes,type='l',col='red',ylim=c(0,1500),xaxt="n")
#axis(side = 1,at=1:20,labels = T)
#axis(side = 2,at=seq(1,1500,by=250))
#lines(seq(2,20,2),b_missing_genes,col='orange')
#lines(seq(2,20,2),false_positive_genes,col='black')
#lines(seq(2,20,2),b_false_positive_genes,col='blue')
#legend(10, 1500, legend=c("number of missing genes", "missing genes Binom", "number of false positives", "false positives Binom"), 
 #      col=c("red", "orange","black","blue"), lty=c(1,1), cex=0.8)



plot(seq(2,20,2),inner_genes/found_genes[10],type='l',col='red',ylim=c(0,1),xaxt="n")
axis(side = 1,at=1:20,labels = T)
axis(side = 2,at=seq(0,1,by=20))
lines(seq(2,20,2),b_inner_genes/found_genes[10],col='orange')

lines(seq(2,20,2),false_positive_genes/(4528-found_genes[10]),col='black') # false postive rate 
lines(seq(2,20,2),false_positive_genes/found_genes,col='black',lty=3) # false discovery rate

lines(seq(2,20,2),b_false_positive_genes/(4528- found_genes[10]),col='blue')
lines(seq(2,20,2),b_false_positive_genes/b_found_genes,col='blue',lty=3)
legend(10, 0.75, legend=c("beta-binomial:power", "binomial:power", "beta-binomial:false postive rate", "Binom:false positive rate "), 
	   col=c("red", "orange","black","blue"), lty=c(1,1), cex=0.65)