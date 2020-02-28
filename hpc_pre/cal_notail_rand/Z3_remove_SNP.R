load("~/cloud/project/F1_snpfilter/5_5Bsnp/cal_no_right_tail.RData") 
# get cal_no_right_tail data.frame

cal_with_label <- cal_no_right_tail[c(grep("^y.*[AH]_rc",names(cal_no_right_tail),value=T), grep("^r.*[AH]_rc",names(cal_no_right_tail),value=T))] 

x <- rowSums(cal_with_label[,c(1:20)[-13]],na.rm = TRUE)
y <- rowSums(cal_with_label[,c(1:20)[-13]+20],na.rm = TRUE)
pval_df <- data.frame(yps=x,ypsrm = x+y)
pres <- apply(pval_df, 1 , function(x) binom.test(x[1],x[2])$p.value)
pres_rank <- rank(pres)

cal_filter_10psnps <- cal_with_label[pres_rank > 66352*0.12 ,] 
hist(which(pres_rank < 66352 * 0.12)) # check the postion of those 

ax <- rowSums(cal_filter_10psnps[,c(1:20)[-13]],na.rm = TRUE)
ay <- rowSums(cal_filter_10psnps[,c(1:20)[-13]+20],na.rm = TRUE)
apval_df <- data.frame(yps=ax,ypsrm = ax+ay)
apres <- apply(apval_df, 1 , function(x) binom.test(x[1],x[2])$p.value)


a <- sample(size=nrow(cal_filter_10psnps), x=1:1000, replace=TRUE )
cal_filter_10psnps <- cbind(a,cal_filter_10psnps)

 # sort by average read depth 
cal_random_filter <- cal_filter_10psnps %>% group_by(a) %>% summarise_all(funs(sum))


source("~/cloud/project/F0_otherRaw/func.R")
source("~/cloud/project/F0_otherRaw/func_lklhratio.R")

# some genes because of python script has "-1" in them. so change them to zero 
# calih <- cbind(geneName=calih[,1],calih[,-1] + ( calih[,-1] < 0 ) * -(calih[,-1]))


cal_random_filter_Ct <- colSums(cal_random_filter[,-1], na.rm = T)
#save(list = c('rep44','depth44','expi','cali','expi_Ct','cali_Ct','paraGetB','neglhBinomial_cisonly','paraGetBB','neglhbetaBinomial_cisonly','dbetabinom'),file = '/cloud/project/hpc_pre/hpc_need_Jan14.RData')

save(list = c('cal_random_filter','cal_random_filter_Ct','paraGetB','neglhBinomial_cisonly','paraGetBB','neglhbetaBinomial_cisonly','dbetabinom',"paraGetB_lkr","paraGetBB_lkr"),file = '~/cloud/project/hpc_pre/cal_notail_rand/hpc_need_filter.RData')
