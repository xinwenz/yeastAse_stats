load("~/cloud/project/F1_snpfilter/5_5Bsnp/cal_no_right_tail.RData") 
# get cal_no_right_tail data.frame

a <- sample(size=nrow(cal_no_right_tail), x=1:1000, replace=TRUE )
cal_with_label <- cbind(a,cal_no_right_tail)

cal_with_label <- cal_with_label[,c("a",grep("^y.*[AH]_rc",names(cal_no_right_tail),value=T), grep("^r.*[AH]_rc",names(cal_no_right_tail),value=T))]  # sort by average read depth 

cal_random <- cal_with_label %>% group_by(a) %>% summarise_all(funs(sum))

 
source("~/cloud/project/F0_otherRaw/func.R")
source("~/cloud/project/F0_otherRaw/func_lklhratio.R")

# some genes because of python script has "-1" in them. so change them to zero 
# calih <- cbind(geneName=calih[,1],calih[,-1] + ( calih[,-1] < 0 ) * -(calih[,-1]))


cal_random_Ct <- colSums(cal_random[,-1], na.rm = T)
#save(list = c('rep44','depth44','expi','cali','expi_Ct','cali_Ct','paraGetB','neglhBinomial_cisonly','paraGetBB','neglhbetaBinomial_cisonly','dbetabinom'),file = '/cloud/project/hpc_pre/hpc_need_Jan14.RData')

save(list = c('cal_random','cal_random_Ct','paraGetB','neglhBinomial_cisonly','paraGetBB','neglhbetaBinomial_cisonly','dbetabinom',"paraGetB_lkr","paraGetBB_lkr"),file = '~/cloud/project/hpc_pre/cal_notail_rand/hpc_need.RData')
