source("~/cloud/project/F1_snpfilter/5_5Bsnp/mapping_res_alys.R") # get calih data.frame

source("~/cloud/project/F0_otherRaw/func.R")
source("~/cloud/project/F0_otherRaw/func_lklhratio.R")

calh_notail <- calh_notail %>% arrange(rowSums(.[-1])) # sort by average read depth 
calh_notail <- calh_notail[,c("group",grep(".*yGc",names(calh_notail),value=T),grep(".*rGc",names(calh_notail),value=T))]
# some genes because of python script has "-1" in them. so change them to zero 
# calih <- cbind(geneName=calih[,1],calih[,-1] + ( calih[,-1] < 0 ) * -(calih[,-1]))



calh_notail_Ct <- colSums(calh_notail[,-1], na.rm = T)
#save(list = c('rep44','depth44','expi','cali','expi_Ct','cali_Ct','paraGetB','neglhBinomial_cisonly','paraGetBB','neglhbetaBinomial_cisonly','dbetabinom'),file = '/cloud/project/hpc_pre/hpc_need_Jan14.RData')

save(list = c('calh_notail','calh_notail_Ct','paraGetB','neglhBinomial_cisonly','paraGetBB','neglhbetaBinomial_cisonly','dbetabinom',"paraGetB_lkr","paraGetBB_lkr"),file = '~/cloud/project/hpc_pre/calh_notail/hpc_need.RData')

