source("~/cloud/project/F1_snpfilter/5_5Bsnp/mapping_res_alys.R") # get calih data.frame

source("~/cloud/project/F0_otherRaw/func.R")
source("~/cloud/project/F0_otherRaw/func_lklhratio.R")

cal_agrg <- cal_agrg %>% arrange(rowSums(.[-1])) # sort by average read depth 
cal_agrg <- cal_agrg[,c("a",grep(".*yGc",names(calh_notail),value=T),grep(".*rGc",names(calh_notail),value=T))]

# some genes because of python script has "-1" in them. so change them to zero 
# calih <- cbind(geneName=calih[,1],calih[,-1] + ( calih[,-1] < 0 ) * -(calih[,-1]))


cal_agrg_Ct <- colSums(cal_agrg[,-1], na.rm = T)
#save(list = c('rep44','depth44','expi','cali','expi_Ct','cali_Ct','paraGetB','neglhBinomial_cisonly','paraGetBB','neglhbetaBinomial_cisonly','dbetabinom'),file = '/cloud/project/hpc_pre/hpc_need_Jan14.RData')

save(list = c('cal_agrg','cal_agrg_Ct','paraGetB','neglhBinomial_cisonly','paraGetBB','neglhbetaBinomial_cisonly','dbetabinom',"paraGetB_lkr","paraGetBB_lkr"),file = '~/cloud/project/hpc_pre/cal_notail_aggr/hpc_need.RData')
