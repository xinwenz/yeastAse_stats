load('~/cloud/project/F0_otherRaw/rep44.RData')

rep44 <- rep44 %>% arrange(rowSums(.[-1]))
depth44 <- colSums(rep44[,-1], na.rm = TRUE)

source("~/cloud/project/F0_otherRaw/func.R")
source("~/cloud/project/F0_otherRaw/func_lklhratio.R")



save(list = c('rep44','depth44','paraGetB','neglhBinomial_cisonly','paraGetBB','neglhbetaBinomial_cisonly','dbetabinom',"paraGetB_lkr","paraGetBB_lkr"),file = '~/cloud/project/hpc_pre/Gier2015/hpc_need_gier.RData')