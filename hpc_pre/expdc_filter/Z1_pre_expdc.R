load('~/cloud/project/F1_snpfilter/7_5Bsnp_exp/expdc.RData')

source("~/cloud/project/F0_otherRaw/func.R")
source("~/cloud/project/F0_otherRaw/func_lklhratio.R")


expdc_Ct <- colSums(expdc[,-1], na.rm = T)

save(list = c('expdc','expdc_Ct','paraGetB','neglhBinomial_cisonly','paraGetBB','neglhbetaBinomial_cisonly','dbetabinom',"paraGetB_lkr","paraGetBB_lkr"),file = '~/cloud/project/hpc_pre/expdc_filter/hpc_need.RData')
