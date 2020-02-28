load('~/cloud/project/F0_otherRaw/hpc_need_Feb4.RData')
source("~/cloud/project/F0_otherRaw/func.R")
source("~/cloud/project/F0_otherRaw/func_lklhratio.R")



save(list = c('exph','exph_Ct_esum','paraGetB','neglhBinomial_cisonly','paraGetBB','neglhbetaBinomial_cisonly','dbetabinom',"paraGetB_lkr","paraGetBB_lkr"),file = '~/cloud/project/hpc_pre/exph/hpc_need_exph.RData')
