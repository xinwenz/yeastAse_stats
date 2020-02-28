source("~/cloud/project/M2_cali_gene_input.R") # get calih data.frame
source("~/cloud/project/F0_otherRaw/func.R")

calih_Ct <- colSums(calih[,-1])
#save(list = c('rep44','depth44','expi','cali','expi_Ct','cali_Ct','paraGetB','neglhBinomial_cisonly','paraGetBB','neglhbetaBinomial_cisonly','dbetabinom'),file = '/cloud/project/hpc_pre/hpc_need_Jan14.RData')

save(list = c('calih','calih_Ct','paraGetB','neglhBinomial_cisonly','paraGetBB','neglhbetaBinomial_cisonly','dbetabinom'),file = '~/cloud/project/hpc_pre/calih_gene/hpc_need.RData')

