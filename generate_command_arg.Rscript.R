load('/cloud/project/otherRaw/run_hpc_Cor_seq.RData')
#cali_C order: 
cali_ordr <- cali_cor_match[,1]
sink('/cloud/project/hpc_pre/cali_hpc_com.txt')
for(i in 1:20) {
  input <- c(cali_ordr[1:i],cali_ordr[1:i] + 20)
  cat(c(formatC(i, width = 2, format = "d", flag = "0"),input,'\n'))
}
sink()

# rep44_order 
tmp <- rep44_match[order(rep44_match[,3],decreasing = T),]
sink('/cloud/project/hpc_pre/rep_hpc_com.txt')
for(j in 1:22) {
  input <- c(tmp[1:j,1],tmp[1:j,2])
  cat(c(formatC(j, width = 2, format = "d", flag = "0"),input,'\n'))
}
sink()

#expi_yps
tmp <- expi_yps_match[order(expi_yps_match[,3],decreasing = T),]
sink('/cloud/project/hpc_pre/expi_yps_hpc_com.txt')
for(j in 1:10) {
  input <- c(tmp[1:j,1],tmp[1:j,2])
  cat(c(formatC(j, width = 2, format = "d", flag = "0"),input,'\n'))
}
sink()

#expi_rm
tmp <- expi_rm_match[order(expi_rm_match[,3],decreasing = T),]
sink('/cloud/project/hpc_pre/expi_rm_hpc_com.txt')
for(j in 1:10) {
  input <- c(tmp[1:j,1],tmp[1:j,2])
  cat(c(formatC(j, width = 2, format = "d", flag = "0"),input,'\n'))
}
sink()

# expi
expi_ordr <- expi_cor_match[,1]
sink('/cloud/project/hpc_pre/expi_true_hpc_com.txt')
for(i in 1:20) {
  input <- c(expi_ordr[1:i],expi_ordr[1:i] + 20)
  cat(c(formatC(i, width = 2, format = "d", flag = "0"),input,'\n'))
}
sink()
