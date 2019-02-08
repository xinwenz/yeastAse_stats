load('/cloud/project/otherRaw/run_hpc_Cor_seq.RData')
#cali_C order:####### 
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


########## exph ##########
sink('~/cloud/project/hpc_pre/exph_true_hpc_com.txt')
i <- 1 
while(i <= 20) {
  input <- c(i,i+20)
  cat(c(formatC(i, width = 3, format = "d", flag = "0"),input,'\n'))
  i <- i + 1
}

for(siz in 2:18) 
  for(t in 1:20) {
    choRep <- sample(1:20,size=siz,replace = F)
    input <- c(choRep,choRep + 20)
    cat(c(formatC(i, width = 3, format = "d", flag = "0"),input,'\n'))
    i <- i + 1 
}

for(k in 1:20) {
  choRep <- (1:20)[-k]
  input <- c(choRep,choRep + 20 )
  cat(c(formatC(i, width = 3, format = "d", flag = "0"),input,'\n'))
  i <- i + 1 
}

choRep <- (1:20)
input <- c(choRep,choRep + 20 )
cat(c(formatC(i, width = 3, format = "d", flag = "0"),input,'\n'))

sink()



####### rep44_ fixed_hybrid 22##############
sink('~/cloud/project/hpc_pre/rep_null_hpc_com.txt')
i <- 1 
while(i <= 22) {
  input <- c(i,i+22)
  cat(c(formatC(i, width = 3, format = "d", flag = "0"),input,'\n'))
  i <- i + 1
}

for(siz in 2:20) 
  for(t in 1:22) {
    choRep <- sample(1:22,size=siz,replace = F)
    input <- c(choRep,choRep + 22)
    cat(c(formatC(i, width = 3, format = "d", flag = "0"),input,'\n'))
    i <- i + 1 
}

for(k in 1:22) {
  choRep <- (1:22)[-k]
  input <- c(choRep,choRep + 22 )
  cat(c(formatC(i, width = 3, format = "d", flag = "0"),input,'\n'))
  i <- i + 1 
}

choRep <- (1:22)
input <- c(choRep,choRep + 22 )
cat(c(formatC(i, width = 3, format = "d", flag = "0"),input,'\n'))

sink()

################exph_yps###########
sink('~/cloud/project/hpc_pre/exph_yps_hpc_com.txt')
i <- 1 
while(i <= 10) {
  input <- c(i,i+10)
  cat(c(formatC(i, width = 3, format = "d", flag = "0"),input,'\n'))
  i <- i + 1
}

for(siz in 2:8) 
  for(t in 1:10) {
    choRep <- sample(1:10,size=siz,replace = F)
    input <- c(choRep,choRep + 10)
    cat(c(formatC(i, width = 3, format = "d", flag = "0"),input,'\n'))
    i <- i + 1 
}

for(k in 1:10) {
  choRep <- (1:10)[-k]
  input <- c(choRep,choRep + 10 )
  cat(c(formatC(i, width = 3, format = "d", flag = "0"),input,'\n'))
  i <- i + 1 
}

cat(c(formatC(i, width = 3, format = "d", flag = "0"),1:20,'\n'))

sink()

sink('~/cloud/project/hpc_pre/exph_ypsRM14A_hpc_com.txt')
tmplt <- (1:10)[-7]
i <- 1 
while(i <= 9) {
  input <- c(tmplt[i],tmplt[i]+10 )
  cat(c(formatC(i, width = 3, format = "d", flag = "0"),input,'\n'))
  i <- i + 1
}

for(siz in 2:7) 
  for(t in 1:9) {
    choRep <- sample(tmplt,size=siz,replace = F)
    input <- c(choRep,choRep + 10)
    cat(c(formatC(i, width = 3, format = "d", flag = "0"),input,'\n'))
    i <- i + 1 
  }

for(k in 1:9) {
  choRep <- tmplt[-k]
  input <- c(choRep,choRep + 10 )
  cat(c(formatC(i, width = 3, format = "d", flag = "0"),input,'\n'))
  i <- i + 1 
}

cat(c(formatC(i, width = 3, format = "d", flag = "0"),tmplt,tmplt+10,'\n'))

sink()


sink('~/cloud/project/hpc_pre/exph_ypsRM14A9A_hpc_com.txt')
tmplt <- (1:10)[-c(7,9)]
i <- 1 
while(i <= 8) {
  input <- c(tmplt[i],tmplt[i]+10 )
  cat(c(formatC(i, width = 3, format = "d", flag = "0"),input,'\n'))
  i <- i + 1
}

for(siz in 2:6) 
  for(t in 1:8) {
    choRep <- sample(tmplt,size=siz,replace = F)
    input <- c(choRep,choRep + 10)
    cat(c(formatC(i, width = 3, format = "d", flag = "0"),input,'\n'))
    i <- i + 1 
  }

for(k in 1:8) {
  choRep <- tmplt[-k]
  input <- c(choRep,choRep + 10 )
  cat(c(formatC(i, width = 3, format = "d", flag = "0"),input,'\n'))
  i <- i + 1 
}

cat(c(formatC(i, width = 3, format = "d", flag = "0"),tmplt,tmplt+10,'\n'))

sink()

#######################  
sink('~/cloud/project/hpc_pre/exph_trueRM14A9A_hpc_com.txt')
tmplt <- (1:20)[-c(7,19)]
i <- 1 
while(i <= 18) {
  input <- c(tmplt[i],tmplt[i]+20)
  cat(c(formatC(i, width = 3, format = "d", flag = "0"),input,'\n'))
  i <- i + 1
}

for(siz in 2:16) 
  for(t in 1:18) {
    choRep <- sample(tmplt,size=siz,replace = F)
    input <- c(choRep,choRep + 20)
    cat(c(formatC(i, width = 3, format = "d", flag = "0"),input,'\n'))
    i <- i + 1 
  }

for(k in 1:18) {
  choRep <- tmplt[-k]
  input <- c(choRep,choRep + 20 )
  cat(c(formatC(i, width = 3, format = "d", flag = "0"),input,'\n'))
  i <- i + 1 
}

cat(c(formatC(i, width = 3, format = "d", flag = "0"),tmplt,tmplt+20,'\n'))

sink()