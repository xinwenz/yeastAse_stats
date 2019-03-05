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



####### rep44_ fixed_hybrid 22######4#####
hymatx1 <- matrix(sample((1:44),size=44,replace = F),ncol=2)
hymatx2 <- matrix(sample((1:44),size=44,replace = F),ncol=2)
hymatx3 <- matrix(sample((1:44),size=44,replace = F),ncol=2)
hymatx4 <- matrix(sample((1:44),size=44,replace = F),ncol=2)

sink('~/cloud/project/hpc_pre/rep_null_hpc_com_hym1.txt')
i <- 1 
for(siz in 1:22) 
  for(t in 1:150) {
    set.seed(siz+t*2)
    choRep <- sample(1:22,size=siz,replace = F)
    input <- c(hymatx1[choRep,1],hymatx1[choRep,2])
    cat(c(formatC(i, width = 4, format = "d", flag = "0"),input,'\n'))
    i <- i + 1 
}
sink()
sink('~/cloud/project/hpc_pre/rep_null_hpc_com_hym1.txt')
i <- 1 
for(siz in 1:22) 
  for(t in 1:150) {
    set.seed(siz+t)
    choRep <- sample(1:22,size=siz,replace = F)
    input <- c(hymatx1[choRep,1],hymatx1[choRep,2])
    cat(c(formatC(i, width = 4, format = "d", flag = "0"),input,'\n'))
    i <- i + 1 
  }
sink()
sink('~/cloud/project/hpc_pre/rep_null_hpc_com_hym2.txt')
i <- 1 
for(siz in 1:22) 
  for(t in 1:150) {
    set.seed(siz+t*2)
    choRep <- sample(1:22,size=siz,replace = F)
    input <- c(hymatx2[choRep,1],hymatx2[choRep,2])
    cat(c(formatC(i, width = 4, format = "d", flag = "0"),input,'\n'))
    i <- i + 1 
  }
sink()
sink('~/cloud/project/hpc_pre/rep_null_hpc_com_hym3.txt')
i <- 1 
for(siz in 1:22) 
  for(t in 1:150) {
    set.seed(siz+t*3)
    choRep <- sample(1:22,size=siz,replace = F)
    input <- c(hymatx3[choRep,1],hymatx3[choRep,2])
    cat(c(formatC(i, width = 4, format = "d", flag = "0"),input,'\n'))
    i <- i + 1 
  }
sink()
sink('~/cloud/project/hpc_pre/rep_null_hpc_com_hym4.txt')
i <- 1 
for(siz in 1:22) 
  for(t in 1:150) {
    set.seed(siz+t*4)
    choRep <- sample(1:22,size=siz,replace = F)
    input <- c(hymatx4[choRep,1],hymatx4[choRep,2])
    cat(c(formatC(i, width = 4, format = "d", flag = "0"),input,'\n'))
    i <- i + 1 
  }
sink()

################exph_yps###########



#############################################################################
hymatx1 <- matrix(sample((1:20)[-c(7,19)],size=18,replace = F),ncol=2)
hymatx2 <- matrix(sample((1:20)[-c(7,19)],size=18,replace = F),ncol=2)
hymatx3 <- matrix(sample((1:20)[-c(7,19)],size=18,replace = F),ncol=2)
hymatx4 <- matrix(sample((1:20)[-c(7,19)],size=18,replace = F),ncol=2)

sink('~/cloud/project/hpc_pre/exph_ypsRM14A9A_hpc_com_hym1.txt')
i <- 1 
for(siz in 1:9) 
  for(t in 1:100) {
    set.seed(siz + t*1)
    choRep <- sample(1:9,size=siz,replace = F)
    input <- c(hymatx1[choRep,1],hymatx1[choRep,2])
    cat(c(formatC(i, width = 3, format = "d", flag = "0"),input,'\n'))
    i <- i + 1 
  }
sink()

sink('~/cloud/project/hpc_pre/exph_ypsRM14A9A_hpc_com_hym2.txt')
i <- 1 
for(siz in 1:9) 
  for(t in 1:100) {
    set.seed(siz+ t*2)
    choRep <- sample(1:9,size=siz,replace = F)
    input <- c(hymatx2[choRep,1],hymatx2[choRep,2])
    cat(c(formatC(i, width = 3, format = "d", flag = "0"),input,'\n'))
    i <- i + 1 
  }
sink()

sink('~/cloud/project/hpc_pre/exph_ypsRM14A9A_hpc_com_hym3.txt')
i <- 1 
for(siz in 1:9) 
  for(t in 1:100) {
    set.seed(siz+t*3)
    choRep <- sample(1:9,size=siz,replace = F)
    input <- c(hymatx3[choRep,1],hymatx3[choRep,2])
    cat(c(formatC(i, width = 3, format = "d", flag = "0"),input,'\n'))
    i <- i + 1 
  }
sink()

sink('~/cloud/project/hpc_pre/exph_ypsRM14A9A_hpc_com_hym4.txt')
i <- 1 
for(siz in 1:9) 
  for(t in 1:100) {
    set.seed(siz + t*4)
    choRep <- sample(1:9,size=siz,replace = F)
    input <- c(hymatx4[choRep,1],hymatx4[choRep,2])
    cat(c(formatC(i, width = 3, format = "d", flag = "0"),input,'\n'))
    i <- i + 1 
  }
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


library(gtools)
null44_scr <- permutations(n = 22,r=2)
sink('~/cloud/project/hpc_pre/rep_null_hpc_com_nostr.txt')
i <- 1 
for(siz in 1:35) 
  for(t in 1:500) {
    set.seed(siz+t*2)
    choRep <- sample(1:nrow(null44_scr),size=siz,replace = F)
    input <- c(null44_scr[choRep,1],null44_scr[choRep,2])
    cat(c(formatC(i, width = 5, format = "d", flag = "0"),input,'\n'))
    i <- i + 1 
  }
sink()

yps_scr <- permutations(n = 18,r=2)
sink('~/cloud/project/hpc_pre/exph_ypsRM14A9A_hpc_com_nostr.txt')
i <- 1 
for(siz in 1:20) 
  for(t in 1:150) {
    set.seed(siz+t*2)
    choRep <- sample(1:nrow(yps_scr),size=siz,replace = F)
    input <- c(yps_scr[choRep,1],yps_scr[choRep,2])
    cat(c(formatC(i, width = 4, format = "d", flag = "0"),input,'\n'))
    i <- i + 1 
  }
sink()



sink('~/cloud/project/hpc_pre/exph_trueRM14A9A_hpc_com.txt')
i <- 1 
tmplt <- (1:20)[-c(7,19)]
for(siz in 1:18) {
  for(t in 1:150) {
    choRep <- sample(tmplt,size=siz,replace = F)
    input <- c(choRep,choRep + 20)
    cat(c(formatC(i, width = 4, format = "d", flag = "0"),input,'\n'))
    i <- i + 1 
  }
}
sink()

