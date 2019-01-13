source("/cloud/project/otherRaw/func.R")

cali_d <- colSums(cali[,grep("^y.*[HA]",names(cali),value=T)],na.rm=T)/colSums(cali[,grep("^r.*[HA]",names(cali),value=T)],na.rm=T)

x115 <- read.table('/cloud/project/otherRaw/E115_summary.txt',header=T,stringsAsFactors = F)
x530 <- read.table('/cloud/project/otherRaw/E530_summary.txt',header = T,stringsAsFactors = F)
x <- rbind(x115,x530)
x <- x[x$RG!='unknown',c('RG','assigned')]

y115 <- read.table('/cloud/project/otherRaw/J115_summary.txt',header=T,stringsAsFactors = F)
y530 <- read.table('/cloud/project/otherRaw/J530_summary.txt',header = T,stringsAsFactors = F)
y <- rbind(y115,y530)
y <- y[y$RG!='unknown',c('RG','assigned')]

xy <- merge.data.frame(x,y,by="RG")
xyadd <- xy$assigned.x + xy$assigned.y
exptotalRead_tmp<- data.frame(RG=xy$RG,CTS=xyadd)  # expression total read counts data 
exR <- exptotalRead_tmp[grep(pattern = "C$",invert = T,exptotalRead_tmp$RG),]

smpNm <- c(paste0('y',gsub(pattern = '_',"",exR$RG)), paste0('r',gsub(pattern = '_',"",exR$RG)))
expi_Ct <- c(exR$CTS * cali_d /(1+cali_d),exR$CTS /(1+cali_d))
names(expi_Ct) <- smpNm ######### smp_C valuable

c120 <- read.table('/cloud/project/otherRaw/cali_demx.txt',header=T,stringsAsFactors = F)
z <- z[grep(pattern = "C$",invert = T,z$RG),]
library(dplyr)
z <- arrange(z,RG)
cali_Ct <- round(c(z$assigned * cali_d /(1+cali_d),z$assigned / (1 + cali_d)))
names(cali_Ct) <- smpNm

save(list = c('expi','cali','expi_Ct','cali_Ct','paraGetB','neglhBinomial_cisonly','paraGetBB','neglhbetaBinomial_cisonly','dbetabinom'),file = '/cloud/project/hpc_pre/hpc_need_Jan12.RData')

####### local test of R program ################ 
names(cali_Ct)[1:10]
names(cali_Ct)[21:30]
