load('/cloud/project/otherRaw/rep44.RData')

source("/cloud/project/otherRaw/func.R")

cali_d <- colSums(cal[-wrong_rows,grep("^y.*[HA]_rc",names(cal),value=T)],na.rm=T)/ colSums(cal[-wrong_rows,grep("^r.*[HA]_rc",names(cal),value=T)],na.rm=T)

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
z <- c120[grep(pattern = "C$",invert = T,c120$RG),c("RG","assigned")]
z <- z[z$RG != "unknown",]
library(dplyr)
z <- arrange(z,RG)
cali_Ct <- round(c(z$assigned * cali_d /(1+cali_d),z$assigned / (1 + cali_d)))
names(cali_Ct) <- smpNm

save(list = c('rep44','depth44','expi','cali','expi_Ct','cali_Ct','paraGetB','neglhBinomial_cisonly','paraGetBB','neglhbetaBinomial_cisonly','dbetabinom'),file = '/cloud/project/hpc_pre/hpc_need_Jan14.RData')

####### local test of R program ########################################################### 
names(cali_Ct)[1:10]
names(cali_Ct)[21:30]

paraGetBB(cali[1903,],names(cali_Ct)[1:19],names(cali_Ct)[21:39],cali_Ct)
paraGetB(cali[1903,],names(cali_Ct)[1:19],names(cali_Ct)[21:39],cali_Ct)

paraGetBB(expi[4,],names(expi_Ct)[1:4],names(expi_Ct)[21:24],expi_Ct)
paraGetB(expi[4,],names(expi_Ct)[1:4],names(expi_Ct)[21:24],expi_Ct)

paraGetBB(repfine[4,],names(depth44)[1:4],names(depth44)[21:24],depth44)
paraGetB(repfine[4,],names(depth44)[1:4],names(depth44)[21:24],depth44)

paraGetBB(expi[expi$ypsGene=='APL5',],names(expi_Ct)[c(8,19,5,3,6,4,9,1,2,20)],names(expi_Ct)[c(10,11,17,14,12,13,15,18,16,7)],expi_Ct)
paraGetB(expi[expi$ypsGene=='APL5',],names(expi_Ct)[c(8,19,5,3,6,4,9,1,2,20)],names(expi_Ct)[c(10,11,17,14,12,13,15,18,16,7)],expi_Ct)


# check likelihood surface, around 0 , it's almost flat, so if starts here, optimizer would fail , and confidence interval fail, because of this too. 
for(i in -20:30) {
ans <- neglhbetaBinomial_cisonly(0.14433,i,
    xHy = unlist(expi[expi$ypsGene=='ALG1',c(8,19,5,3,6,4,9,1,2,20)+1]),
    nHy =unlist(expi[expi$ypsGene=='ALG1',c(8,19,5,3,6,4,9,1,2,20)+1] + expi[expi$ypsGene=='ALG1',c(10,11,17,14,12,13,15,18,16,7)+1]),
    C1 = unlist(expi_Ct[c(8,19,5,3,6,4,9,1,2,20)+1]),
    C2 = unlist(expi_Ct[c(8,19,5,3,6,4,9,1,2,20)+1]))
  print(c(i,ans))
}

neglhBinomial_cisonly(0.14433,
  xHy = unlist(expi[expi$ypsGene=='ALG1',c(8,19,5,3,6,4,9,1,2,20)+1]),
nHy =unlist(expi[expi$ypsGene=='ALG1',c(8,19,5,3,6,4,9,1,2,20)+1] + expi[expi$ypsGene=='ALG1',c(10,11,17,14,12,13,15,18,16,7)+1]),
                C1 = unlist(expi_Ct[c(8,19,5,3,6,4,9,1,2,20)+1]),
                C2 = unlist(expi_Ct[c(8,19,5,3,6,4,9,1,2,20)+1]))
