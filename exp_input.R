
########################## get yps genome mapping expression #############
setwd("/cloud/project//bowtie2_yps128_5/")
fileNum <- sapply(strsplit(x=list.files(path="."),split = "_"),"[",1)
fileLab <- rep(c("A","H","C"),times=10)
smpID <- paste0(fileNum,fileLab)
file <- paste0(fileNum,"_", fileLab,".pilecount")

exp_yps <- read.table(file[1],header=T,stringsAsFactors = F)
colnames(exp_yps)[-(1:2)] <- paste0("y",smpID[1],"_",c("rf","rc","af","ac")) # coverage,ref,refcount, alter fer, alter count


for (i in 2:30) {
  print(file[i])
  print(smpID[i])
  tmp <- read.table(file[i],header=T,stringsAsFactors = F)
  colnames(tmp)[-c(1,2)] <- paste0("y",smpID[i],"_",c("rf","rc","af","ac"))
  exp_yps <- merge.data.frame(exp_yps,tmp,by.x=c(1,2),by.y=c(1,2),all=F)  # if all=T, then has 87381 snps here; if all=F:42647
}


colnames(exp_yps)[c(1,2)] <- c("ypsChrom","ypsPosit")
###################### get rm genome mapping expression #############
setwd("/cloud/project//bowtie2_rm11_B/")
fileNum <- sapply(strsplit(x=list.files(path="."),split = "_"),"[",1)
fileLab <- rep(c("A","H","C"),times=10)
smpID <- paste0(fileNum,fileLab)
file <- paste0(fileNum,"_", fileLab,".pilecount")

exp_rm <- read.table(file[1],header=T,stringsAsFactors = F)
colnames(exp_rm)[-(1:2)] <- paste0("r",smpID[1],"_",c("rf","rc","af","ac")) # coverage,ref,refcount, alter fer, alter count

for (i in 2:30) {
  print(file[i])
  print(smpID[i])
  tmp <- read.table(file[i],header=T,stringsAsFactors = F)
  colnames(tmp)[-c(1,2)] <- paste0("r",smpID[i],"_",c("rf","rc","af","ac"))
  exp_rm <- merge.data.frame(exp_rm,tmp,by.x=c(1,2),by.y=c(1,2),all=F)  # if all=T, then has 87381 snps here; if all=F
}

colnames(exp_rm)[c(1,2)] <- c("rmChrom","rmPosit")

#### merge yps_rm data based on posion map ##########
setwd('/cloud/project/otherRaw/')
load("posmap_52filter.RData")

mer1 <- merge.data.frame(posmap_52filter,exp_yps,by = c(1,2),sort = F,all = F)
mer2 <- merge.data.frame(mer1,exp_rm,by.x = c(3,4),by.y = c(1,2), sort = F, all = F )
expf <- mer2[,c(3,4,1,2,5:249)] # change order of columns 
#rm(list=ls()[-c(1,2,3)])
j <- sapply(expf,is.factor)
expf[j] <- lapply(expf[j],as.character)


