setwd("/cloud/project")

######### 1 : read in from yps128_5 cali #############
#setwd("/cloud/project/cali_yps128_5/")
setwd("/cloud/project/cali_yps_5/")
fileNum <- sapply(strsplit(x=list.files(path="./"),split = "_"),"[",2)
fileLab <- rep(c("A","H","C"),times=10)
smpID <- paste0(fileNum,fileLab)
file <- paste0("D_",fileNum,"_", fileLab,".pilecount")

cal_yps <- read.table(file[1],header=T)
colnames(cal_yps)[-(1:2)] <- paste0("y",smpID[1],"_",c("rf","rc","af","ac")) # coverage,ref,refcount, alter fer, alter count


for (i in 2:30) {
  print(file[i])
  print(smpID[i])
  tmp <- read.table(file[i],header=T)
  colnames(tmp)[-c(1,2)] <- paste0("y",smpID[i],"_",c("rf","rc","af","ac"))
  cal_yps <- merge.data.frame(cal_yps,tmp,by.x=c(1,2),by.y=c(1,2),all=T,sort=F)  # if all=T, then has 87381 snps here; if all=F:42647
}

colnames(cal_yps)[c(1,2)] <- c("ypsChrom","ypsPosit")


######### 2 : read in from  rm11_B cali #############
#setwd("/cloud/project/cali_rm11_B/")
setwd("/cloud/project/cali_rm11_B/")
fileNum <- sapply(strsplit(x=list.files(path="."),split = "_"),"[",2)
fileLab <- rep(c("A","H","C"),times=10)
smpID <- paste0(fileNum,fileLab)
file <- paste0('D_',fileNum,"_", fileLab,".pilecount")

cal_rm <- read.table(file[1],header=T)
colnames(cal_rm)[-(1:2)] <- paste0("r",smpID[1],"_",c("rf","rc","af","ac")) # coverage,ref,refcount, alter fer, alter count

for (i in 2:30) {
  print(file[i])
  print(smpID[i])
  tmp <- read.table(file[i],header=T)
  colnames(tmp)[-c(1,2)] <- paste0("r",smpID[i],"_",c("rf","rc","af","ac"))
  cal_rm <- merge.data.frame(cal_rm,tmp,by.x=c(1,2),by.y=c(1,2),all=T,sort=F)  # if all=T, then has 87381 snps here; if all=F
}

colnames(cal_rm)[c(1,2)] <- c("rmChrom","rmPosit")


######### 3 : merge to get cal (data frame) ##########
setwd('/cloud/project/otherRaw')
load("posmap52.RData")
mer1 <- merge.data.frame(posmap52,cal_yps,by = c(1,2),sort = F,all.x=T)
mer2 <- merge.data.frame(mer1,cal_rm,by.x = c(3,4),by.y = c(1,2), sort = F, all.x=T)
cal <- mer2[,c(3,4,1,2,5:249)]
j <- sapply(cal,is.factor)
cal[j] <- lapply(cal[j],as.character)

x <- colSums(cal[,grep("^y.*[HAC]_rc",names(cal),value=T)],na.rm=T)
y <- colSums(cal[,grep("^r.*[HAC]_rc",names(cal),value=T)],na.rm=T)
xy <- data.frame(x,y)
ggplot(xy,aes(x=x,y=y)) + geom_point() + labs(title="inlcude co-culture" , x = "total read counts mapped to Yps128 allele" , y = "total read counts mapped to Rm11-1a allele" ) + scale_y_continuous(labels = scales::scientific) + scale_x_continuous(labels = scales::scientific) + geom_abline(slope = 1,intercept = 0) + theme(text=element_text(size=15))

