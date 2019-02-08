setwd("~/cloud/project/cali_yrGroup")
fileNum <- sapply(strsplit(x=list.files(path="./"),split = "_"),"[",2)
fileLab <- rep(c("A","H"),times=10)
smpID <- paste0(fileNum,fileLab)
file <- paste0("D_",fileNum,"_", fileLab,".yr.gCount")

cali <- read.table(file[1],header=F)
colnames(cali)[-1] <- paste0(smpID[1],"_",c("yGc","rGc","comGc")) # coverage,ref,refcount, alter fer, alter count



for (i in 2:20) {
  print(file[i])
  print(smpID[i])
  tmp <- read.table(file[i],header=F)
  colnames(tmp)[-1] <- paste0(smpID[i],"_",c("yGc","rGc","comGc"))
  cali <- merge.data.frame(cali,tmp,by.x=1,by.y=1,all=T,sort=F)  
}

colnames(cali)[1] <- "group"
calh <- cali[,grep(".*comGc",names(cali),invert = T,value=T)]

x <- colSums(calh[,grep(".*_yGc",names(calh),value=T)],na.rm=T)
y <- colSums(calh[,grep(".*_rGc",names(calh),value=T)],na.rm=T)
xy <- data.frame(x,y)
ggplot(xy,aes(x=x,y=y)) + geom_point() + labs(title="inlcude co-culture" , x = "total read counts mapped to Yps128 allele" , y = "total read counts mapped to Rm11-1a allele" ) + scale_y_continuous(labels = scales::scientific) + scale_x_continuous(labels = scales::scientific) + geom_abline(slope = 1,intercept = 0) + theme(text=element_text(size=15))

