######### 1 : read in from yps128_5 cali #############
#setwd("/cloud/project/cali_yps128_5/")
setwd("~/cloud/project/F1_snpfilter/cali_yps_5/")
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
#setwd("~/cloud/project/cali_rm11_B/")
setwd("~/cloud/project/F1_snpfilter/cali_rm_B/")
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
setwd('~/cloud/project/F0_otherRaw')
load("posmap52.RData")
mer1 <- merge.data.frame(posmap52,cal_yps,by = c(1,2),sort = F,all.x=T)
mer2 <- merge.data.frame(mer1,cal_rm,by.x = c(3,4),by.y = c(1,2), sort = F, all.x=T)
cal <- mer2[,c(3,4,1,2,5:249)]
j <- sapply(cal,is.factor)
cal[j] <- lapply(cal[j],as.character)

yr <- colSums(cal[,grep("^y.*[HA]_rc",names(cal),value=T)],na.rm=T)
ya <- colSums(cal[,grep("^y.*[HA]_ac",names(cal),value=T)],na.rm=T)
yra <- data.frame(x=yr,y=ya, type="YPS128")

rr <- colSums(cal[,grep("^r.*[HA]_rc",names(cal),value=T)],na.rm=T)
ra <- colSums(cal[,grep("^r.*[HA]_ac",names(cal),value=T)],na.rm=T)
rra <- data.frame(x=ra,y=rr, type="RM11-1a")

df = rbind(yra,rra)

ggplot(df) + 
  geom_point(aes(x=x,y=y,fill=type),size=2.2,color='black', pch=21) + 
  scale_fill_manual(values=c('white', 'black')) + 
  #geom_point(data=rra,aes(x=ra,y=rr), color='black') + 
  geom_abline(slope = 1,intercept = 0) +
  labs(title="Mapping bias of 20 replicates using one reference" , x = "total read counts mapped to YPS128 allele \n(million counts)" , y = "total read counts mapped to RM11-1a allele \n(million counts)", fill='genome used as reference' ) + 
  #scale_y_continuous(labels = scales::scientific) +
  scale_y_continuous(labels=seq(0,2.5,by=0.25),breaks=seq(0,2.5*10^6,by=0.25*10^6)) +
  #scale_x_continuous(labels = scales::scientific) +
  scale_x_continuous(labels=seq(0,2.5,by=0.25),breaks=seq(0,2.5*10^6,by=0.25*10^6)) + 
  geom_abline(slope = 1,intercept = 0) + 
  theme_bw() +
  theme(axis.text.x=element_text(angle=0, hjust=1)) +
  theme(
    panel.grid = element_line(color="grey85"),
    legend.position = c(0.8, 0.3),
    legend.background = element_rect( size=0.5, linetype="solid",colour ="black"),
    legend.title = element_text(face='bold'),
    legend.text = element_text(size=12,face="plain"),
    plot.title = element_text(size = 16, face = "bold",hjust=.5),
    axis.text.x = element_text(colour="grey20",size=13,hjust=.5,vjust=.5,face="plain"),
    axis.text.y = element_text(colour="grey20",size=13,hjust=.5,vjust=.5,face="plain"),
    axis.title.x = element_text(colour="grey20",size=13,hjust=.5,vjust=0,face="bold"), 
    axis.title.y = element_text(colour="grey20",size=13,hjust=.5,vjust=.5,face="bold"))

ggsave("~/cloud/project/figures/fig_S1_dnabias.pdf", device = "pdf", width = 16, height = 12, units = "cm",dpi=300)
