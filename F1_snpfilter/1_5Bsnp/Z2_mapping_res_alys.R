setwd("~/cloud/project/snpfilter/1_5Bsnp/mapping_res/")
fileNum <- sapply(strsplit(x=list.files(path="./"),split = "_"),"[",2)
fileLab <- rep(c("A","H"),times=10)
smpID <- paste0(fileNum,fileLab)
file <- paste0("D_",fileNum,"_", fileLab,".yr.gCount")

cali_1 <- read.table(file[1],header=F)
colnames(cali_1)[-1] <- paste0(smpID[1],"_",c("yGc","rGc","comGc")) # coverage,ref,refcount, alter fer, alter count



for (i in 2:20) {
    print(file[i])
    print(smpID[i])
    tmp <- read.table(file[i],header=F)
    colnames(tmp)[-1] <- paste0(smpID[i],"_",c("yGc","rGc","comGc"))
    cali_1 <- merge.data.frame(cali_1,tmp,by.x=1,by.y=1,all=T,sort=F)  
}

colnames(cali_1)[1] <- "group"
calh_1 <- cali_1[,grep(".*comGc",names(cali_1),invert = T,value=T)]

x <- colSums(calh_1[,grep(".*_yGc",names(calh_1),value=T)],na.rm=T)
y <- colSums(calh_1[,grep(".*_rGc",names(calh_1),value=T)],na.rm=T)
xy <- data.frame(x,y)

ggplot(xy,aes(x=x,y=y)) + 
    geom_point(size=3) +
    labs(title="Mapping of 20 replicates using two references \n with all positions " , x = "total read counts mapped to Yps128 allele (million counts)" , y = "total read counts mapped to Rm11-1a allele (million counts)" ) + 
    #scale_y_continuous(labels = scales::scientific) +
   coord_cartesian(ylim=c(0,1.2*10^6),xlim=c(0,1.2*10^6)) +
    scale_y_continuous(labels=seq(0,1.2,by=0.2),breaks=seq(0,1.2*10^6,by=0.2*10^6)) +
    #scale_x_continuous(labels = scales::scientific) +
    scale_x_continuous(labels=seq(0,1.2,by=0.2),breaks=seq(0,1.2*10^6,by=0.2*10^6)) + 
    geom_abline(slope = 1,intercept = 0) + 
    theme_bw() +
    theme(axis.text.x=element_text(angle=0, hjust=1)) +
    theme(
        panel.grid = element_line(color="grey85"),  
        legend.position = c(0.8, 0.3),
        legend.background = element_rect( size=0.5, linetype="solid",colour ="black"),
        legend.title = element_text(face='bold'),
        legend.text = element_text(size=14,face="plain"),
        plot.title = element_text(size = 18, face = "bold",hjust=.5),
        axis.text.x = element_text(colour="grey20",size=13,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=13,hjust=.5,vjust=.5,face="plain"),
        axis.title.x = element_text(colour="grey20",size=15,hjust=.5,vjust=0,face="bold"), 
        axis.title.y = element_text(colour="grey20",size=15,hjust=.5,vjust=.5,face="bold"))

