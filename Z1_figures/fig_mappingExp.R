setwd("~/cloud/project/expi_yrGroup")
fileNum <- sapply(strsplit(x=list.files(path="./"),split = "_"),"[",1)
fileLab <- rep(c("A","H"),times=10)
smpID <- paste0(fileNum,fileLab)
file <- paste0(fileNum,"_", fileLab,".yr.geneCount")

expi <- read.table(file[1],header=F)
colnames(expi)[-1] <- paste0(smpID[1],"_",c("yGc","rGc","comGc")) # coverage,ref,refcount, alter fer, alter count



for (i in 2:20) {
    print(file[i])
    print(smpID[i])
    tmp <- read.table(file[i],header=F)
    colnames(tmp)[-1] <- paste0(smpID[i],"_",c("yGc","rGc","comGc"))
    expi <- merge.data.frame(expi,tmp,by.x=1,by.y=1,all=F,sort=F)  
}

colnames(expi)[1] <- "geneName"
exph <- expi[,grep(".*comGc",names(expi),invert = T,value=T)]

x <- colSums(exph[,grep(".*_yGc",names(exph),value=T)],na.rm=T)
y <- colSums(exph[,grep(".*_rGc",names(exph),value=T)],na.rm=T)
xy <- data.frame(x,y)

ggplot(xy,aes(x=x,y=y)) + 
    geom_point(size=3) +
    labs(title="Expression read counts of 20 replicates using two references \n without suspected positions " , x = "total read counts mapped to Yps128 allele (million counts)" , y = "total read counts mapped to Rm11-1a allele (million counts)" ) + 
    #scale_y_continuous(labels = scales::scientific) +
    scale_y_continuous(labels=seq(0,8,by=0.5),breaks=seq(0,8*10^6,by=0.5*10^6)) +
    #scale_x_continuous(labels = scales::scientific) +
    scale_x_continuous(labels=seq(0,8,by=0.5),breaks=seq(0,8*10^6,by=0.5*10^6)) + 
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
