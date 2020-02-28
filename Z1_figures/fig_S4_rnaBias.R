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
    geom_point(size=2.2) + 
    labs(title="Expression read counts of 20 replicates using two \n references without suspected positions" , x = "total read counts mapped to YPS128 allele \n (million counts)" , y = "total read counts mapped to RM11-1a allele \n (million counts)" ) + 
    coord_cartesian(ylim=c(2.5*10^6,8*10^6),xlim=c(2.5*10^6,8*10^6)) +
    scale_y_continuous(labels=seq(2.5,8,by=0.5),breaks=seq(2.5*10^6,8*10^6,by=0.5*10^6)) +
    #scale_x_continuous(labels = scales::scientific) +
    scale_x_continuous(labels=seq(2.5,8,by=0.5),breaks=seq(2.5*10^6,8*10^6,by=0.5*10^6)) + 
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

ggsave("~/cloud/project/figures/fig_S4_rnabias.pdf", device = "pdf", width = 16, height = 12, units = "cm",dpi=300)

exph <- exph %>% arrange(rowSums(.[-1])) # sort by average read depth 
exph <- exph[,c("geneName",grep(".*yGc",names(exph),value=T),grep(".*rGc",names(exph),value=T))]
# some genes because of python script has "-1" in them. so change them to zero 
exph <- cbind(geneName=exph[,1],exph[,-1] + ( exph[,-1] < 0 ) * -(exph[,-1]))
