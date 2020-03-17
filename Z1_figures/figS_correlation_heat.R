library(reshape2)
library(ggplot2)
library(viridis)

setwd("~/cloud/project/M3_expi_yrGroup")
fileNum <- sapply(strsplit(x=list.files(path="./"),split = "_"),"[",1)
fileLab <- rep(c("A","H"),times=10)
smpID <- paste0(fileNum,fileLab)
file <- paste0(fileNum,"_", fileLab,".yr.geneCount")

expi <- read.table(file[1],header=F)
colnames(expi)[-1] <- paste0(smpID[1],"_",c("yGc","rGc","comGc")) 

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
ggplot(xy,aes(x=x,y=y)) + geom_point() + labs(title="Expression read counts" , x = "total Expression read counts mapped to Yps128 allele" , y = "total read counts mapped to Rm11-1a allele" ) + scale_y_continuous(labels = scales::scientific) + scale_x_continuous(labels = scales::scientific) + geom_abline(slope = 1,intercept = 0) + theme(text=element_text(size=15))

exph <- exph %>% arrange(rowSums(.[-1])) # sort by average read depth 
exph <- exph[,c("geneName",grep(".*yGc",names(exph),value=T),grep(".*rGc",names(exph),value=T))]
# some genes because of python script has "-1" in them. so change them to zero 
exph <- cbind(geneName=exph[,1],exph[,-1] + ( exph[,-1] < 0 ) * -(exph[,-1]))
expfc <- exph
colnames(expfc)[-1] <- paste0(c(rep("y",20),rep("r",20)), gsub("_.*",'',colnames(exph)[-1]))

melted_cormat <- melt(cor(expfc[,-1],method='pearson'))
ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
    scale_fill_viridis(alpha = 1, begin = 0, end = 1, direction = 1,discrete = FALSE, option = "B") +
    geom_tile() + labs(title=" pearson correlation of expression profiles",x='40 expression profiles', y= '40 expression profiles' ) + 
    coord_fixed(ratio = 1) +
    theme_bw() +
    theme(axis.text.x=element_text(angle=90, hjust=1)) +
    theme(
        panel.grid = element_line(color="grey85"),  
        plot.title = element_text(size = 18, face = "bold",hjust=.5),
        axis.text.x = element_text(colour="grey20",size=8,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=8,hjust=.5,vjust=.5,face="plain"),
        axis.title.x = element_text(colour="grey20",size=13,hjust=.5,vjust=0,face="bold"), 
        axis.title.y = element_text(colour="grey20",size=13,hjust=.5,vjust=.5,face="bold"))

ggsave("~/cloud/project/figures/fig_S5_correlation.pdf", device = "pdf", width = 16, height = 15, units = "cm",dpi=300)

