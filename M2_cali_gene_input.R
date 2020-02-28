setwd("~/cloud/project/M2_cali_gene_yrGroup")
fileNum <- sapply(strsplit(x=list.files(path="./"),split = "_"),"[",2)
fileLab <- rep(c("A","H"),times=10)
smpID <- paste0(fileNum,fileLab)
file <- paste0("D_",fileNum,"_", fileLab,".yr.gCount")

cali_gene <- read.table(file[1],header=F)
colnames(cali_gene)[-1] <- paste0(smpID[1],"_",c("yGc","rGc","comGc")) # coverage,ref,refcount, alter fer, alter count



for (i in 2:20) {
    print(file[i])
    print(smpID[i])
    tmp <- read.table(file[i],header=F)
    colnames(tmp)[-1] <- paste0(smpID[i],"_",c("yGc","rGc","comGc"))
    cali_gene <- merge.data.frame(cali_gene,tmp,by.x=1,by.y=1,all=F,sort=F)  
}

colnames(cali_gene)[1] <- "geneName"
calih <- cali_gene[,grep(".*comGc",names(cali_gene),invert = T,value=T)]

x <- colSums(calih[,grep(".*_yGc",names(calih),value=T)],na.rm=T)
y <- colSums(calih[,grep(".*_rGc",names(calih),value=T)],na.rm=T)
xy <- data.frame(x,y)
ggplot(xy,aes(x=x,y=y)) + geom_point() + labs(title="DNA read counts" , x = "total DNA read counts mapped to Yps128 allele" , y = "total DNA counts mapped to Rm11-1a allele" ) + scale_y_continuous(labels = scales::scientific) + scale_x_continuous(labels = scales::scientific) + geom_abline(slope = 1,intercept = 0) + theme(text=element_text(size=15))

calih <- calih %>% arrange(rowSums(.[-1])) # sort by average read depth 
calih <- calih[,c("geneName",grep(".*yGc",names(calih),value=T),grep(".*rGc",names(calih),value=T))]
# some genes because of python script has "-1" in them. so change them to zero 
calih <- cbind(geneName=calih[,1],calih[,-1] + ( calih[,-1] < 0 ) * -(calih[,-1]))
