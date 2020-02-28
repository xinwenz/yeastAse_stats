setwd("~/cloud/project/F1_snpfilter/7_5Bsnp_exp/mapping_res/")
fileNum <- sapply(strsplit(x=list.files(path="./"),split = "_"),"[",1)
fileLab <- rep(c("A","H"),times=10)
smpID <- paste0(fileNum,fileLab)
file <- paste0(fileNum,"_", fileLab,".yr.geneCount")

exp_da <- read.table(file[1],header=F)
colnames(exp_da)[-1] <- paste0(smpID[1],"_",c("yGc","rGc","comGc")) # coverage,ref,refcount, alter fer, alter count



for (i in 2:20) {
    print(file[i])
    print(smpID[i])
    tmp <- read.table(file[i],header=F)
    colnames(tmp)[-1] <- paste0(smpID[i],"_",c("yGc","rGc","comGc"))
    exp_da <- merge.data.frame(exp_da,tmp,by.x=1,by.y=1,all=F,sort=F)  
}

colnames(exp_da)[1] <- "geneName"
exp_db <- exp_da[,grep(".*comGc",names(exp_da),invert = T,value=T)]

x <- colSums(exp_db[,grep(".*_yGc",names(exp_db),value=T)],na.rm=T)
y <- colSums(exp_db[,grep(".*_rGc",names(exp_db),value=T)],na.rm=T)
xy <- data.frame(x,y)
ggplot(xy,aes(x=x,y=y)) + geom_point() + labs(title="Expression read counts" , x = "total Expression read counts mapped to Yps128 allele" , y = "total read counts mapped to Rm11-1a allele" ) + scale_y_continuous(labels = scales::scientific) + scale_x_continuous(labels = scales::scientific) + geom_abline(slope = 1,intercept = 0) + theme(text=element_text(size=15))

exp_db <- exp_db %>% arrange(rowSums(.[-1])) # sort by average read depth 
expdb <- exp_db[,c("geneName",grep(".*yGc",names(exph),value=T),grep(".*rGc",names(exph),value=T))]
# some genes because of python script has "-1" in them. so change them to zero 
expdc <- cbind(geneName=expdb[,1],expdb[,-1] + ( expdb[,-1] < 0 ) * -(expdb[,-1]))
save(expdc, file="~/cloud/project/F1_snpfilter/7_5Bsnp_exp/expdc.RData")
