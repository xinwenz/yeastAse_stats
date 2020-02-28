library(dplyr)

setwd('~/cloud/project/F0_otherRaw')
# remove mummer unmatch snps   
yr_pos <- read.table(header=T,file='yr5B.posmap')
names(yr_pos)[1:6] <- c('ypsChrom','ypsPosit','ypsN','rmChrom','rmPosit','rmN')
ry_pos <- read.table(header=T,file='ryB5.posmap')
names(ry_pos)[1:6] <- c('rmChrom','rmPosit','rmN','ypsChrom','ypsPosit','ypsN')
ry_pos_exchg <- ry_pos[,c(4,5,6,1,2,3,7)]
res <- merge.data.frame(yr_pos,ry_pos_exchg,by.x = c(1,2,4,5),by.y=c(1,2,4,5),all=T,sort=F)
res_good <- res[which(res$drct.x == res$drct.y),] # 86351

#write.table(res_good[,c(1,2)],file="~/cloud/project/snpfilter/yps128_5_snpls_351",row.names = F,quote=F,col.names = F)
#write.table(res_good[,c(3,4)],file="~/cloud/project/snpfilter/rm11_B_snpls_351",row.names = F,quote=F,col.names = F)

### go to HPC , do a samtools count , in order to filter some SNPs without any evidence ### 
### hpc data saved to #####
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
  cal_yps <- merge.data.frame(cal_yps,tmp,by.x=c(1,2),by.y=c(1,2),all=T,sort=F)
}

colnames(cal_yps)[c(1,2)] <- c("ypsChrom","ypsPosit")

#####
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

mer1 <- merge.data.frame(res_good,cal_yps,by = c(1,2),sort = F,all.x=T)
mer2 <- merge.data.frame(mer1,cal_rm,by.x = c(3,4),by.y = c(1,2), sort = F, all.x=T)
cal <- mer2[,c(3,4,1,2,5:250)]
j <- sapply(cal,is.factor)
cal[j] <- lapply(cal[j],as.character)

x <- colSums(cal[,grep("^y.*[HAC]_rc",names(cal),value=T)],na.rm=T)
y <- colSums(cal[,grep("^r.*[HAC]_rc",names(cal),value=T)],na.rm=T)
xy <- data.frame(x,y)

ggplot(xy,aes(x=x,y=y)) + geom_point() + labs(title="inlcude co-culture" , x = "total read counts mapped to Yps128 allele" , y = "total read counts mapped to Rm11-1a allele" ) + scale_y_continuous(labels = scales::scientific) + scale_x_continuous(labels = scales::scientific) + geom_abline(slope = 1,intercept = 0) + theme(text=element_text(size=15))

## use YPS as reference 
x <- colSums(cal[,grep("^y.*[HA]_rc",names(cal),value=T)],na.rm=T)
y <- colSums(cal[,grep("^y.*[HA]_ac",names(cal),value=T)],na.rm=T)
xy1 <- data.frame(x,y)
xy1 <- cbind(xy1, genome="Yps128")

y <- colSums(cal[,grep("^r.*[HA]_rc",names(cal),value=T)],na.rm=T)
x <- colSums(cal[,grep("^r.*[HA]_ac",names(cal),value=T)],na.rm=T)
xy2 <- data.frame(x,y)
xy2 <- cbind(xy2, genome="Rm11-1a")

xyxy <- rbind(xy1,xy2)

ggplot(xyxy,aes(x=x,y=y, fill=genome )) + 
  geom_point(size=3, pch=21) +
  scale_fill_manual("genome used as reference",values=c("black","grey90")) + 
  labs(title="Mapping bias of 20 replicates using one genome as reference only" , x = "total read counts mapped to Yps128 allele (million counts)" , y = "total read counts mapped to Rm11-1a allele (million counts)" ) + 
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
    legend.text = element_text(size=14,face="plain"),
    plot.title = element_text(size = 18, face = "bold",hjust=.5),
    axis.text.x = element_text(colour="grey20",size=13,hjust=.5,vjust=.5,face="plain"),
    axis.text.y = element_text(colour="grey20",size=13,hjust=.5,vjust=.5,face="plain"),
    axis.title.x = element_text(colour="grey20",size=15,hjust=.5,vjust=0,face="bold"), 
    axis.title.y = element_text(colour="grey20",size=15,hjust=.5,vjust=.5,face="bold"))


######### filter snps with no evidence ######## 
evd_rcd_yA <- c()  #check yps_altnative hits match with rm 
evd_rcd_rA <- c()
for(i in 1:nrow(cal)) {
  if(cal[i,"ypsN.x"] == '.') {
    alt_ans <- "IST"
  } else if(cal[i,"rmN.x"] == '.'){
    alt_ans <- "DEL" 
  } else {
    alt_ans <- cal[i,"rmN.x"]
  }  
  
  ck1 <- sum(startsWith(as.character(cal[i,grep("^y.*af",names(cal),value=T)]), alt_ans),na.rm=T) 
  
  if(cal[i,"rmN.y"] == '.') {
    alt_ans <- "IST"
  } else if(cal[i,"ypsN.y"] == '.') {
    alt_ans <- "DEL" 
  } else {
    alt_ans <- cal[i,"ypsN.y"]
  }
  
  ck2 <- sum(startsWith(as.character(cal[i,grep("^r.*af",names(cal),value=T)]),alt_ans),na.rm = T)
  
  evd_rcd_yA <- c(evd_rcd_yA,ck1)
  evd_rcd_rA <- c(evd_rcd_rA,ck2)
} 

no_evd_rows <- which(evd_rcd_yA < 29 | evd_rcd_rA < 29 ) 

cal_evd <- cal[-no_evd_rows,]

save(cal_evd, file = "~/cloud/project/F1_snpfilter/3_5Bsnp/cal_evd.RData")
 #### split blocks yps 
mydisc <- function(chr1,pos1,chr2,pos2) {
  if(chr1 != chr2) {
    return(Inf)
  }else{
    return(abs(as.numeric(pos1)-as.numeric(pos2)))
  }
}

tf2groupName <- function(x) {
  ans <- vector(length=length(x))
  k <- 1 
  ans[1] <- 1 
  for(i in 2:length(x)) {
    if(x[i] == TRUE) {
      k <- k+1
      ans[i] <- k
    }else{
      ans[i] <- ans[i-1]
    }
  }
  return(ans)
}

#### yps_block ###  
yps_block <- cal_evd[,1:4] %>% arrange(ypsChrom,ypsPosit)
yps_block_tf <- vector(length=nrow(yps_block))
yps_block_tf[1] <- TRUE  
for(i in 2:nrow(yps_block)) {
  yps_block_tf[i] <- 
    mydisc(yps_block[i-1,"ypsChrom"], yps_block[i-1,"ypsPosit"],
             yps_block[i,"ypsChrom"], yps_block[i,"ypsPosit"]) > 700000
} 

yps_rm_69293_group <- cbind(yps_block,yps_block_tf,gN= tf2groupName(yps_block_tf))
options(scipen=999)

write.table(unique(yps_rm_69293_group[,c(1,2)]),file="~/cloud/project/snpfilter/3_5Bsnp/yps128_5_snpls_evd",row.names = F,quote=F,col.names = F)

write.table(unique(yps_rm_69293_group[,c(3,4)]),file="~/cloud/project/snpfilter/3_5Bsnp//rm11_B_snpls_evd",row.names = F,quote=F,col.names = F)

write.table(unique(yps_rm_69293_group[,c(1,2,6)]),file="~/cloud/project/snpfilter/3_5Bsnp/yps128_5_snpls_evd_group",row.names = F,quote=F,col.names = F)

write.table(unique(yps_rm_69293_group[,c(3,4,6)]),file="~/cloud/project/snpfilter/3_5Bsnp/rm11_B_snpls_evd_group",row.names = F,quote=F,col.names = F)



######### geneName add ########## 
add_gene_name_snpls <- function() {
genN <- read.table("~/cloud/project/F1_snpfilter/yps128_5_snp_gene.txt",header = F)
overlap_gene <- c()
for(i in 2:nrow(genN)) {
  if(genN[i,2] == genN[i-1,2] ){
    print(unname(genN[i-1,]))
    print(unname(genN[i,]))
    overlap_gene <- c(overlap_gene,i-1,i)
  }
}
geneNu <- genN[-overlap_gene,]


yps_rm_487_gene <- merge.data.frame(yps_block,geneNu,by.x = c(1,2),by.y= c(1,2),sort=F,all=F)
names(yps_rm_487_gene)[5] <- "geneN"


write.table(yps_rm_487_gene,file="~/cloud/project/F1_snpfilter/3_5Bsnp/yr_gene_487",row.names = F,quote=F,col.names = F)

write.table(unique(yps_rm_487_gene[,c(1,2)]),file="~/cloud/project/F1_snpfilter/3_5Bsnp/yps128_5_snpls_u487",row.names = F,quote=F,col.names = F)

write.table(unique(yps_rm_487_gene[,c(3,4)]),file="~/cloud/project/F1_snpfilter/3_5Bsnp/rm11_B_snpls_u487",row.names = F,quote=F,col.names = F)

write.table(unique(yps_rm_487_gene[,c(1,2,5)]),file="~/cloud/project/F1_snpfilter/3_5Bsnp/yps128_5_snpls_u487_gene",row.names = F,quote=F,col.names = F)

write.table(unique(yps_rm_487_gene[,c(3,4,5)]),file="~/cloud/project/F1_snpfilter/3_5Bsnp/rm11_B_snpls_u487_gene",row.names = F,quote=F,col.names = F)
}
