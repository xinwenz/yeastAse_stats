load("~/cloud/project/F1_snpfilter/5_5Bsnp/cal_no_right_tail.RData") 
# get cal_no_right_tail data.frame

cal_a <- cal_no_right_tail[c("ypsChrom","ypsPosit", grep("^y.*[AH]_rc",names(cal_no_right_tail),value=T), grep("^r.*[AH]_rc",names(cal_no_right_tail),value=T))] 

x <- rowSums(cal_a[,c(3:22)[-13]],na.rm = TRUE)
y <- rowSums(cal_a[,c(3:22)[-13]+20],na.rm = TRUE)
cal_pval_df <- data.frame(yps=x,ypsrm = x+y)
cal_pres <- apply(cal_pval_df, 1 , function(x) binom.test(x[1],x[2])$p.value)
cal_p <- cbind(cal_a[,c(1:2)], dna_sum=cal_pval_df,dna_pres=cal_pres)
#pres_rank <- rank(pres)

#cal_filter_10psnps <- cal_with_label[pres_rank > 66352*0.12 ,] 
#hist(which(pres_rank < 66352 * 0.12)) # check the postion of those 

#ax <- rowSums(cal_filter_10psnps[,c(1:20)[-13]],na.rm = TRUE)
#ay <- rowSums(cal_filter_10psnps[,c(1:20)[-13]+20],na.rm = TRUE)
#apval_df <- data.frame(yps=ax,ypsrm = ax+ay)
#apres <- apply(apval_df, 1 , function(x) binom.test(x[1],x[2])$p.value)

###########  load RNA expression SNPS data ################# ###################
########################## get yps genome mapping expression #############
setwd("~/cloud/project/old_files/bowtie2_yps128_5/")
fileNum <- sapply(strsplit(x=list.files(path="."),split = "_"),"[",1)
fileLab <- rep(c("A","H","C"),times=10)
smpID <- paste0(fileNum,fileLab)
file <- paste0(fileNum,"_", fileLab,".pilecount")

exp_yps <- read.table(file[1],header=T,stringsAsFactors = F)
colnames(exp_yps)[-(1:2)] <- paste0("y",smpID[1],"_",c("rf","rc","af","ac")) # coverage,ref,refcount, alter fer, alter count


for (i in 2:30) {
    print(file[i])
    print(smpID[i])
    tmp <- read.table(file[i],header=T,stringsAsFactors = F)
    colnames(tmp)[-c(1,2)] <- paste0("y",smpID[i],"_",c("rf","rc","af","ac"))
    exp_yps <- merge.data.frame(exp_yps,tmp,by.x=c(1,2),by.y=c(1,2),all=F)  # if all=T, then has 87381 snps here; if all=F:42647
}


colnames(exp_yps)[c(1,2)] <- c("ypsChrom","ypsPosit")
###################### get rm genome mapping expression #############
setwd("~/cloud/project/old_files/bowtie2_rm11_B/")
fileNum <- sapply(strsplit(x=list.files(path="."),split = "_"),"[",1)
fileLab <- rep(c("A","H","C"),times=10)
smpID <- paste0(fileNum,fileLab)
file <- paste0(fileNum,"_", fileLab,".pilecount")

exp_rm <- read.table(file[1],header=T,stringsAsFactors = F)
colnames(exp_rm)[-(1:2)] <- paste0("r",smpID[1],"_",c("rf","rc","af","ac")) # coverage,ref,refcount, alter fer, alter count

for (i in 2:30) {
    print(file[i])
    print(smpID[i])
    tmp <- read.table(file[i],header=T,stringsAsFactors = F)
    colnames(tmp)[-c(1,2)] <- paste0("r",smpID[i],"_",c("rf","rc","af","ac"))
    exp_rm <- merge.data.frame(exp_rm,tmp,by.x=c(1,2),by.y=c(1,2),all=F)  # if all=T, then has 87381 snps here; if all=F
}

colnames(exp_rm)[c(1,2)] <- c("rmChrom","rmPosit")

#### merge yps_rm data based on posion map ##########
setwd('~/cloud/project/F0_otherRaw/')
load("posmap52.RData")

mer1 <- merge.data.frame(posmap52,exp_yps,by = c(1,2),sort = F,all = F)
mer2 <- merge.data.frame(mer1,exp_rm,by.x = c(3,4),by.y = c(1,2), sort = F, all = F )
expf <- mer2[,c(3,4,1,2,5:249)] # change order of columns 
#rm(list=ls()[-c(1,2,3)])
j <- sapply(expf,is.factor)
expf[j] <- lapply(expf[j],as.character)


###################  

exp_a <- expf[c("ypsChrom","ypsPosit", grep("^y.*[AH]_rc",names(expf),value=T), grep("^r.*[AH]_rc",names(expf),value=T))] 
# remove -1 
exp_a[,-c(1,2)] <- lapply(exp_a[,-c(1,2)], function(x) pmin(x,0)*(-1) + x )

x <- rowSums(exp_a[,c(3:22)[-c(7,19)]],na.rm = TRUE)
y <- rowSums(exp_a[,c(3:22)[-c(7,19)] +20],na.rm = TRUE)
exp_pval_df <- data.frame(yps=x,ypsrm = x+y)
exp_pres <- apply(exp_pval_df, 1 , function(x) binom.test(x[1],x[2])$p.value)
exp_p <- cbind(exp_a[,c(1:2)],rna_sum=exp_pval_df, rna_pres=exp_pres)

#### merge DNA and RNA p-values #### 
calexp_p <- merge.data.frame(cal_p, exp_p, by = c(1,2), sort = F,all = F)######

tDR <- calexp_p$dna_pres < 0.05
dim(tDR) <- c(nrow(calexp_p),1)
rownames(tDR) <- rownames(calexp_p$ypsPosit)
tDR <- cbind(tDR, calexp_p$rna_pres < 0.05)
colnames(tDR) <- colnames(c( "DNA", "RNA" ))  
cnt <- table(DNA=tDR[,1],RNA=tDR[,2])
fisher.test(cnt)
##### 
