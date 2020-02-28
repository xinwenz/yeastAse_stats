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

setwd('~/cloud/project/otherRaw')
# remove mummer unmatch snps   
yr_pos <- read.table(header=T,file='yr5B.posmap')
names(yr_pos)[1:6] <- c('ypsChrom','ypsPosit','ypsN','rmChrom','rmPosit','rmN')
ry_pos <- read.table(header=T,file='ryB5.posmap')
names(ry_pos)[1:6] <- c('rmChrom','rmPosit','rmN','ypsChrom','ypsPosit','ypsN')
ry_pos_exchg <- ry_pos[,c(4,5,6,1,2,3,7)]
res <- merge.data.frame(yr_pos,ry_pos_exchg,by.x = c(1,2,4,5),by.y=c(1,2,4,5),all=T,sort=F)
res_good <- res[which(res$drct.x == res$drct.y),] # 86351

## rm11_B_bad.pos is from nnp/1_0_0_pion_freebayes vcf file. 
rm11_bad_pos <- read.table("~/cloud/project/snpfilter/2_5Bsnp/rm11_B_bad.pos", header=F, stringsAsFactors = F)
yps128_bad_pos <- read.table("~/cloud/project/snpfilter/2_5Bsnp/yps128_5_bad.pos", header=F,stringsAsFactors = F)
colnames(rm11_bad_pos) <- c("rmChrom","rmPosit")
colnames(yps128_bad_pos) <- c("ypsChrom","ypsPosit")

# rm11_B_bad_position:  within 125 bp , delete them, 125 = 75 + 50 , 
# 50 is the longest indel in vcf file: rm11_B.vcf
# yps128_5_bad_position; within 78 bp, dlelete them, 78 = 75 + 3 
# 3 is  the longest indel in vcf file: yps128_5.vcf

yps128_bad_upper <- data.frame(ypsChrom = yps128_bad_pos$ypsChrom, ypsPosit= yps128_bad_pos$ypsPosit + 78, ypsFlag=1:nrow(yps128_bad_pos) ) 
yps128_bad_lower <- data.frame(ypsChrom = yps128_bad_pos$ypsChrom, ypsPosit = yps128_bad_pos$ypsPosit - 78, ypsFlag=1:nrow(yps128_bad_pos))
yps128_bad_all <- rbind(yps128_bad_upper, yps128_bad_lower)


res1 <- merge.data.frame(res_good,yps128_bad_all,all=T,by=c("ypsChrom","ypsPosit"))
decision <- vector(mode = "integer", length=nrow(res1))
for(i in 1:max(yps128_bad_all$ypsFlag)) {
  lower_bound <- which(res1$ypsFlag == i)[1]
  upper_bound <- which(res1$ypsFlag == i)[2]
  decision[lower_bound:upper_bound] <- i
  print(i)
}
res2 <- res1[-which(decision>0), ]


rm11_bad_upper <- data.frame(rmChrom = rm11_bad_pos$rmChrom, rmPosit = rm11_bad_pos$rmPosit + 125, rmFlag=1:nrow(rm11_bad_pos))
rm11_bad_lower <- data.frame(rmChrom = rm11_bad_pos$rmChrom, rmPosit = rm11_bad_pos$rmPosit - 125, rmFlag= 1 :nrow(rm11_bad_pos))
rm11_bad_all <- rbind(rm11_bad_upper, rm11_bad_lower)

res3 <- merge.data.frame(res2,rm11_bad_all,all=T,by=c("rmChrom","rmPosit"))

decision2 <- vector(mode = "integer", length=nrow(res3))
for(i in 1:max(rm11_bad_all$rmFlag)) {
  lower_bound <- which(res3$rmFlag == i)[1]
  upper_bound <- which(res3$rmFlag == i)[2]
  decision2[lower_bound:upper_bound] <- i
  print(i)
}
res4 <- res3[-which(decision2>0), ]



yps_block <- res4[,c(3,4,1,2)] %>% arrange(ypsChrom,ypsPosit)
yps_block_tf <- vector(length=nrow(yps_block))
yps_block_tf[1] <- TRUE  
for(i in 2:nrow(yps_block)) {
  yps_block_tf[i] <- 
    mydisc(yps_block[i-1,"ypsChrom"], yps_block[i-1,"ypsPosit"],
           yps_block[i,"ypsChrom"], yps_block[i,"ypsPosit"]) > 700000
} 

yps_rm_84349_group <- cbind(yps_block,yps_block_tf,gN= tf2groupName(yps_block_tf))

options(scipen=999)

write.table(unique(yps_rm_84349_group[,c(1,2)]),file="~/cloud/project/snpfilter/2_5Bsnp/yps128_5_snpls_bys",row.names = F,quote=F,col.names = F)

write.table(unique(yps_rm_84349_group[,c(3,4)]),file="~/cloud/project/snpfilter/2_5Bsnp/rm11_B_snpls_bys",row.names = F,quote=F,col.names = F)

write.table(unique(yps_rm_84349_group[,c(1,2,6)]),file="~/cloud/project/snpfilter/2_5Bsnp/yps128_5_snpls_bys_group",row.names = F,quote=F,col.names = F)

write.table(unique(yps_rm_84349_group[,c(3,4,6)]),file="~/cloud/project/snpfilter/2_5Bsnp/rm11_B_snpls_bys_group",row.names = F,quote=F,col.names = F)
