setwd('~/cloud/project/otherRaw')
# remove mummer unmatch snps   
yr_pos <- read.table(header=T,file='yr5B.posmap')
names(yr_pos)[1:6] <- c('ypsChrom','ypsPosit','ypsN','rmChrom','rmPosit','rmN')
ry_pos <- read.table(header=T,file='ryB5.posmap')
names(ry_pos)[1:6] <- c('rmChrom','rmPosit','rmN','ypsChrom','ypsPosit','ypsN')
ry_pos_exchg <- ry_pos[,c(4,5,6,1,2,3,7)]
res <- merge.data.frame(yr_pos,ry_pos_exchg,by.x = c(1,2,4,5),by.y=c(1,2,4,5),all=T,sort=F)
res_good <- res[which(res$drct.x == res$drct.y),] # 86351




yps_block <- res_good[,1:4] %>% arrange(ypsChrom,ypsPosit)
yps_block_tf <- vector(length=nrow(yps_block))
yps_block_tf[1] <- TRUE  
for(i in 2:nrow(yps_block)) {
    yps_block_tf[i] <- 
        mydisc(yps_block[i-1,"ypsChrom"], yps_block[i-1,"ypsPosit"],
               yps_block[i,"ypsChrom"], yps_block[i,"ypsPosit"]) > 700000
} 

yps_rm_86351_group <- cbind(yps_block,yps_block_tf,gN= tf2groupName(yps_block_tf))

write.table(unique(yps_rm_86351_group[,c(1,2)]),file="~/cloud/project/snpfilter/1_5Bsnp/yps128_5_snpls_nofilter",row.names = F,quote=F,col.names = F)

write.table(unique(yps_rm_86351_group[,c(3,4)]),file="~/cloud/project/snpfilter/1_5Bsnp/rm11_B_snpls_nofilter",row.names = F,quote=F,col.names = F)

write.table(unique(yps_rm_86351_group[,c(1,2,6)]),file="~/cloud/project/snpfilter/1_5Bsnp/yps128_5_snpls_nofilter_group",row.names = F,quote=F,col.names = F)

write.table(unique(yps_rm_86351_group[,c(3,4,6)]),file="~/cloud/project/snpfilter/1_5Bsnp/rm11_B_snpls_nofilter_group",row.names = F,quote=F,col.names = F)
