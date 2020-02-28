setwd('~/cloud/project/otherRaw')
# filter out qulified snps 
# generate snp list for yps and rm 
# condense insertion and deletion into 1 pos. Make rm and yps fair #  
yr_pos <- read.table(header=T,file='yr5B.posmap')
names(yr_pos)[1:6] <- c('ypsChrom','ypsPosit','ypsN','rmChrom','rmPosit','rmN')

ry_pos <- read.table(header=T,file='ryB5.posmap')
names(ry_pos)[1:6] <- c('rmChrom','rmPosit','rmN','ypsChrom','ypsPosit','ypsN')
ry_pos_exchg <- ry_pos[,c(4,5,6,1,2,3,7)]

res <- merge.data.frame(yr_pos,ry_pos_exchg,by.x = c(1,2,4,5),by.y=c(1,2,4,5),all=T,sort=F)
res_good <- res[which(res$drct.x == res$drct.y),] # 86351
res_bad <- res[is.na(res$drct.x != res$drct.y),] # 17620

# combine yps same row. 
combn <- function(df,strn) {
  nr <- nrow(df)
  mpool <- df[1,]
  ans  <- data.frame()
  if(strn=='yps'){pos <- c(1,2)}else{pos <- c(3,4)}
  for (i in 2:nr) {
     #print(pos)
     #print(df[i,pos])
     #print(mpool)
  if(all(df[i,pos] == mpool[1,pos])){
     #print("true the same")
     mpool <- rbind(mpool,df[i,]) 
   }else{
      # deal with mpool
      ans <- rbind(ans,mpool[ceiling(nrow(mpool)/2),])
      # re put mpool
      mpool <- df[i,]
    }
   }
  ans <- rbind(ans,mpool[ceiling(nrow(mpool)/2),])
  return(ans)
}  


library(dplyr)
res4 <- res_good %>% arrange(ypsChrom,ypsPosit)
res5 <- combn(res4,"yps")
res6 <- res5 %>% arrange(rmChrom,rmPosit)
res7 <- combn(res6,"rm")

rm(list=grep("res[456]",ls(),value=T))
w1 <- res7[,c(1,2)]
write.table(w1,file="yps_snpls",row.names = F,quote=F,col.names = F)
write.table(res7[,c(3,4)],file="rm_snpls",row.names = F,quote=F,col.names = F)

posmap52 <- res7[,c(1,2,3,4,5,6,8,9,10)]
names(posmap772752)[9] <- "drct"
save(posmap52,file = "/cloud/project/otherRaw/posmap52.RData")
