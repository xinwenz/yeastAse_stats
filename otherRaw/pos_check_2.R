library(dplyr)
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

#### measure distance #### 
mydisc <- function(chr1,pos1,chr2,pos2) {
    if(chr1 != chr2) {
        return(Inf)
    }else{
        return(abs(as.numeric(pos1)-as.numeric(pos2)))
    }
}



res_good_y <- res_good %>% arrange(ypsChrom,ypsPosit)
yps_dis <- vector(length=nrow(res_good_y))
yps_dis[1] <- mydisc(res_good_y[1,"ypsChrom"], res_good_y[1,"ypsPosit"],
       res_good_y[2,"ypsChrom"], res_good_y[2,"ypsPosit"]) > 75  
for(i in 2: (nrow(res_good_y)-1)) {
    tmp <- 
       ( mydisc(res_good_y[i-1,"ypsChrom"], res_good_y[i-1,"ypsPosit"],
                res_good_y[i,"ypsChrom"], res_good_y[i,"ypsPosit"]) > 75 ) &
       ( mydisc(res_good_y[i+1,"ypsChrom"], res_good_y[i+1,"ypsPosit"],
                res_good_y[i,"ypsChrom"], res_good_y[i,"ypsPosit"]) > 75 )
    yps_dis[i] <- tmp
} 
yps_dis[nrow(res_good_y)] <- mydisc(res_good_y[i,"ypsChrom"], res_good_y[i,"ypsPosit"],
                                    res_good_y[i+1,"ypsChrom"], res_good_y[i+1,"ypsPosit"]) > 75  
res_good_y <- cbind(res_good_y,yps_dis)

res_good_yr <- res_good_y %>% arrange(rmChrom,rmPosit)
rm_dis <- vector(length=nrow(res_good_yr))
rm_dis[1] <- mydisc(res_good_yr[1,"rmChrom"], res_good_yr[1,"rmPosit"],
            res_good_yr[2,"rmChrom"], res_good_yr[2,"rmPosit"]) > 75
for(j in 2: (nrow(res_good_yr)-1)) {
    tmp2 <- 
        ( mydisc(res_good_yr[j-1,"rmChrom"], res_good_yr[j-1,"rmPosit"],
                 res_good_yr[j,"rmChrom"], res_good_yr[j,"rmPosit"]) > 75 ) &
        ( mydisc(res_good_yr[j+1,"rmChrom"], res_good_yr[j+1,"rmPosit"],
                 res_good_yr[j,"rmChrom"], res_good_yr[j,"rmPosit"]) > 75 )
    rm_dis[j] <- tmp2
} 
rm_dis[nrow(res_good_yr)] <-  mydisc(res_good_yr[j+1,"rmChrom"], res_good_yr[j+1,"rmPosit"],
                                     res_good_yr[j,"rmChrom"], res_good_yr[j,"rmPosit"]) > 75
res_good_yr <- cbind(res_good_yr,rm_dis)
