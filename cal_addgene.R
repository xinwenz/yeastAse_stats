snpGeneName <- read.table(header=F,file="/cloud/project/otherRaw/yps128_5_snp_geneName.bed",stringsAsFactors = F)
colnames(snpGeneName) <- c("ypsChrom","ypsPosit","ypsTript")
calg <- merge.data.frame(snpGeneName,cal[-wrong_rows,],by = c(1,2),sort = F,all=F) 
calh <- cbind(ypsGene=calg[,3],calg[,grep("[HA]_rc$",names(cal),value=T)])

x <- colSums(calh[,grep("^y.*[HA]_rc",names(calh),value=T)],na.rm=T)
y <- colSums(calh[,grep("^r.*[HA]_rc",names(calh),value=T)],na.rm=T)
xy <- data.frame(x,y)
ggplot(xy,aes(x=x,y=y)) + geom_point() + labs(title=" Total DNA Read counts in 20 hybrid samples" , x = "total read counts mapped to Yps128 allele" , y = "total read counts mapped to Rm11-1a allele" ) + scale_y_continuous(labels = scales::scientific) + scale_x_continuous(labels = scales::scientific) + geom_abline(slope = 1,intercept = 0) + theme(text=element_text(size=15))

library(MASS)
midmean <- function(x) {
  x[is.na(x)] <- 0 
  #x <- x[x >=  2 ]
  if(length(x) == 0) {return(0)}
  round(fitdistr(x,"poisson")$estimate)
  #max(x)
  #res <- sum(x,rm.na=T)
  #if(is.na(res)) {res <- 1}
  #return(res)
}

cali <- calh %>% group_by(ypsGene) %>% summarise_all(funs(midmean))
tmp <- calh %>% group_by(ypsGene) %>% summarise(n=n())
hist(tmp$n,breaks=50) ## number of snps, a gene have. 
#something is not quite right with aggregate
#exp_mid <- aggregate(. ~ ypsTript,exp_cvmatch_allmatch_tranID,max)


x <- colSums(cali[,grep("^y.*[HA]_rc",names(cali),value=T)],na.rm=T)
y <- colSums(cali[,grep("^r.*[HA]_rc",names(cali),value=T)],na.rm=T)
xy <- data.frame(x,y)
ggplot(xy,aes(x=x,y=y)) + geom_point() + labs(title=" Total Expresssion Read counts in 20 hybrid samples" , x = "total read counts mapped to Yps128 allele" , y = "total read counts mapped to Rm11-1a allele" ) + scale_y_continuous(labels = scales::scientific) + scale_x_continuous(labels = scales::scientific) + geom_abline(slope = 1,intercept = 0) + theme(text=element_text(size=15))

cali_d <- x/y
