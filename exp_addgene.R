snpGeneName <- read.table(header=F,file="/cloud/project/otherRaw/yps128_5_snp_geneName.bed",stringsAsFactors = F)
colnames(snpGeneName) <- c("ypsChrom","ypsPosit","ypsTript")
expg <- merge.data.frame(snpGeneName,expf,by = c(1,2),sort = F,all=F) 
exph <- cbind(ypsGene=expg[,3],expg[,grep("[HA]_rc$",names(expf),value=T)])

x <- colSums(exph[,grep("^y.*[HA]_rc",names(exph),value=T)],na.rm=T)
y <- colSums(exph[,grep("^r.*[HA]_rc",names(exph),value=T)],na.rm=T)
xy <- data.frame(x,y)
ggplot(xy,aes(x=x,y=y)) + geom_point() + labs(title=" Total Expresssion Read counts in 20 hybrid samples using YPS as reference and correct SNPS" , x = "total read counts mapped to Yps128 allele" , y = "total read counts mapped to Rm11-1a allele" ) + scale_y_continuous(labels = scales::scientific) + scale_x_continuous(labels = scales::scientific) + geom_abline(slope = 1,intercept = 0) + theme(text=element_text(size=15))



library(MASS)
midmean <- function(x) {
  #x <- x[!is.na(x)]
  #x <- x[x >=  2 ]
  if(length(x) == 0) {return(0)}
  round(fitdistr(x,"poisson")$estimate)
  #max(x)
  #res <- sum(x,rm.na=T)
  #if(is.na(res)) {res <- 1}
  #return(res)
}

expi <- exph %>% group_by(ypsGene) %>% summarise_all(funs(midmean))
names(expi)[-1] <- sub("(.*)_rc","\\1",names(expi[,-1]))

       
tmp <- exph %>% group_by(ypsGene) %>% summarise(n=n())
hist(tmp$n,breaks=50) ## number of snps, a gene have. 
#something is not quite right with aggregate
#exp_mid <- aggregate(. ~ ypsTript,exp_cvmatch_allmatch_tranID,max)


x <- colSums(expi[,grep("^y.*[HA]",names(expi),value=T)],na.rm=T)
y <- colSums(expi[,grep("^r.*[HA]",names(expi),value=T)],na.rm=T)
xy <- data.frame(x,y)
ggplot(xy,aes(x=x,y=y)) + geom_point() + labs(title=" Total Expresssion Read counts in 20 hybrid samples" , x = "total read counts mapped to Yps128 allele" , y = "total read counts mapped to Rm11-1a allele" ) + scale_y_continuous(labels = scales::scientific) + scale_x_continuous(labels = scales::scientific) + geom_abline(slope = 1,intercept = 0) + theme(text=element_text(size=15))

