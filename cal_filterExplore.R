rm(list = grep(pattern="cal.*",ls(),value=T,invert = T))
library(ggplot2)
yaN <- grep("^y.*af",names(cal),value=T)
rrN <- grep("^r.*rf",names(cal),value=T)
raN <- grep("^r.*af",names(cal),value=T)
yrN <- grep("^y.*rf",names(cal),value=T)

## check if the evidence from reads match with the alignment result. 
evd_rcd_yA <- c()
evd_rcd_rA <- c()
for(i in 1:nrow(cal)) {
  if(cal[i,"ypsN.x"] == '.') {
    alt_ans <- "IST"
  } else if(cal[i,"rmN.x"] == '.'){
      alt_ans <- "DEL" 
  } else {
    alt_ans <- cal[i,"rmN.x"]
  }  
  
  ck1 <- sum(startsWith(as.character(cal[i,yaN]), alt_ans),na.rm=T) 
  
  if(cal[i,"rmN.y"] == '.') {
    alt_ans <- "IST"
  } else if(cal[i,"ypsN.y"] == '.') {
    alt_ans <- "DEL" 
  } else {
    alt_ans <- cal[i,"ypsN.y"]
  }
  
  ck2 <- sum(startsWith(as.character(cal[i,raN]),alt_ans),na.rm = T)
  
  evd_rcd_yA <- c(evd_rcd_yA,ck1)
  evd_rcd_rA <- c(evd_rcd_rA,ck2)
} 

wrong_rows <- which(evd_rcd_yA < 29 | evd_rcd_rA < 29 ) 
# ref count smaller than alternative count.  Probablliy meands the ref info in wrong in the assembly#####  

### Prepare to check allle matching ####
x <- colSums(cal[,grep("^y.*[HA]_rc",names(cal),value=T)],na.rm=T)
y <- colSums(cal[,grep("^y.*[HA]_ac",names(cal),value=T)],na.rm=T)
xy <- data.frame(x,y)
ggplot(xy,aes(x=x,y=y)) + geom_point() + labs(title=" Total Read counts in 20 hybrid samples using YPS as reference genome" , x = "total read counts mapped to Yps128 allele" , y = "total read counts mapped to Rm11-1a allele" ) + scale_y_continuous(labels = scales::scientific) + scale_x_continuous(labels = scales::scientific) + geom_abline(slope = 1,intercept = 0) + theme(text=element_text(size=15))

x <- colSums(cal[,grep("^r.*[HA]_ac",names(cal),value=T)],na.rm=T)
y <- colSums(cal[,grep("^r.*[HA]_rc",names(cal),value=T)],na.rm=T)
xy <- data.frame(x,y)
ggplot(xy,aes(x=x,y=y)) + geom_point() + labs(title=" Total Read counts in 20 hybrid samples using RM11  as reference" , x = "total read counts mapped to Yps128 allele" , y = "total read counts mapped to Rm11-1a allele" ) + scale_y_continuous(labels = scales::scientific) + scale_x_continuous(labels = scales::scientific) + geom_abline(slope = 1,intercept = 0) + theme(text=element_text(size=15))

x <- colSums(cal[,grep("^y.*[HA]_rc",names(cal),value=T)],na.rm=T)
y <- colSums(cal[,grep("^r.*[HA]_rc",names(cal),value=T)],na.rm=T)
xy <- data.frame(x,y)
ggplot(xy,aes(x=x,y=y)) + geom_point() + labs(title=" Total Read counts taking only ref count " , x = "total read counts mapped to Yps128 allele" , y = "total read counts mapped to Rm11-1a allele" ) + scale_y_continuous(labels = scales::scientific) + scale_x_continuous(labels = scales::scientific) + geom_abline(slope = 1,intercept = 0) + theme(text=element_text(size=15))




x <- colSums(cal[-wrong_rows,grep("^y.*[HA]_rc",names(cal),value=T)],na.rm=T)
y <- colSums(cal[-wrong_rows,grep("^r.*[HA]_rc",names(cal),value=T)],na.rm=T)
xy <- data.frame(x,y)
ggplot(xy,aes(x=x,y=y)) + geom_point() + labs(title=" Total Read counts in 20 hybrid samples for correct SNPs and taking only ref count" , x = "total read counts mapped to Yps128 allele" , y = "total read counts mapped to Rm11-1a allele" ) + scale_y_continuous(labels = scales::scientific) + scale_x_continuous(labels = scales::scientific) + geom_abline(slope = 1,intercept = 0) + theme(text=element_text(size=15))

cali_d <- x/y


# filter some posmap 
posmap_52filter <- cal[-wrong_rows,1:9]
save(posmap_52filter,file = "/cloud/project/otherRaw/posmap_52filter.RData")
