setwd("~/research/exp")

########################## get yps genome mapping expression #############
setwd("./star_yps128_5/")
fileNum <- sapply(strsplit(x=list.files(path="."),split = "_"),"[",1)
fileLab <- rep(c("A","H","C"),times=10)
smpID <- paste0(fileNum,fileLab)
file <- paste0(fileNum,"_", fileLab,".pilecount")

exp_yps <- read.table(file[1],header=T,stringsAsFactors = F)
colnames(exp_yps)[-(1:2)] <- paste0("y",smpID[1],"_",c("cv","rf","rc","af","ac")) # coverage,ref,refcount, alter fer, alter count


for (i in 2:30) {
  print(file[i])
  print(smpID[i])
  tmp <- read.table(file[i],header=T,stringsAsFactors = F)
  colnames(tmp)[-c(1,2)] <- paste0("y",smpID[i],"_",c("cv","rf","rc","af","ac"))
  exp_yps <- merge.data.frame(exp_yps,tmp,by.x=c(1,2),by.y=c(1,2),all=F)  # if all=T, then has 87381 snps here; if all=F:42647
}


colnames(exp_yps)[c(1,2)] <- c("ypsChrom","ypsPosit")

###################### get rm genome mapping expression #############
setwd("../star_rm11_B/")
fileNum <- sapply(strsplit(x=list.files(path="."),split = "_"),"[",1)
fileLab <- rep(c("A","H","C"),times=10)
smpID <- paste0(fileNum,fileLab)
file <- paste0(fileNum,"_", fileLab,".pilecount")

exp_rm <- read.table(file[1],header=T,stringsAsFactors = F)
colnames(exp_rm)[-(1:2)] <- paste0("r",smpID[1],"_",c("cv","rf","rc","af","ac")) # coverage,ref,refcount, alter fer, alter count

for (i in 2:30) {
  print(file[i])
  print(smpID[i])
  tmp <- read.table(file[i],header=T,stringsAsFactors = F)
  colnames(tmp)[-c(1,2)] <- paste0("r",smpID[i],"_",c("cv","rf","rc","af","ac"))
  exp_rm <- merge.data.frame(exp_rm,tmp,by.x=c(1,2),by.y=c(1,2),all=F)  # if all=T, then has 87381 snps here; if all=F
}

colnames(exp_rm)[c(1,2)] <- c("rmChrom","rmPosit")

#### merge yps_rm data based on posion map ##########
setwd('../')
posmap <- read.table(header=T,file='yr.posmap',stringsAsFactors = F)
names(posmap)[1:4] <- c('ypsChrom','ypsPosit','rmChrom','rmPosit')

mer1 <- merge.data.frame(posmap,exp_yps,by = c(1,2),sort = F,all = F)
mer2 <- merge.data.frame(mer1,exp_rm,by.x = c(3,4),by.y = c(1,2), sort = F, all = F )
exp <- mer2[,c(3,4,1,2,5:305)] # change order of columns 

#####################tr
############# prepare to check allele matching #### 
library('dplyr')
setwd("~/research/exp/")

Mode <- function(x) {
  if(length(x) == 0) {return('None')}
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

get_alt <- function(ref, drt) { 
  if (drt == 1) {return(ref)}
  else if (drt == -1) {
    switch(ref,"A"="T","T"="A","C"="G","G"="C")
  } else {stop("drt not 1 or -1 ")}
}

get_true_alt <- function(ref,drt,myvalue){
  cld1 <- get_alt(ref,drt)
  cld2 <- Mode(myvalue[myvalue!='None'])
  if(cld2 == 'DEL' && cld1 != 'DEL') {return(cld2)}
  else {return(cld1)}
}

isbad <- function(x,real_value) {
  res <- x != real_value &  x != "None"
  return(names(x)[res])
}


######### 2: compare good/bad coverage NUCLETIDE  the nucletide is the same ########### 
# check nucletide to see , if yalt is rref && if ralt is yref ##, if in one experiment is not the same, change to couunt to NA. 

yref_key <- grep("^y.*rf$",names(exp),value=T)
yref_value <- exp[,yref_key]
sum(apply(X=yref_value,1,FUN = function(x) min(x)==max(x)))  # 6788: yps reference mathches 

rref_key <- grep("^r.*rf$",names(exp),value=T)
rref_value <- exp[,rref_key]
sum(apply(X=rref_value,1,FUN=function(x) min(x) == max(x)))  #  6788: rm reference mathches 

yalt_key <- grep("^y.*af$",names(exp),value=T)
tmp <- apply(X=exp,1,FUN = function(x) {true_value <- get_true_alt(x['r10A_rf'],as.integer(x['drct']),x[yalt_key]); y <- x[yalt_key]; y <- y[y != 'None'];all(y == true_value)})
sum(tmp)  # 4376 : alternative not match all,
yaltNotMatch <-  names(tmp)[which(tmp == FALSE)]# there're about -2400 has exp that doesn't match, I will transfer thmem to NA 

ralt_key <- grep("^r.*af$",names(exp),value=T)
tmp2 <- apply(X=exp,1,FUN = function(x) {true_value <- get_true_alt(x['y10A_rf'],as.integer(x['drct']),x[ralt_key]); y <- x[ralt_key]; y <- y[y != 'None']; all( y == true_value)})
#tmp2 <- apply(X=exp_cv_good[2241,],1,FUN = function(x) {print(x['y10A_rf']);print(x['drct']);print(x[ralt_key])})
sum(tmp2) # 1475 , alternative not match all match with true value 
raltNotMatch <- names(tmp2)[which(tmp2 == FALSE)]  # 3000 doesn't match 

exp_aleRev <-exp 

for(i in 1:nrow(exp)) {
  x <- unlist(exp[i,])
  #print(x)
  true_value_yalt <- get_true_alt(x['r10A_rf'],as.integer(x['drct']),x[yalt_key])
  posTochange_yalt <- isbad(x[yalt_key],true_value_yalt)
  countPosTochange_yalt <- gsub('f','c',posTochange_yalt)
  #print(countPosTochange_yalt)
  exp_aleRev[i,countPosTochange_yalt] <- NA
  
  true_value_ralt <- get_true_alt(x['y10A_rf'],as.integer(x['drct']),x[ralt_key])
  posTochange_ralt <- isbad(x[ralt_key],true_value_ralt)
  countPosTochange_ralt <- gsub('f','c',posTochange_ralt)
  #print(countPosTochange_ralt)
  exp_aleRev[i,countPosTochange_ralt] <- NA
  
  if(i %% 100 == 0){print(i)}
}

############### some plots before any filetering ###########
######### 3. draw plot of snps without filtering only hybrid ##############
library("ggplot2")
######## fig S4 ##########
yrefC_key <- grep("^y.*[HA]_rc$",names(exp_aleRev),value=T)
yrefC_value <- exp_aleRev[,yrefC_key]

yaltC_key <- grep("^y.*[HA]_ac$",names(exp_aleRev),value=T)
yaltC_value <- exp_aleRev[,yaltC_key]

x <- colSums(yrefC_value)
y <- colSums(yaltC_value,na.rm = T)
xy <- data.frame(x,y)
ggplot(xy,aes(x=x,y=y)) + geom_point() + labs(title=" Total Read counts in 20 hybrid samples" , x = "total read counts mapped to Yps128 allele" , y = "total read counts mapped to Rm11-1a allele" ) + scale_y_continuous(labels = scales::scientific) + scale_x_continuous(labels = scales::scientific) + geom_abline(slope = 1,intercept = 0) + theme(text=element_text(size=15))
#plot(x,y)
#abline(a=0,b=1)

######### fig S4 ###########
rrefC_key <- grep("^r.*[HA]_rc$",names(exp_aleRev),value=T)
rrefC_value <- exp_aleRev[,rrefC_key]

raltC_key <- grep("^r.*[HA]_ac$",names(exp_aleRev),value=T)
raltC_value <- exp_aleRev[,raltC_key]

y <- colSums(rrefC_value)
x <- colSums(raltC_value,na.rm = T)
xy <- data.frame(x,y)
ggplot(xy,aes(x=x,y=y)) + geom_point() + labs(title=" Total Read counts in 20 hybrid samples" , x = "total read counts mapped to Yps128 allele" , y = "total read counts mapped to Rm11-1a allele" ) + scale_y_continuous(labels = scales::scientific) + scale_x_continuous(labels = scales::scientific) + geom_abline(slope = 1,intercept = 0) + theme(text=element_text(size=15))


####### check allele match, and filter out the ones that doesn't match ####33  
### filter snps of allele matching #####
yrefC_key <- grep(pattern="y.*[HA]_rc",x=names(exp_aleRev),value=T)
yrefC_value <- exp_aleRev[,yrefC_key]

raltC_key <- grep(pattern= "r.*[HA]_ac" , x = names(exp_aleRev),value=T)
raltC_value <- exp_aleRev[,raltC_key]

yaltC_key <- grep(pattern="y.*[HA]_ac",x=names(exp_aleRev),value=T)
yaltC_value <- exp_aleRev[,yaltC_key]

rrefC_key <- grep(pattern="r.*[HA].rc",x=names(exp_aleRev),value=T)
rrefC_value <- exp_aleRev[,rrefC_key]



goodnum <- nrow(yrefC_value)

match_0  <- vector(,goodnum)
match_5  <- vector(,goodnum)
match_10 <- vector(,goodnum)
match_15 <- vector(,goodnum)


for( i in 1:goodnum){
  tmp1 <- yrefC_value[i,] == raltC_value[i,]
  tmp2 <- yaltC_value[i,] == rrefC_value[i,]
  match_0[i] <- sum(tmp1,tmp2,na.rm=T)
  
  tmp1 <- abs( yrefC_value[i,] - raltC_value[i,] ) < (yrefC_value[i,] + raltC_value[i,])/2 * 0.05
  tmp2 <- abs( yaltC_value[i,] - rrefC_value[i,] ) < (yaltC_value[i,] + rrefC_value[i,])/2 * 0.05
  match_5[i] <- sum(tmp1,tmp2,na.rm=T)
  
  tmp1 <- abs( yrefC_value[i,] - raltC_value[i,] ) < (yrefC_value[i,] + raltC_value[i,])/2 * 0.1
  tmp2 <- abs( yaltC_value[i,] - rrefC_value[i,] ) < (yaltC_value[i,] + rrefC_value[i,])/2 * 0.1
  match_10[i] <- sum(tmp1,tmp2,na.rm=T)
  
  tmp1 <- abs( yrefC_value[i,] - raltC_value[i,] ) < (yrefC_value[i,] + raltC_value[i,])/2 * 0.15
  tmp2 <- abs( yaltC_value[i,] - rrefC_value[i,] ) < (yaltC_value[i,] + rrefC_value[i,])/2 * 0.15
  match_15[i] <- sum(tmp1,tmp2,na.rm=T)
  
  
  if(i %% 100 == 0){print(i)}
}

plot(density(match_0),col='orange',ylim=c(0,0.3))
lines(density(match_5),col='green')
lines(density(match_10),col='blue')
lines(density(match_15),col='black')

######  green line to 27..  5% error rate #############  
retain_key <- grep(pattern = ".*[AH].*c$",x=names(exp_aleRev), value = T)

exp_fine <- exp_aleRev[match_10 >= 27,c("ypsChrom","ypsPosit",retain_key)] 

##### check the new filtered data frame ############
yrefC_key <- grep("^y.*[HA]_rc$",names(exp_fine),value=T)
yrefC_value <- exp_fine[,yrefC_key]

yaltC_key <- grep("^y.*[HA]_ac$",names(exp_fine),value=T)
yaltC_value <- exp_fine[,yaltC_key]

x <- colSums(yrefC_value)
y <- colSums(yaltC_value,na.rm = T)
xy <- data.frame(x,y)
ggplot(xy,aes(x=x,y=y)) + geom_point() + labs(title=" Total Read counts in 20 hybrid samples" , x = "total read counts mapped to Yps128 allele" , y = "total read counts mapped to Rm11-1a allele" ) + scale_y_continuous(labels = scales::scientific) + scale_x_continuous(labels = scales::scientific) + geom_abline(slope = 1,intercept = 0) + theme(text=element_text(size=15))

rrefC_key <- grep("^r.*[HA]_rc$",names(exp_fine),value=T)
rrefC_value <- exp_fine[,rrefC_key]

raltC_key <- grep("^r.*[HA]_ac$",names(exp_fine),value=T)
raltC_value <- exp_fine[,raltC_key]

y <- colSums(rrefC_value)
x <- colSums(raltC_value,na.rm = T)
xy <- data.frame(x,y)
ggplot(xy,aes(x=x,y=y)) + geom_point() + labs(title=" Total Read counts in 20 hybrid samples" , x = "total read counts mapped to Yps128 allele" , y = "total read counts mapped to Rm11-1a allele" ) + scale_y_continuous(labels = scales::scientific) + scale_x_continuous(labels = scales::scientific) + geom_abline(slope = 1,intercept = 0) + theme(text=element_text(size=15))

###### combine two reference #########  
ypsAllele <- (yrefC_value + raltC_value) / 2 
names(ypsAllele) <- sapply(strsplit(yrefC_key,split = "_"),"[",1)
rmAllele <- (rrefC_value + yaltC_value) /2 
names(rmAllele)  <- sapply(strsplit(rrefC_key,split = "_"), "[", 1)

x <- colSums(ypsAllele,na.rm= T)
y <- colSums(rmAllele ,na.rm = T)
xy <- data.frame(x,y)
ggplot(xy,aes(x=x,y=y)) + geom_point() + labs(title=" Total Read counts in 20 hybrid samples" , x = "total read counts mapped to Yps128 allele" , y = "total read counts mapped to Rm11-1a allele" ) + scale_y_continuous(labels = scales::scientific) + scale_x_continuous(labels = scales::scientific) + geom_abline(slope = 1,intercept = 0) + theme(text=element_text(size=15))

exp_fine2 <- cbind(exp_fine[,c(1,2)],ypsAllele,rmAllele)

###### add gene name, and find how much gene does it cover ############# 
transId <- read.table(header=F,file="yps128_5_snp_geneName.bed",stringsAsFactors = F)
colnames(transId) <- c("ypsChrom","ypsPosit","ypsTript")
exp_fine3 <- merge.data.frame(transId,exp_fine2,by = c(1,2),sort = F,all=F) # keep all genes and snps, no matter have data or not . 
exp_fine3 <- exp_fine3[,-c(1,2)]

library(MASS)
midmean <- function(x) {
  #x <- x[!is.na(x)]
  #x <- x[x >=  2 ]
  #if(length(x) == 0) {return(0)}
  #fitdistr(x,"poisson")$estimate
  #max(x)
  res <- sum(x,rm.na=T)
  #if(is.na(res)) {res <- 1}
  #return(res)
}

exp_fine4 <- exp_fine3 %>% group_by(ypsTript) %>% summarise_all(funs(midmean))
tmp <- exp_fine3 %>% group_by(ypsTript) %>% summarise(n=n())
hist(tmp$n,breaks=50) ## number of snps, a gene have. 
#something is not quite right with aggregate
#exp_mid <- aggregate(. ~ ypsTript,exp_cvmatch_allmatch_tranID,max)

