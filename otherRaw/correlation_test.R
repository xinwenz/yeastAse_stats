library(reshape2)
library(ggplot2)

load('/cloud/project/otherRaw/rep44.RData')

melted_cormat <- melt(cor(expi[,-1],method='pearson'))
ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + labs(title=" pearson correlation of expression counts" )


# cali data has low correlation between each replicates.. May be the coverage is so low, that all numbers centered around 25...  
cal_melt_pear <- melt(cor(cali[,-1],method="pearson"))
ggplot(cal_melt_pear,aes(x=Var1,y=Var2,fill=value)) + geom_tile() + labs(title=" pearson correlation of calibraration DNA counts" )


merep44_melt <- melt(cor(rep44[,-1],method='pearson'))
ggplot(data = rep44_melt, aes(x=Var1, y=Var2, fill=value)) + geom_tile() +  labs(title=" pearson correlation of 44 replicate expression counts" )
cor_rep44 <- cor(rep44[,-1],method='pearson')


find_best_match <- function(m) {
  n_row <- nrow(m)
  n_col <- ncol(m)
  usedcol <-c ()
  usedrow <- c()
  ans <- c()
  
  i <- 1
  while( i <= n_row -1  ){
    if (i %in% usedrow) {i <- i+1 ; next}
    search_col <- ((i+1):n_col)[!(( (i+1) : n_col ) %in% usedcol)]
    #print(search_col)
    ele <- search_col[which.max(m[i,search_col])]
    usedcol <- c(usedcol,ele)
    usedrow <- c(usedrow,ele)
    #print(ele)
    #print(usedcol)
    #print(c(i,ele,m[i,ele]))
    ans <- rbind(ans, c(i,ele,m[i,ele]))
    i <- i +1 
  }
  return(ans)
}

# three null data set 's match seq 
rep44_match <- find_best_match(cor_rep44) # 0.998
expi_yps_match <- find_best_match(cor(expi[,2:21])) # 0.9913 
expi_rm_match <- find_best_match(cor(expi[,22:41])) # 0.9914

cali_cor <- c()
for(i in 1:20) {
  tmp <- c(i,cor(cali[,i+1],cali[,i+21]))
  cali_cor <- rbind(cali_cor,tmp)
}
cali_cor_match <- cali_cor[order(cali_cor[,2],decreasing = T),] # 0.504


expi_cor <- c()
for(i in 1:20) {
  tmp <- c(i,cor(expi[,i+1],expi[,i+21]))
  expi_cor <- rbind(expi_cor,tmp)
}
expi_cor_match <- expi_cor[order(expi_cor[,2],decreasing = T),] # 0.972

save(list=c("rep44_match","expi_yps_match","expi_rm_match","cali_cor_match","expi_cor_match"),file='/cloud/project/otherRaw/run_hpc_Cor_seq.RData')
