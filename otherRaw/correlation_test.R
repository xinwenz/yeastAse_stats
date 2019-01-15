library(reshape2)
library(ggplot2)

load('/cloud/project/otherRaw/rep44.RData')

melted_cormat <- melt(cor(expi[,-1],method='pearson'))
ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + labs(title=" pearson correlation of expression counts" )


# cali data has low correlation between each replicates.. May be the coverage is so low, that all numbers centered around 25...  
cal_melt_pear <- melt(cor(cali[,-1],method="pearson"))
ggplot(cal_melt_pear,aes(x=Var1,y=Var2,fill=value)) + geom_tile() + labs(title=" pearson correlation of calibraration DNA counts" )


rep44_melt <- melt(cor(rep44[,-1],method='pearson'))
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


rep44_match <- find_best_match(cor_rep44)
