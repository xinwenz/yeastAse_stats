library(reshape2)
library(ggplot2)

a <- round(cor(expi[,-1]),3)
write.csv(a,file='a.csv')

melted_cormat <- melt(cor(expi[,-1],method='spearman'))
ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()


# cali data has low correlation between each replicates.. May be the coverage is so low, that all numbers centered around 25...  
cal_melt_pear <- melt(cor(cali[,-1],method="pearson"))
ggplot(cal_melt_pear,aes(x=Var1,y=Var2,fill=value)) + geom_tile()

cal_melt_kend <- melt(cor(cali[,-1],method="kendall"))
ggplot(cal_melt_kend,aes(x=Var1,y=Var2,fill=value)) + geom_tile()


rep48_melt <- melt(cor(repfine[,-1],method='spearman'))
ggplot(data = rep48_melt, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()
