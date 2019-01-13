library(reshape2)
library(ggplot2)

a <- round(cor(expi[,-1]),3)
write.csv(a,file='a.csv')

melted_cormat <- melt(cor(expi[,-1]))
ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()

b <- cor(cali[,-1])
cal_melt <- melt(cor(cali[,-1]))
ggplot(cal_melt,aes(x=Var1,y=Var2,fill=value)) + geom_tile()

rep48_melt <- melt(cor(repfine[,-1]))
ggplot(data = rep48_melt, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()
