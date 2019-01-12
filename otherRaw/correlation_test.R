a <- round(cor(expi[,-1]),3)
write.csv(a,file='a.csv')
library(reshape2)
melted_cormat <- melt(a)
head(melted_cormat)
library(ggplot2)
ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()
