x <- colSums(expf[,grep("^y.*[HA]_rc",names(expf),value=T)],na.rm=T)
y <- colSums(expf[,grep("^y.*[HA]_ac",names(expf),value=T)],na.rm=T)
xy <- data.frame(x,y)
ggplot(xy,aes(x=x,y=y)) + geom_point() + labs(title=" total Expression Read counts in 20 hybrid samples using YPS as reference genome" , x = "total read counts mapped to Yps128 allele" , y = "total read counts mapped to Rm11-1a allele" ) + scale_y_continuous(labels = scales::scientific) + scale_x_continuous(labels = scales::scientific) + geom_abline(slope = 1,intercept = 0) + theme(text=element_text(size=15))



x <- colSums(expf[,grep("^r.*[HA]_ac",names(expf),value=T)],na.rm=T)
y <- colSums(expf[,grep("^r.*[HA]_rc",names(expf),value=T)],na.rm=T)
xy <- data.frame(x,y)
ggplot(xy,aes(x=x,y=y)) + geom_point() + labs(title=" Total Expresssion Read counts in 20 hybrid samples using RM11  as reference" , x = "total read counts mapped to Yps128 allele" , y = "total read counts mapped to Rm11-1a allele" ) + scale_y_continuous(labels = scales::scientific) + scale_x_continuous(labels = scales::scientific) + geom_abline(slope = 1,intercept = 0) + theme(text=element_text(size=15))

x <- colSums(expf[,grep("^y.*[HAC]_rc",names(expf),value=T)],na.rm=T)
y <- colSums(expf[,grep("^r.*[HAC]_rc",names(expf),value=T)],na.rm=T)
xy <- data.frame(x,y)
ggplot(xy,aes(x=x,y=y)) + geom_point() + labs(title=" Total Expression Read counts in 20 hybrid samples for correct SNPs and taking only ref count" , x = "total read counts mapped to Yps128 allele" , y = "total read counts mapped to Rm11-1a allele" ) + scale_y_continuous(labels = scales::scientific) + scale_x_continuous(labels = scales::scientific) + geom_abline(slope = 1,intercept = 0) + theme(text=element_text(size=15))


