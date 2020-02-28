n1 <-  which(calh_notail_2650$ecisL > 0 )
n2 <- which(calh_notail_2660$ecisL > 0 )

calh_notail_2650[122,]
calh_notail_2660[122,]
x <- unlist(calh_notail[122,grep(pattern = "*_yGc", names(calh_notail), value=T)])
y <- unlist(calh_notail[122,grep(pattern = "*_rGc", names(calh_notail), value=T)])
xy <- data.frame(x=x,y=y)
ggplot(xy) + geom_point(aes(x,y)) + geom_abline(slope = 1,intercept = 0)



rd1 <- rnorm(n=1, mean =1 , sd =2 )
rd2 <- rnorm(n=100, mean = 1, sd = 2)
rd3 <- rnorm(n=1000, mean =1 , sd = 2)

library('bbmle')

maxlik <- function(md,vx) {
    - sum(  dnorm(vx, mean = md, sd =2, log =T) ) 
}


d1 = data.frame( vx = rd1)
tmp1 <- mle2(minuslogl = maxlik, start= list(md = 0.5), data = d1 )
tmp1_fix <- mle2(minuslogl = maxlik, fixed = list(md = 0), start= list(md = 0.5), data = d1 )

coef(tmp1)
confint(tmp1)
logLik(tmp1) - logLik(tmp1_fix) #1.9


d2 = data.frame( vx = rd2)
tmp2 <- mle2(minuslogl = maxlik, start= list(md = 0.5), data = d2 )
tmp2_fix <- mle2(minuslogl = maxlik, start= list(md = 0.5), fixed = list(md=0), data = d2 )

coef(tmp2)
confint(tmp2)
logLik(tmp2) - logLik(tmp2_fix) #8.5


d3 = data.frame( vx = rd3)
tmp3 <- mle2(minuslogl = maxlik, start= list(md = 0.5), data = d3)
tmp3_fix <- mle2(minuslog = maxlik, start=list(md=0.5), fixed = list(md=0), data = d3) 

coef(tmp3)
confint(tmp3)
logLik(tmp3) - logLik(tmp3_fix) # 117.0235


rd6 <- rnorm(n = 1000000 , mean =1 ,sd = 2)
d6 = data.frame( vx = rd6)
tmp6 <- mle2(minuslogl = maxlik, start= list(md = 0.5), data = d6)
tmp6_fix <- mle2(minuslogl = maxlik, start= list(md=0.5), fixed= list(md=0), data =d6) 

coef(tmp6)
confint(tmp6)
logLik(tmp6) - logLik(tmp6_fix) # 125121.5  

1.03 - 0.95
1.13 - 0.89
1.138 - 0.597 
1.0038 - 0.9959
