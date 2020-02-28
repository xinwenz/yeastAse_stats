library(bbmle)
neglhBinomial_cisonly <- function(log_ec,xHy,nHy,C1,C2){
    ec <- 2^log_ec
    pHy <- C1 * ec / (C1 * ec + C2 )
    result <- -sum(dbinom(x = xHy, size = nHy, prob = pHy, log = TRUE))
    return(result)
}

paraGetB <- function(oneGeneRow,x_key,y_key,csource) {
    xHy=unlist(round(oneGeneRow[x_key]))
    nHy=unlist(round(oneGeneRow[y_key] + oneGeneRow[x_key]))
    
    C1 <- csource[x_key]
    C2 <- csource[y_key]
    
    expCis <- list(xHy=xHy,
                   nHy=nHy,
                   C1 = C1,
                   C2 = C2)
    
    fit.binom.cis <-  mle2(
        minuslogl = neglhBinomial_cisonly,
        optimizer = "nlminb",
        lower     = list(log_ec = -20),
        start     = list(log_ec =  0),
        upper     = list(log_ec =  20),
        data      = expCis)
    
    tmp1 <- coef(fit.binom.cis)
    
    tmp2 <- tryCatch({
        confint(fit.binom.cis,level=0.95)
    },error=function(e){
        c(NA,NA)
    })
    
    c(tmp1["log_ec"],tmp2[1],tmp2[2])
}

test_count_row <- c(1,700,2,700)
x_key <- c(1,2)
y_key <- c(3,4)
csource <- rep(2000,4)
paraGetB(test_count_row,x_key,y_key,csource=csource)
paraGetBB(test_count_row,x_key,y_key,csource=csource)

test_count_row <- c(1,100,600,2,600,100)
x_key <- c(1,2,3)
y_key <- c(4,5,6)
csource <- rep(2000,6)
paraGetB(test_count_row,x_key,y_key,csource=csource)
paraGetBB(test_count_row,x_key,y_key,csource=csource)

myx <-  c(1,150,550,100000)
myy <-  c(3,550,150, 300000)
myn <-  myx + myy
opt <- sum(myx) /  sum(myn)

dbinom(x = myx, size = myn, prob = 1/4, log = T)
sum(dbinom(x = myx, size = myn, prob = 1/4, log = T))
