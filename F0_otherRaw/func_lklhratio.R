library(bbmle)
require("foreach")
require("iterators")

dbetabinom <- function(k, n, a , b , log = TRUE) {
    if (log == TRUE) {
        return(lchoose(n = n, k = k) + lbeta(k + a, n-k+ b) - lbeta(a, b))
    } else {
        return(choose(n = n, k = k) * (beta(k+ a, n-k+ b) /  beta(a, b)))
    }
}


neglhbetaBinomial_cisonly <<- function(log_ec,log_rHy,xHy,nHy,C1,C2) {
    ec  <- 2^log_ec
    rHy <- 2^log_rHy
    a <- C1/rHy * ec /(ec + 1 )
    b <- C2/rHy /(ec + 1)
    result <- -sum(dbetabinom(k=xHy, n=nHy, a=a, b=b,log = TRUE))
    return(result)
}

paraGetBB_lkr <- function(oneGeneRow ,x_key,y_key,csource) { 
    xHy=unlist(round(oneGeneRow[x_key]))
    nHy=unlist(round(oneGeneRow[y_key] + oneGeneRow[x_key]))
    
    C1 <- csource[x_key]
    C2 <- csource[y_key]
    
    expCis <- list(xHy=xHy,
                   nHy=nHy,
                   C1=C1,
                   C2=C2)
    
    #print(expCis$C1)
    #print(expCis$C2)
    #print(expCis$xHy)
    #print(expCis$nHy)
    
    # split line """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    fit.beta.cis.nb <- mle2(
        minuslogl = neglhbetaBinomial_cisonly,
        method="SANN",
        #method = "nlminb",
        #optimizer = "nlminb",
        #lower     = list(log_ec = -20,log_rHy = -20),
        start     = list(log_ec = 0, log_rHy =15),
        #upper     = list(log_ec =  20,log_rHy =  30),
        data      = expCis
    )
    
    
   fit.beta.cis.nb.ec0 <- mle2(
       minuslogl = neglhbetaBinomial_cisonly,
       method="SANN",
       #method = "nlminb",
       #optimizer = "nlminb",
       fix = list(log_ec = 0 ), 
       #lower     = list(log_rHy = -20),
       start     = list(log_rHy =15),
       #upper     = list(log_rHy =  30),
       data      = expCis
   )
   
   #fit.beta.cis.nb.r0 <- mle2(
    #   minuslogl = neglhbetaBinomial_cisonly,
    #   method="SANN",
       #method = "nlminb",
       #optimizer = "nlminb",
    #   fix = list(log_rHy = -20 ), 
       #lower     = list(log_rHy = -20),
     #  start     = list(log_ec =0),
       #upper     = list(log_rHy =  30),
      # data      = expCis
   #)
   tmp1 <- coef(fit.beta.cis.nb)
   tmp2 <- coef(fit.beta.cis.nb.ec0)
   #tmp3 <- coef(fit.beta.cis.nb.r0)
   #lnR <-  logLik(fit.beta.cis.nb) - logLik(fit.beta.cis.nb.base)
   #pval <- pchisq(2 * lnR, df=1, lower.tail = FALSE)
   
   c(tmp1, logLik(fit.beta.cis.nb), tmp2, logLik(fit.beta.cis.nb.ec0))
   
}


###  Also test binomial model ##########  
neglhBinomial_cisonly <- function(log_ec,xHy,nHy,C1,C2){
    ec <- 2^log_ec
    pHy <- C1 * ec / (C1 * ec + C2 )
    result <- -sum(dbinom(x = xHy, size = nHy, prob = pHy, log = TRUE))
    return(result)
}

paraGetB_lkr <- function(oneGeneRow,x_key,y_key,csource) {
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
        method = "SANN",
        #lower     = list(log_ec = -20),
        start     = list(log_ec =  0),
        #upper     = list(log_ec =  20),
        data      = expCis)
    
    
    baselike <- (-1) * neglhBinomial_cisonly(log_ec = 0, xHy=expCis$xHy, nHy=expCis$nHy, C1 = expCis$C1, C2= expCis$C2)
    
    
    
    tmp1 <- coef(fit.binom.cis)
    #lnR <- logLik(fit.binom.cis) - baselike
    #pval <- pchisq(2 * lnR, df=1, lower.tail = FALSE)
    
    c(tmp1, logLik(fit.binom.cis), baselike)
}


test_count_row <- c(71,67,67,69, c(71,67,67,69)*1.5)
x_key <- c(1,2,3)
y_key <- c(4,5,6)
csource <- rep(20000,6)
paraGetB_lkr(test_count_row,x_key,y_key,csource=csource)
paraGetBB_lkr(test_count_row,x_key,y_key,csource=csource)

lamd12 <- sum(test_count_row)
allexp <- sum(csource)/2
lambda12 <- lamd12/allexp
theta <-  6.6
2 ^ theta * lambda12 / (1 + 2 ^ theta * lambda12)
