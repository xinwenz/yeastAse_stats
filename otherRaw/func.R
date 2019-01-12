library(bbmle)
require("foreach")
require("doParallel")
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
  a <- C1 * ec / (rHy * (ec + 1 ))
  b <- C2 * 1 / (rHy * (ec + 1))
  result <- -sum(dbetabinom(k=xHy, n=nHy, a=a, b=b,log = TRUE))
  return(result)
}

paraGetBB <- function(oneGeneRow ,x_key,y_key,csource) { 
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
    #method = "L-BFGS-B",
    optimizer = "nlminb",
    lower     = list(log_ec = -20,log_rHy = -20),
    start     = list(log_ec =   0,log_rHy =   0),
    upper     = list(log_ec =  20,log_rHy =  20),
    data      = expCis)
  
  
  tmp1 <- tryCatch({
    coef(fit.beta.cis.nb)
  },error=function(e){
    rbind(log_ec=NA,log_rHy=NA)
  })
  
  tmp2 <- tryCatch({
    confint(fit.beta.cis.nb, level = 0.95)
  },error=function(e){
    rbind(log_ec=c(NA,NA),log_rHy=c(NA,NA))
  })
  
  if(typeof(tmp2)=="S4") {tmp2 <-rbind(log_ec=c(NA,NA),log_rHy=c(NA,NA))}
  # split line """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  if(is.na(tmp2["log_ec",1])) {
    fit.beta.cis.LB <- mle2(
      minuslogl = neglhbetaBinomial_cisonly,
      method = "L-BFGS-B",
      #optimizer = "nlminb",
      lower     = list(log_ec = -20,log_rHy = -20),
      start     = list(log_ec =   0,log_rHy =   0),
      upper     = list(log_ec =  20,log_rHy =  20),
      data      = expCis)
    
    tmp1 <- tryCatch({
      coef(fit.beta.cis.LB)
    },error=function(e){
      rbind(log_ec=NA,log_rHy=NA)
    })
    
    tmp2 <- tryCatch({
      confint(fit.beta.cis.LB,level=0.95)
    },error=function(e){
      rbind(log_ec=c(NA,NA),log_rHy=c(NA,NA))
    })
    
    if(typeof(tmp2)=="S4") {tmp2 <-rbind(log_ec=c(NA,NA),log_rHy=c(NA,NA))}
  }
  
  # """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  if(is.na(tmp2["log_ec",1])) {
    fit.beta.cis.cg <- mle2(
      minuslogl = neglhbetaBinomial_cisonly,
      method = "CG",
      #optimizer = "nlminb",
      #lower     = list(log.ec = -20,log.rHy = -20),
      start     = list(log_ec = 0,log_rHy = 0),
      #upper     = list(log.ec =  20,log.rHy =  20),
      data      = expCis)
    
    tmp1 <- tryCatch({
      coef(fit.beta.cis.cg)
    },error=function(e){
      rbind(log_ec=NA,log_rHy=NA)
    })
    
    tmp2 <- tryCatch({
      confint(fit.beta.cis.cg,level = 0.95)
    },error=function(e){
      rbind(log_ec=c(NA,NA),log_rHy=c(NA,NA))
    })
    
    if(typeof(tmp2)=="S4") {tmp2 <-rbind(log_ec=c(NA,NA),log_rHy=c(NA,NA))}
  }
  #""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  c(tmp1["log_ec"],tmp2["log_ec",1],tmp2["log_ec",2],
    tmp1["log_rHy"],tmp2["log_rHy",1],tmp2["log_rHy",2])
}


###  Also test binomial model ##########  
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