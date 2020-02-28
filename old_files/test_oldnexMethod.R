library(bbmle)

dbetabinom <- function(k, n, a , b , log = TRUE) {
    if (log == TRUE) {
        return(lchoose(n = n, k = k) + lbeta(k + a, n-k+ b) - lbeta(a, b))
    } else {
        return(choose(n = n, k = k) * (beta(k+ a, n-k+ b) /  beta(a, b)))
    }
}


beta_cis_new  <<- function(log_ec,log_rHy,xHy,nHy,C1,C2) {
    ec  <- 2^log_ec
    rHy <- 2^log_rHy
    a <- C1/rHy * ec /(ec + 1 )
    b <- C2/rHy /(ec + 1)
    result <- -sum(dbetabinom(k=xHy, n=nHy, a=a, b=b,log = TRUE))
    return(result)
}


beta_cis_old  <<- function(log_ec,log_rHy,xHy,nHy,C1,C2) {
    ec  <- 2^log_ec
    rHy <- 2^log_rHy
    a <- 1/rHy * ec /(ec + 1 )
    b <- 1/rHy /(ec + 1)
    result <- -sum(dbetabinom(k=xHy, n=nHy, a=a, b=b,log = TRUE))
    return(result)
}


bino_cis_new <- function(log_ec,xHy,nHy,C1,C2){
    ec <- 2^log_ec
    pHy <- C1 * ec / (C1 * ec + C2 )
    result <- -sum(dbinom(x = xHy, size = nHy, prob = pHy, log = TRUE))
    return(result)
}

bino_cis_old <- function(log_ec,xHy,nHy,C1,C2){
    ec <- 2^log_ec
    pHy <- 1 * ec / (1 + ec )
    result <- -sum(dbinom(x = xHy, size = nHy, prob = pHy, log = TRUE))
    return(result)
}

paraGetBB <- function(oneGeneRow,x_key,y_key,csource) { 
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
        minuslogl = beta_cis_new,
        #method="SANN",
        method = "L-BFGS-B",
        #optimizer = "nlminb",
        lower     = list(log_ec = -20,log_rHy = -20),
        start     = list(log_ec =   0,log_rHy =15),
        upper     = list(log_ec =  20,log_rHy =  30),
        data      = expCis
    )
    
    
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
    if( is.na(tmp2["log_ec",1]) )    {
        fit.beta.cis.LB <- mle2(
            minuslogl = beta_cis_new,
            method = "L-BFGS-B",
            #optimizer = "nlminb",
            lower     = list(log_ec = -20,log_rHy = -20),
            start     = list(log_ec =   0,log_rHy = 20),
            upper     = list(log_ec =  20,log_rHy =  30),
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
    if(is.na(tmp2["log_ec",1]) )     {
        fit.beta.cis.cg <- mle2(
            minuslogl = beta_cis_new,
            #method = "CG",
            optimizer = "nlminb",
            lower     = list(log.ec = -20,log.rHy = -20),
            start     = list(log_ec = 0,log_rHy = 0),
            upper     = list(log.ec =  20,log.rHy =  20),
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

paraGetBB_old <- function(oneGeneRow,x_key,y_key,csource) { 
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
        minuslogl = beta_cis_old,
        #method="SANN",
        method = "L-BFGS-B",
        #optimizer = "nlminb",
        lower     = list(log_ec = -20,log_rHy = -20),
        start     = list(log_ec =   0,log_rHy =15),
        upper     = list(log_ec =  20,log_rHy =  30),
        data      = expCis
    )
    
    
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
    if( is.na(tmp2["log_ec",1]) )    {
        fit.beta.cis.LB <- mle2(
            minuslogl = beta_cis_old,
            method = "L-BFGS-B",
            #optimizer = "nlminb",
            lower     = list(log_ec = -20,log_rHy = -20),
            start     = list(log_ec =   0,log_rHy = 20),
            upper     = list(log_ec =  20,log_rHy =  30),
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
    if(is.na(tmp2["log_ec",1]) )     {
        fit.beta.cis.cg <- mle2(
            minuslogl = beta_cis_old,
            #method = "CG",
            optimizer = "nlminb",
            lower     = list(log.ec = -20,log.rHy = -20),
            start     = list(log_ec = 0,log_rHy = 0),
            upper     = list(log.ec =  20,log.rHy =  20),
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
        minuslogl = bino_cis_new,
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
paraGetB_old <- function(oneGeneRow,x_key,y_key,csource) {
    xHy=unlist(round(oneGeneRow[x_key]))
    nHy=unlist(round(oneGeneRow[y_key] + oneGeneRow[x_key]))
    
    C1 <- csource[x_key]
    C2 <- csource[y_key]
    
    expCis <- list(xHy=xHy,
                   nHy=nHy,
                   C1 = C1,
                   C2 = C2)
    
    fit.binom.cis <-  mle2(
        minuslogl = bino_cis_old,
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

x_key <- c(1,2,3,4,5)
y_key <- c(6,7,8,9,10)
#oneGeneRow <- c(8,22,29,44,1000,12,18,31,36,1000)
oneGeneRow <- c(10,20,30,40,10,10,20,30,40,10)
csource <- c(100,200,300,400,100,100,200,300,400,100)

paraGetBB(oneGeneRow,x_key,y_key,csource)
paraGetBB_old(oneGeneRow,x_key,y_key,csource)

paraGetB(oneGeneRow,x_key,y_key,csource)
paraGetB_old(oneGeneRow,x_key,y_key,csource)



ptoec <- function(p) {
    log2(p/(1-p))
}
