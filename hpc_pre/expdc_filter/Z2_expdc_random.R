sink('~/cloud/project/hpc_pre/expdc_filter/expdc_hybridRM14A9A_hpc_com.txt')
i <- 1 
tmplt <- (1:20)[-c(7,19)]
for(siz in 1:18) {
    for(t in 1:20) {
        choRep <- sample(tmplt,size=siz,replace = F)
        input <- c(choRep,choRep + 20)
        cat(c(formatC(i, width = 3, format = "d", flag = "0"),input,'\n'))
        i <- i + 1 
    }
}
sink()


yps_scr <- permutations(n = 18,r=2)
sink('~/cloud/project/hpc_pre/expdc_filter/exph_ypsRM14A9A_hpc_com.txt')
i <- 1 
for(siz in 1:20) 
    for(t in 1:25) {
        set.seed(siz+t*2)
        choRep <- sample(1:nrow(yps_scr),size=siz,replace = F)
        input <- c( tmplt[yps_scr[choRep,1]],tmplt[yps_scr[choRep,2]] )
        cat(c(formatC(i, width = 4, format = "d", flag = "0"),input,'\n'))
        i <- i + 1 
    }
sink()
