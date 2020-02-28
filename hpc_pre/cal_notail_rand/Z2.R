sink('~/cloud/project/hpc_pre/cal_notail_rand/NO4A_hpc_com.txt')

i <- 1 
tmplt <- (1:20)[-13]
for(siz in 1:19) {
    for(t in 1:5) {
        choRep <- sample(tmplt,size=siz,replace = F)
        input <- c(choRep,choRep + 20)
        cat(c(formatC(i, width = 3, format = "d", flag = "0"),input,'\n'))
        i <- i + 1 
    }
}

sink()


sink('~/cloud/project/hpc_pre/cal_notail_rand/NO4A_yps_com.txt')
yps_scr <- permutations(n = 19,r=2)
tmplt <- (1:20)[-13]
i <- 1 
for(siz in 1:20) 
    for(t in 1:5) {
        set.seed(siz+t*2)
        choRep <- sample(1:nrow(yps_scr),size=siz,replace = F)
        input <- c( tmplt[yps_scr[choRep,1]],tmplt[yps_scr[choRep,2]] )
        cat(c(formatC(i, width = 3, format = "d", flag = "0"),input,'\n'))
        i <- i + 1 
    }
sink()

sink('~/cloud/project/hpc_pre/cal_notail_rand/NO4A_rm_com.txt')
rm_scr <- permutations(n = 19,r=2)
tmplt <- (1:20)[-13] + 20
i <- 1 
for(siz in 1:20) 
    for(t in 1:5) {
        set.seed(siz+t*2)
        choRep <- sample(1:nrow(rm_scr),size=siz,replace = F)
        input <- c( tmplt[rm_scr[choRep,1]],tmplt[rm_scr[choRep,2]] )
        cat(c(formatC(i, width = 3, format = "d", flag = "0"),input,'\n'))
        i <- i + 1 
    }
sink()
