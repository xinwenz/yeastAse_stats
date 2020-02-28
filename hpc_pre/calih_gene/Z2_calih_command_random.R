sink('~/cloud/project/hpc_pre/calih_gene/calih_hpc_com.txt')

i <- 1 
tmplt <- (1:20)
for(siz in 1:20) {
    for(t in 1:150) {
        choRep <- sample(tmplt,size=siz,replace = F)
        input <- c(choRep,choRep + 20)
        cat(c(formatC(i, width = 4, format = "d", flag = "0"),input,'\n'))
        i <- i + 1 
    }
}

sink()

sink('~/cloud/project/hpc_pre/calih_gene/calih_RM4A_hpc_com.txt')

i <- 1 
tmplt <- (1:20)[-13]
for(siz in 1:19) {
    for(t in 1:150) {
        choRep <- sample(tmplt,size=siz,replace = F)
        input <- c(choRep,choRep + 20)
        cat(c(formatC(i, width = 4, format = "d", flag = "0"),input,'\n'))
        i <- i + 1 
    }
}

sink()