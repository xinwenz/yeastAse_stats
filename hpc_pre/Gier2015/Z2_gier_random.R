library(gtools)
null44_scr <- permutations(n = 22,r=2)
sink('~/cloud/project/hpc_pre/Gier2015/rep_null_hpc_com_nostr.txt')
i <- 1 
for(siz in 1:20) 
    for(t in 1:100) {
        set.seed(siz+t*2)
        choRep <- sample(1:nrow(null44_scr),size=siz,replace = F)
        input <- c(null44_scr[choRep,1],null44_scr[choRep,2])
        cat(c(formatC(i, width = 4, format = "d", flag = "0"),input,'\n'))
        i <- i + 1 
    }
sink()
