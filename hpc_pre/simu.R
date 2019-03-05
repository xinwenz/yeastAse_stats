## ### 
# NB(Cr,p)
###  constant P ##### 
did_base <- sample(x = exph$`10A_yGc`[50:4200], size = 2000,replace = FALSE)
u <- did_base/sum(did_base)
N1 = 500000
N2 = N1
quantile(u*N1)

simu_gen <- function(p) {
    ans <- c()
    for(i in 1:100) {
        tmp <- rnbinom(2000, size = u * N1 * (1-p)/ p , prob = 1- p)
        ans <- cbind(ans,tmp)
    }
    ans <- data.frame(ans)
    colnames(ans) <- paste0("rep_",1:100)
    ans <- ans %>% arrange(rowSums(.))
    return(ans)
}

simu_p1_null <- simu_gen(p=0.1)
simu_p4_null <- simu_gen(p=0.4)
simu_p8_null <- simu_gen(p=0.8)
simu_pc_null <- simu_gen(p=runif(2000,min=0,max=.8))

simu_p1_depth <- colSums(simu_p1_null)
simu_p4_depth <- colSums(simu_p4_null)
simu_p8_depth <- colSums(simu_p8_null)
simu_pc_depth <- colSums(simu_pc_null)

save(list = c(ls(pattern="simu_p.*"),ls(pattern="paraGet*"),ls(pattern="neglh*"),"dbetabinom"),file = "~/cloud/project/hpc_pre/hpc_need_simu.RData" )
#############################3 command line args ################  

simu_src <- permutations(n = 100,r=2)
sink('~/cloud/project/hpc_pre/simu_hpc_com.txt')
i <- 1 
for(siz in 1:25) 
    for(t in 1:150) {
        choRep <- sample(1:nrow(simu_src),size=siz,replace = F)
        input <- c(simu_src[choRep,1],simu_src[choRep,2])
        cat(c(formatC(i, width = 4, format = "d", flag = "0"),input,'\n'))
        i <- i + 1 
    }
sink()

################### analyse data ##############
load_rep("hpc_simu_p1")
for(i in 1:3750) {
  mydataName <- paste0("simu_p1_",formatC(i,width=4,flag="0"),".RData")
  load(mydataName)
}

fpo_plot(nm="simu_p1_",eachRep = 150 ,nRep = 25 ,repGenes = 2000,note='p=0.1')


##############
load_rep("hpc_simu_p4")
for(i in 1:3750) {
  mydataName <- paste0("simu_p4_",formatC(i,width=4,flag="0"),".RData")
  load(mydataName)
}

fpo_plot(nm="simu_p4_",eachRep = 150 ,nRep = 25 ,repGenes = 2000,note='p=0.4')

############ 
load_rep("hpc_simu_p8")
for(i in 1:3750) {
  mydataName <- paste0("simu_p8_",formatC(i,width=4,flag="0"),".RData")
  load(mydataName)
}

fpo_plot(nm="simu_p8_",eachRep = 150 ,nRep = 25 ,repGenes = 2000,note='p=0.8')

########### 
load_rep("hpc_simu_p8")
for(i in 1:3750) {
  mydataName <- paste0("simu_p8_",formatC(i,width=4,flag="0"),".RData")
  load(mydataName)
}

fpo_plot(nm="simu_p8_",eachRep = 150 ,nRep = 25 ,repGenes = 2000,note='p=0.8')

##########
load_rep("hpc_simu_pc")
for(i in 1:3750) {
  mydataName <- paste0("simu_pc_",formatC(i,width=4,flag="0"),".RData")
  load(mydataName)
}

fpo_plot(nm="simu_pc_",eachRep = 150 ,nRep = 25 ,repGenes = 2000,note='p=runif(0 - 0.8)')