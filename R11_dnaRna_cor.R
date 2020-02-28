load("~/cloud/project/R5_exphTrueRM14A9A_ci/exph_trueRM14A9A_esum_2000.RData")
load("~/cloud/project/R5_exphTrueRM14A9A_liko/exph_trueRM14A9A_liko_450.RData")
load("~/cloud/project/R5_exphTrueRM14A9A_liko/exph_trueRM14A9A_liko_250.RData")
load("~/cloud/project/R1_hpc_null_nostr_liko/gier_null_1000.RData")
load("~/cloud/project/R6_yps149_esum_liko/exph_yps_0500.RData")
load("~/cloud/project/R6_yps149_esum_liko/exph_yps_0490.RData")

load("~/cloud/project/R9_cal_notail_aggr_liko/cal_agrg_liko_475.RData")
load("~/cloud/project/R9_cal_yps_liko/cal_agrg_yps_liko_500.RData")
load("~/cloud/project/R5_exphTrueRM14A9A_liko/exph_trueRM14A9A_liko_450.RData")
load("~/cloud/project/R9_calh_RM4A_gene_liko/calh_notail_liko_2850.RData")


lamd12 <- rowSums(exph[,-1])
allexp <- sum(colSums(exph[,-1]))/2
lambda12 <- lamd12/allexp
theta <-  exph_trueRM14A9A_liko_450$BBf_rHy
ggplot() + geom_histogram(aes(x=theta))
o1 <- 2 ^ theta * lambda12 
ggplot() + geom_histogram(aes(x=o1)) + xlim(c(0,2))
p1 <- 2 ^ theta * lambda12 / (1 + 2 ^ theta * lambda12)
ggplot() + geom_histogram(aes(x=p1))

theta2 <- exph_yps_0500$BBf_rHy
ggplot() + geom_histogram(aes(x=theta2))
o2 <- 2^ theta2 * lambda12 
ggplot() + geom_histogram(aes(x=o2)) + xlim(c(0,10))
p2 <- 2 ^ theta2 * lambda12 / (1 + 2 ^ theta2 * lambda12)   
ggplot() + geom_histogram(aes(x=p2)) 



lamd12_44 <- rowSums(rep44[,-1])
allexp_44 <- sum(colSums(rep44[,-1]))/2
lambda12_44 <- lamd12_44/allexp_44
theta_44 <-  gier_null_2000$BBf_rHy
ggplot() + geom_histogram(aes(x=theta_44))
o_44 <- 2 ^ theta_44 * lambda12_44 
ggplot() + geom_histogram(aes(x=o_44)) + xlim(c(0,2))
p_44 <- 2 ^ theta_44 * lambda12_44 / (1 + 2 ^ theta_44 * lambda12_44)
ggplot() + geom_histogram(aes(x=p_44))

a <- which(exph_trueRM14A9A_esum_2000$ecisL > 0)

plot(exph_trueRM14A9A_esum_2000$log.ecis[a],p[a])

#mean(p)
#[1] 0.1431714
#> median(p)
#[1] 0.04808939
#> p
#######################################################
expf_snps_name <- merge.data.frame(expf,geneNu,by.x = c(1,2),by.y= c(1,2),sort=F,all=F)

useful_names <- grep("^[yr].*[HA]_rc", names(expf_snps_name), value=T)
tmp_rna <- expf_snps_name[,c("V3",useful_names)]

save(tmp_rna, file="~/cloud/project/RNAsnpsWithname.RData")

pect <- c( 0.1634535, 0.1701593, 0.1550712, 0.1687622, 0.1598212, 0.1723945,
           0.2922604, 0.1684828 ,0.1830120 ,0.1891590, 0.1693210 ,0.1587035,
           0.1768650, 0.1788209, 0.1961442, 0.1958648, 0.1606594 ,0.1930707,
           0.1399832, 0.1690416)


table(sample(-1:1, 100, replace = TRUE), sample(-1:1, 100, replace = TRUE))


tDR <- calh_notail_liko_2850$pval1 < 0.05
dim(tDR) <- c(4759,1)
rownames(tDR) <- rownames(calh_notail_liko_2850)
tDR <- cbind(tDR, exph_trueRM14A9A_liko_450[rownames(tDR),]$pval1 < 0.05)
colnames(tDR) <- colnames(c( "DNA", "RNA" ))  
cnt <- table(DNA=tDR[,1],RNA=tDR[,2])
fisher.test(cnt)
