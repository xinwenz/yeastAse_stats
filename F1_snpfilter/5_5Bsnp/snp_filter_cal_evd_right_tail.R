load("~/cloud/project/F1_snpfilter/3_5Bsnp/cal_evd.RData")
attach(cal_evd)

ckdis_len <- function(x) {
    t1 <- paste0("y",x,"_rc")
    t2 <- paste0("r",x, "_rc")

    p <<- ggplot(cal_evd) +
    geom_histogram(aes(get(t1)), binwidth = 1, alpha=0.2, fill="red") +
    #geom_histogram(aes(y10A_ac), binwidth = 1, alpha=0.2, fill = "blue") + 
    geom_histogram(aes(get(t2)), binwidth = 1, alpha=0.2, fill="yellow") +
    #geom_histogram(aes(r10A_ac), binwidth = 1, alpha = 0.2, fill= "green") + 
    scale_x_continuous(breaks =seq(0,100,by=5)) 
    #ylim(0,7000)
    
    #print(p)
    n1 <<- which( rank(get(t1), ties.method = "min") > nrow(cal_evd) * 0.9975  )
    n2 <<- which( rank(get(t2), ties.method = "min") > nrow(cal_evd) * 0.9975  )
    n12 <<- union(n1,n2)
    
    #print(p)
    #print(length(n1))
    #print(length(n2))
    #print(paste("obsv_vs_expect:", length(intersect(n1,n2)), 0.0025 * 0.0025 * nrow(cal_evd) ) )
    #print(min(get(t1)[n1]))
    #print(min(get(t2)[n2]))
    
}


ckdis_len("10A")
ckdis_len("10H")
ckdis_len("11A")
ckdis_len("11H")
ckdis_len("13A")
ckdis_len("13H")
ckdis_len("14A")
ckdis_len("14H")
ckdis_len("15A")
ckdis_len("15H")
ckdis_len("3A")
ckdis_len("3H")
ckdis_len("4A") # too much MA, median small than 5
ckdis_len("4H")
ckdis_len("5A")
ckdis_len("5H")
ckdis_len("6A")
ckdis_len("6H")
ckdis_len("9A")
ckdis_len("9H")

items <- paste0(rep(c(10,11,13,14,15,3,4,5,6,9),each=2),rep(c("A","H"),10))
right_tail_n <- c()
for(i in items) {
    if(i == "4A") {next}
    ckdis_len(i)
    right_tail_n <- union(right_tail_n,n12)
}

cal_no_right_tail <- cal_evd[-right_tail_n, ]
save(cal_no_right_tail,file="~/cloud/project/F1_snpfilter/5_5Bsnp/cal_no_right_tail.RData")
################################### test the genome-wide coverage ########## 
maindf <- cal_no_right_tail  # %>% filter(ypsChrom != "Seg23" & ypsChrom != "Seg6"  ) 
windsiz <- 50
maxlab <- floor(nrow(maindf)/windsiz)
remnum <- nrow(maindf) - maxlab * windsiz  
segL <- c ( rep(1: maxlab , each= windsiz ), rep(maxlab+1 , remnum ) )

maindf_segL <- cbind(segL, maindf)
expnames <- grep("^[yr].*[HA]_rc", names(maindf), value=T)

maindf_segL_res <- maindf_segL %>% group_by(segL) %>% summarize_at(.vars = expnames , .funs=c(ms="sum"), na.rm = T )

maindf_ypsChrom_res <- maindf_segL %>% group_by(segL) %>% summarize(ypsChromm = min(ypsChrom))

maindf_ypsChrom_segL_res <- merge.data.frame(maindf_segL_res,maindf_ypsChrom_res)

yrn <- grep("^y.*[HA]_rc_ms",names(maindf_ypsChrom_segL_res),value=T)
#ya <- colSums(cal[,grep("^y.*[HA]_ac",names(cal),value=T)],na.rm=T)
#yra <- data.frame(x=yr,y=ya, type="YPS128")

rrn <- grep("^r.*[HA]_rc_ms",names(maindf_ypsChrom_segL_res),value=T)


for (i in 1:20) {
    pp <- ggplot(maindf_ypsChrom_segL_res) + 
        geom_bar(aes(x=segL, y= get(yrn[i]), color=ypsChromm, fill=ypsChromm), stat = "identity") +
        geom_bar(aes(x=segL, y= (-1) * get(rrn[i]), color=ypsChromm, fill=ypsChromm), stat = "identity") +
        scale_color_viridis(discrete = TRUE, option = "C") +
        scale_fill_viridis(discrete = TRUE, option = "C") + 
        ggtitle(yrn[i])
    print(i)
    print(pp)
}















yps_block <- cal_no_right_tail[,1:4] %>% arrange(ypsChrom,ypsPosit)


genN <- read.table("~/cloud/project/F1_snpfilter/yps128_5_snp_gene.txt",header = F)
overlap_gene <- c()
    for(i in 2:nrow(genN)) {
        if(genN[i,2] == genN[i-1,2] ){
            print(unname(genN[i-1,]))
            print(unname(genN[i,]))
            overlap_gene <- c(overlap_gene,i-1,i)
        }
    }
geneNu <- genN[-overlap_gene,]


yps_rm_notail_gene <- merge.data.frame(yps_block,geneNu,by.x = c(1,2),by.y= c(1,2),sort=F,all=F)
names(yps_rm_notail_gene)[5] <- "geneN"
    
    
    write.table(yps_rm_notail_gene,file="~/cloud/project/F1_snpfilter/5_5Bsnp/yr_gene_notail",row.names = F,quote=F,col.names = F)
    
    write.table(unique(yps_rm_notail_gene[,c(1,2)]),file="~/cloud/project/F1_snpfilter/5_5Bsnp/yps128_5_snpls_notail",row.names = F,quote=F,col.names = F)
    
    write.table(unique(yps_rm_notail_gene[,c(3,4)]),file="~/cloud/project/F1_snpfilter/5_5Bsnp/rm11_B_snpls_notail",row.names = F,quote=F,col.names = F)
    
    write.table(unique(yps_rm_notail_gene[,c(1,2,5)]),file="~/cloud/project/F1_snpfilter/5_5Bsnp/yps128_5_snpls_notail_gene",row.names = F,quote=F,col.names = F)
    
    write.table(unique(yps_rm_notail_gene[,c(3,4,5)]),file="~/cloud/project/F1_snpfilter/5_5Bsnp/rm11_B_snpls_notail_gene",row.names = F,quote=F,col.names = F)
    