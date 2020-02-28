load("~/cloud/project/F1_snpfilter/5_5Bsnp/cal_no_right_tail.RData") 
# get cal_no_right_tail data.frame

cal_with_label <- cal_no_right_tail[c(grep("^y.*[AH]_rc",names(cal_no_right_tail),value=T), grep("^r.*[AH]_rc",names(cal_no_right_tail),value=T))] 

x <- rowSums(cal_with_label[,c(1:20)[-13]],na.rm = TRUE)
y <- rowSums(cal_with_label[,c(1:20)[-13]+20],na.rm = TRUE)
pval_df <- data.frame(yps=x,ypsrm = x+y)
pres <- apply(pval_df, 1 , function(x) binom.test(x[1],x[2])$p.value)
pres_rank <- rank(pres)

cal_filter_10psnps <- cal_no_right_tail[pres_rank > 66352*0.12 ,] 
#hist(which(pres_rank < 66352 * 0.12)) # check the postion of those 

#ax <- rowSums(cal_filter_10psnps[,c(1:20)[-13]],na.rm = TRUE)
#ay <- rowSums(cal_filter_10psnps[,c(1:20)[-13]+20],na.rm = TRUE)
#apval_df <- data.frame(yps=ax,ypsrm = ax+ay)
#apres <- apply(apval_df, 1 , function(x) binom.test(x[1],x[2])$p.value)

yps_block <- cal_filter_10psnps[,1:4] %>% arrange(ypsChrom,ypsPosit)
add_gene_name_snpls <- function() {
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
    
yps_rm_gene <- merge.data.frame(yps_block,geneNu,by.x = c(1,2),by.y= c(1,2),sort=F,all=F)
names(yps_rm_gene)[5] <- "geneN"
    
write.table(unique(yps_rm_gene[,c(1,2)]),file="~/cloud/project/F1_snpfilter/7_5Bsnp_exp/yps128_5_snpls",row.names = F,quote=F,col.names = F)
    
write.table(unique(yps_rm_gene[,c(3,4)]),file="~/cloud/project/F1_snpfilter/7_5Bsnp_exp/rm11_B_snpls",row.names = F,quote=F,col.names = F)
    
write.table(unique(yps_rm_gene[,c(1,2,5)]),file="~/cloud/project/F1_snpfilter/7_5Bsnp_exp/yps128_5_snpls_gene",row.names = F,quote=F,col.names = F)
    
write.table(unique(yps_rm_gene[,c(3,4,5)]),file="~/cloud/project/F1_snpfilter/7_5Bsnp_exp/rm11_B_snpls_gene",row.names = F,quote=F,col.names = F)
}
add_gene_name_snpls()
