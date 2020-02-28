# rm11_seg_length
# from HPC  :  bioawk -c fastx '{printf length($seq); printf "," }' rm11_B.fasta

rm11_seg <- c(1483673,1024423,1056886,1075814,920235,730198,801060,552808,431223,118427,256308,207699,26410,919793,788098,38684,669304,534885,311180)

yps128_seg <- c(1495194,1082331,1073056,928731,26007,25770,932551,808601,536555,576146,441150,319342,1051,1140,1215,1096,1428,1862,1173,1528,289440,238034,87421,572538,36647,454265,766963,667599,723488)

scer64_seg <- c(230218,813184,316620,1531933,576874,270161,1090940,562643,439888,745751,666816,1078177,924431,784333,1091291,948066,85779)
    
    
rm11_seg <- sort(rm11_seg,decreasing = T)
b <- data.frame(rank=1:length(rm11_seg), cum = cumsum(rm11_seg)/10^6, genome = "rm11")

yps128_seg <- sort(yps128_seg,decreasing = T)
a <- data.frame(rank=1:length(yps128_seg), cum = cumsum(yps128_seg)/10^6, genome = "yps128")

scer64_seg <- sort(scer64_seg,decreasing = T)
c <- data.frame(rank=1:length(scer64_seg), cum = cumsum(scer64_seg)/10^6, genome = "scer64")

ans <- rbind(a,b,c)

library(ggplot2)
ggplot(ans,aes(x=rank, y= cum, linetype=genome)) + 
  geom_point(aes(color=genome)) + 
  geom_line(aes(color=genome), lwd = 1) + 
  scale_linetype_manual(values=c("solid","solid",'dashed')) + 
  scale_color_manual(values=c('cyan', 'orange', 'black')) + 
  labs(title="Assembly Contiguity", y="cumulative assembly length (Mbp)", x="contig length rank") + 
  theme_bw() +
  theme(axis.text.x=element_text(angle=0, hjust=1)) +
  #coord_cartesian(xlim=c(1,20), ylim=c(0, 0.8)) +
  #scale_x_discrete(labels=1:20) + 
  #scale_y_continuous(breaks = seq(0,0.8,by=0.05)) + 
  theme(
    panel.grid = element_line(color="grey85"),  
    legend.position = c(0.875, 0.5),
    legend.background =  element_rect( size=0.5, linetype="solid",colour ="darkblue"),
    legend.title = element_blank(),
    legend.text = element_text(size=14,face="plain"),
    plot.title = element_text(size = 18, face = "bold",hjust=.5),
    axis.text.x = element_text(colour="grey20",size=13,hjust=.5,vjust=.5,face="plain"),
    axis.text.y = element_text(colour="grey20",size=13,hjust=.5,vjust=.5,face="plain"),
    axis.title.x = element_text(colour="grey20",size=15,hjust=.5,vjust=0,face="bold"), 
    axis.title.y = element_text(colour="grey20",size=15,hjust=.5,vjust=.5,face="bold"))

ggsave("fig_contig.pdf", device = "pdf", width = 15, height = 15, units = "cm",dpi=300)

       