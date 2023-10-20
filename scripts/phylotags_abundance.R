#source("../scripts/phylotags_preprocessing.R")
microbio_selected$samplename <- gsub("^\\s+|\\s+$", "", rownames(microbio_selected)) 

# 16S rRNA
ssu_normcounts <- counts(ssu_dds.wald,normalized=TRUE) %>% 
  as.data.frame()  %>%
  rownames_to_column(var="asv") %>% 
  as_tibble()

ssu_sig_asv.nw_ob_norm <- ssu_normcounts %>%
  filter(asv %in% rownames(ssu_sigtab.nw_ob)) 

ssu_gathered_sigASV_nw_ob <- ssu_sig_asv.nw_ob_norm %>%
  gather(colnames(ssu_sig_asv.nw_ob_norm)[2:115], key = "samplename", value = "normalized_counts")
ssu_gathered_sigASV_nw_ob <- inner_join(microbio_selected %>% 
                                           mutate(samplename=rownames(microbio_selected)),
                                         ssu_gathered_sigASV_nw_ob, by="samplename")#inner_join(microbio_selected,ssu_gathered_sigASV_nw_ob)


# dnaK  
dnaK_normcounts <- counts(dnaK_dds.wald,normalized=TRUE) %>% 
  as.data.frame()  %>%
  rownames_to_column(var="asv") %>% 
  as_tibble()

dnaK_sig_asv.nw_ob_norm <- dnaK_normcounts %>%
  filter(asv %in% rownames(dnaK_sigtab.nw_ob)) 

dnaK_gathered_sigASV_nw_ob <- dnaK_sig_asv.nw_ob_norm %>%
  gather(colnames(dnaK_sig_asv.nw_ob_norm)[2:115], key = "samplename", value = "normalized_counts")
dnaK_gathered_sigASV_nw_ob <- inner_join(microbio_selected %>% 
                                           mutate(samplename=rownames(microbio_selected)),
                                         dnaK_gathered_sigASV_nw_ob, by="samplename")#inner_join(microbio_selected,dnaK_gathered_sigASV_nw_ob)

# gyrB
gyrB_normcounts <- counts(gyrB_dds.wald,normalized=TRUE) %>% 
  as.data.frame()  %>%
  rownames_to_column(var="asv") %>% 
  as_tibble()

gyrB_sig_asv.nw_ob_norm <- gyrB_normcounts %>%
  filter(asv %in% rownames(gyrB_sigtab.nw_ob)) 

gyrB_gathered_sigASV_nw_ob <- gyrB_sig_asv.nw_ob_norm %>%
  gather(colnames(gyrB_sig_asv.nw_ob_norm)[2:81], key = "samplename", value = "normalized_counts")
gyrB_gathered_sigASV_nw_ob <- inner_join(microbio_selected %>% 
                                           mutate(samplename=rownames(microbio_selected)),
                                         gyrB_gathered_sigASV_nw_ob, by="samplename")#inner_join(microbio_selected,dnaK_gathered_sigASV_nw_ob)


############
# PLOTS ----
############

#ordenar de acuerdo a foldChanges
#dnaK_sig_asv.nw_ow_norm <- dnaK_sig_asv.nw_ow_norm[match(rownames(dnaK_sigtab.nw_ow[order(dnaK_sigtab.nw_ow$log2FoldChange),]),dnaK_sig_asv.nw_ow_norm$asv),]


#dnaK_sig_asv.nw_ob_norm <- dnaK_sig_asv.nw_ob_norm[match(rownames(dnaK_sigtab.nw_ob[order(dnaK_sigtab.nw_ob$log2FoldChange),]),dnaK_sig_asv.nw_ob_norm$asv),]
#dnaK_gathered_sigASV_nw_ob <- inner_join(microbio_selected,dnaK_gathered_sigASV_nw_ob)

#Plot and for labels order according to Foldchange and sepcies names
#dnaK_gathered_sigASV_nw_ob$asv <- factor(dnaK_gathered_sigASV_nw_ob$asv,levels=rownames(dnaK_sigtab.nw_ob[order(dnaK_sigtab.nw_ob$log2FoldChange),]))

#label with species and trimm the prefix

#plot_labels[c(2,3)]<- "CAG-83"
#plot_labels[9]<-"Gemmiger.sp" 
#pdf("/lsalazar/proyectos/phylotags/reports/DA_lean-ow_samples.pdf",height=5.5,width=5.5)
#ggplot(dnaK_gathered_sigASV_nw_ob, aes(x=normalized_counts , y = reorder(asv,normalized_counts), label = samplename, color=bmi_class)) + 
#  geom_point(size = 3, alpha = 0.8) +
#  scale_color_manual(values=c("#FFF600", "#FE5A1D", "#FF0000"))+
#  geom_line(size = 1.2, aes(group = asv, color = bmi_class)) +
#  geom_text_repel(data          = subset(dnaK_gathered_sigASV_nw_ow, normalized_counts > 10),
#                  size          = 2,
#                  box.padding   = 0.35,
#                  point.padding = 0.5,
                  #force         = 50,
#                  segment.size  = 0.2,
#                  segment.color = "grey50",
#                  direction     = "y")+ 
#  scale_x_log10() +
#  scale_y_discrete(labels = plot_labels) +
#  xlab("log10 Normalized Counts") +
#  ylab("ASVs") +
#  ggtitle("Significant Differential Abundant ASVs Lean vs Overweight  ") +
#  theme_bw() +
#  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#  theme(plot.title = element_text(hjust = 0.5))

show.DEgenes<-function(ASV=NULL, tdata=NULL, cov=NULL, cg=NULL, q=NULL){#tdata=sp.dge.keep$counts,cov=sp.targets, q=NULL){ 
  gene<-which(rownames(tdata)==ASVid) 
  y<-as.numeric(tdata[gene,])
  y<-log2(y+1) #for log2 count data
  interaction.plot(droplevels(cov$line),droplevels(cov$treatment),y,
                   ylim=c(min(y),max(y)),
                   ylab="Log2(Counts Per Million)",
                   xlab="Line",
                   legend=FALSE,
                   col=c("blue","red"))
  points(rep(1:4, each=6)[cov$treatment=="par"],y[cov$treatment=="par"],pch=19,col="red")
  points(rep(1:4, each=6)[cov$treatment=="ctl"],y[cov$treatment=="ctl"],pch=17,col="blue") 
  title(main=paste(cg,"(", FBid, ")", sep=" "), sub=paste("FDR=", signif(q , digits=4)))
}

# INTERACTION PLOT
#https://ggplot2tutor.com/tutorials/interaction_plot
# dnaK
plot_labels_dnaK = sapply(strsplit(as.character(dnaK_sigtab.nw_ob[rownames(dnaK_sigtab.nw_ob[order(dnaK_sigtab.nw_ob$log2FoldChange),]),]$species), 
                              "\\__"), `[`, 2)

fig5B <- ggplot(dnaK_gathered_sigASV_nw_ob, aes(asv, log10(normalized_counts +1))) +
  geom_jitter(aes(color=bmi_class), position=position_jitter(0.2), size = 2) +
  scale_color_manual(values =  c("#00AFBB", "#E7B800", "#FC4E07"))+
  ylab("Log10 Normalized counts") +
  xlab("") +
  scale_x_discrete(labels = plot_labels_dnaK) +
  ggtitle("dnaK") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust = 1), 
        plot.margin = margin(0,0,0,1,"cm")) +
  guides(color=guide_legend(title="BMI class"))
dev.off()

# gyrB
plot_labels_gyrB = sapply(strsplit(as.character(gyrB_sigtab.nw_ob[rownames(gyrB_sigtab.nw_ob[order(gyrB_sigtab.nw_ob$log2FoldChange),]),]$species),
                                   "\\__"), `[`, 2)

gyrB_gathered_sigASV_nw_ob$asv <- factor(gyrB_gathered_sigASV_nw_ob$asv, levels=rownames(gyrB_sigtab.nw_ob[order(gyrB_sigtab.nw_ob$log2FoldChange),]))

fig5A <- ggplot(gyrB_gathered_sigASV_nw_ob, aes(asv, log10(normalized_counts +1))) +
  geom_jitter(aes(color=bmi_class), position=position_jitter(0.2), size = 2) +
  scale_color_manual(values =  c("#00AFBB", "#E7B800", "#FC4E07"))+
  ylab("Log10 Normalized counts") +
  xlab("") +
  scale_x_discrete(labels = plot_labels_gyrB) +
  ggtitle("gyrB") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust = 1), 
        plot.margin = margin(0,0,0,1,"cm")) +
  guides(color=guide_legend(title="BMI class"))

# Fig5 paper
plot_grid(fig5A, fig5B, ncol=1, nrow=2, labels="AUTO")
