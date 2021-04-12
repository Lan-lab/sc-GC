library(dplyr)
library(ggplot2)

#######read the celllphoneDB output results
all_means = read.table("cellphoneDB/out/means.txt", header=T, stringsAsFactors = F, sep="\t", comment.char = '', check.names=F)
all_pval=read.table("cellphoneDB/out/pvalues.txt",sep = "\t",header = T, stringsAsFactors = F, comment.char = '', check.names = F)

######read the gene pairs that need to display
gene_pair<-read.table("dotplot_gene.txt",header = F,sep = "\t",as.is = T)
gene_pair<-as.character(gene_pair$V1)

selected_columns<-c("Mφ_APOE|Fib_1","Mφ_APOE|Endo_1","Mφ_APOE|SMC_1",
                    "DC_LAMP3|Fib_1","DC_LAMP3|Endo_1","DC_LAMP3|SMC_1",
                    "pDC_LILRA4|Fib_1","pDC_LILRA4|Endo_1","pDC_LILRA4|SMC_1",
                    "Tumor|Fib_1","Tumor|Endo_1","Tumor|SMC_1","Tumor|Tumor",
                    "Fib_1|Mφ_APOE","Fib_1|DC_LAMP3","Fib_1|pDC_LILRA4","Fib_1|Tumor",
                    "Endo_1|Mφ_APOE","Endo_1|DC_LAMP3","Endo_1|pDC_LILRA4","Endo_1|Tumor",
                    "SMC_1|Mφ_APOE","SMC_1|DC_LAMP3","SMC_1|pDC_LILRA4","SMC_1|Tumor")


sel_pval = all_pval[match(unique(gene_pair),all_pval$interacting_pair), selected_columns]
sel_means = all_means[match(unique(gene_pair),all_means$interacting_pair), selected_columns]

df_names = expand.grid(unique(gene_pair), selected_columns)
pval = unlist(sel_pval)
pval[pval==0] = 0.0009
plot.data = cbind(df_names,pval)
pr = unlist(as.data.frame(sel_means))
pr[pr==0] = 1
plot.data = cbind(plot.data,log2(pr))
colnames(plot.data) = c('pair', 'clusters', 'pvalue', 'mean')
my_palette <- colorRampPalette(c("black", "blue", "yellow", "red"), alpha=TRUE)(n=399)
my_palette = colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(50)[-c(25:35)]
#my_palette <- rev(colorRampPalette(brewer.pal(11,"RdBu"))(n=399))

pdf(file = "dotplot.pdf",width = 15,height = 20)
ggplot(plot.data,aes(x=clusters,y=pair)) +
  geom_point(aes(size=-log10(pvalue),color=mean)) +
  scale_color_gradientn('Log2 mean (Molecule 1, Molecule 2)', colors=my_palette) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),legend.position="none",
        axis.text=element_text(size=14, colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, size=5,colour = "black"),
        axis.text.y = element_text(size=5, colour = "black"),
        axis.title=element_blank(),
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))+
  scale_y_discrete(limits = rev(unique(gene_pair)))
dev.off()


###################heatmap##############
allcell<- readRDS(file = "R_all_cells.rds")
cluster_ave <- AverageExpression(allcell, return.seurat = TRUE,assays = "RNA",slot = "data")

######read the gene pairs that need to display
list1<-read.table("heatmap_gene.txt",header = F,sep = "\t",as.is = T)
clu1<-c("Fib_1","Endo_1","SMC_1")
clu2<-c("Endo_1","Fib_1","SMC_1")


input<- subset(cluster_ave,idents=clu1)
pdf(file = "heatmap_1.pdf",width = 15,height = 5)
pheatmap(as.matrix(input@assays$RNA@scale.data)[list1$V1,clu1], show_rownames = T, show_colnames = F, 
         scale="row", cluster_cols = F,
        cluster_rows = F, fontsize_row = 12, fontsize_col = 12,
         family = "Arial",
         color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(20)[-(10:12)])
dev.off()

input<- subset(cluster_ave,idents=clu2)
pdf(file = "heatmap_2.pdf",width = 15,height = 5)
pheatmap(as.matrix(input@assays$RNA@scale.data)[list1$V2,clu2], show_rownames = T, show_colnames = F, 
         scale="row", cluster_cols = F,
         cluster_rows = F, fontsize_row = 12, fontsize_col = 12,
         family = "Arial",
         color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(20)[-(10:12)])
dev.off()