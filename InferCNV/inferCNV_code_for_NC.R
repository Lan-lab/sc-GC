# load packages
library(infercnv)
library(Seurat)
# gene_order_file
gene_order_file <- read.table(file = './gencode_v21_gen_pos.complete.txt',header = F, sep = '\t')
gene_order_file$V1 <- unlist(strsplit(as.character(as.matrix(gene_order_file$V1)),'\\|'))[seq(1, 2*nrow(gene_order_file),2)]

# create Output_Directory
out_dir <- 'PATH_TO_OUT_DIR'
if(!dir.exists(out_dir)){
  dir.create(out_dir,recursive = T)
}

# load seurat dataset of Patient GC02
load("./GC02_for_inferCNV.RData")

## prepare files for infercnv
exp <- GC02_for_inferCNV
exp@meta.data$seurat_clusters <- exp@active.ident
meta.data <- exp@meta.data
meta.data$samples <- rownames(meta.data)
new.level <- levels(meta.data$seurat_clusters)
exp <- GetAssayData(object = exp, slot = "counts")
exp <- as.data.frame(exp)
exp <- exp[,match(colnames(exp),rownames(meta.data))]
match.result <-match(rownames(exp),gene_order_file$V1)
big.index <- which(is.na(match.result)==FALSE)
small.index <- match.result[big.index]
gene_order_file.tmp <- gene_order_file[small.index, ]
rownames(gene_order_file.tmp) <- gene_order_file.tmp[,1]
gene_order_file.tmp <- gene_order_file.tmp[,-1]
colnames(gene_order_file.tmp) <- c('chr','start','stop')

# save gene_order_file
write.table(gene_order_file.tmp, file = paste0(out_dir, '/gene_order_tmp.txt'), sep = '\t', col.names = F, quote = F)

# save expression.matrix
exp <- exp[big.index, ]
expression.matrix <- exp
saveRDS(expression.matrix, file = paste0(out_dir,'/exp.rds'))

# save annotation_file
meta.data$cluster.names <- meta.data$seurat_clusters
dataframe1 <- meta.data[!meta.data$seurat_clusters %in% c("GC02_Epi"),]
dataframe2 <- meta.data[meta.data$seurat_clusters %in% c("GC02_Epi"),]
meta.data <- rbind(dataframe1, dataframe2)
annotation_file <- data.frame(samples = meta.data$samples, clusters.names = meta.data$cluster.names)
rownames(annotation_file) <- annotation_file$samples
write.table(annotation_file, file = paste0(out_dir, '/annotation_file.txt'), sep = '\t', col.names = F, row.names = F, quote = F)


# create the infercnv object
infercnv_obj = CreateInfercnvObject(raw_counts_matrix = expression.matrix,
                                    annotations_file=paste0(out_dir, '/annotation_file.txt'),
                                    gene_order_file=paste0(out_dir, '/gene_order_tmp.txt'),
                                    ref_group_names= new.level[!new.level %in% c("GC02_Epi")])
# perform infercnv operations to reveal cnv signal
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir=paste0(out_dir,'/result'), 
                             cluster_by_groups=TRUE, 
                             plot_steps=FALSE,
                             denoise=TRUE,
                             HMM=TRUE,
                             HMM_type = 'i3',
                             num_threads = 50,
                             scale_data = T
)   




