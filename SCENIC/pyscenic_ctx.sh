pyscenic ctx ./Data/exp_matrix.adjacencies.tsv \
./Data/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather \
./Data/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather \
--annotations_fname ./Data/motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
--expression_mtx_fname ./Data/DC_exp.csv \
--output ./Data/exp_matrix.motifs.csv \
--num_workers 8 