pyscenic ctx ./exp_matrix.adjacencies.tsv \
./hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather \
./hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather \
--annotations_fname ./motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
--expression_mtx_fname ./exp_matrix.csv \
--output ./exp_matrix.motifs.csv \
--num_workers 8 