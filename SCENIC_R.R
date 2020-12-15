############################################################
## SCENIC on CPI + CD
## edward_kim@college.harvard.edu - June 2020
#############################################################
​

##############################
# 0 - Load librairies
##############################


nextflow run aertslab/SCENICprotocol \
    -profile docker \
    --loom_input expr_mat_tiny.loom \
    --loom_output pyscenic_integrated-output.loom \
    --TFs test_TFs_tiny.txt \
    --motifs motifs.tbl \
    --db *feather \
    --thr_min_genes 1
    
