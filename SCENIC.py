# bash script
wget https://resources.aertslab.org/cistarget/zsync_curl
chmod a+x zsync_curl
ZSYNC_CURL="${PWD}/zsync_curl"
echo "${ZSYNC_CURL}"


feather_database_url='https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-500bp-upstream-7species.mc9nr.feather'
"${ZSYNC_CURL}" "${feather_database_url}.zsync"


# download database ranking the whole human genome based on regulatory features
feather_database="${feather_database_url##*/}"
wget "${feather_database_url}"

# download motif annotation
wget https://resources.aertslab.org/cistarget/motif2tf/motifs-v9-nr.hgnc-m0.001-o0.0.tbl

#lastly a list of transcription factors is required for the network inference step
git clone https://github.com/aertslab/pySCENIC.git

install docker
 conda install -c conda-forge docker 

 wget https://raw.githubusercontent.com/aertslab/SCENICprotocol/master/example/test_TFs_tiny.txt
# Motif to TF annotation database:
wget https://raw.githubusercontent.com/aertslab/SCENICprotocol/master/example/motifs.tbl
# Ranking databases:
wget https://raw.githubusercontent.com/aertslab/SCENICprotocol/master/example/genome-ranking.feather
# Finally, get a tiny sample expression matrix (loom format):
wget https://raw.githubusercontent.com/aertslab/SCENICprotocol/master/example/expr_mat_tiny.loom



nextflow run aertslab/SCENICprotocol \
    -profile docker \
    --loom_input expr_mat_tiny.loom \
    --loom_output pyscenic_integrated-output.loom \
    --TFs test_TFs_tiny.txt \
    --motifs motifs.tbl \
    --db *feather \
    --thr_min_genes 1




    