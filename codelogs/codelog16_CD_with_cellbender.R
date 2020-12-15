############################################################
## Computational Integration Analysis
## edward_kim@college.harvard.edu - December 2020
#############################################################

##############################
# 0 - Load librairies
##############################
library(dplyr)
library(ggplot2)
library(Seurat)
library(future)
library(future.apply)
#library(fgsea)
#library(harmony)
#library(scran)
#library(pheatmap)

##############################
# 1 - Definitions & Settings
##############################
'%ni%' = Negate('%in%')
plan("multicore", workers = 10)
options(future.globals.maxSize = 20000 * 1024^2)

############################## 
# 2 - Source file 
##############################
source("scripts/Stat111functions.R")
source("scripts/new_clustering_functions.R")

##############################
# eval barcode rank plots
##############################
dirnames <- dir("data/CD_martin",pattern="GSM")

## Key functions
# output is a list of Seurat objects
load_data <- function(dirname) {
	s <- Read10X(data.dir = paste("data/CD_martin",dirname,sep="/"))
	s1 <- CreateSeuratObject(counts = s, project = dirname, min.cells = 0, min.features = 0)
	s1
}
# output is a saved barcode rank plot
plot_brp <- function(nUMIs,lower.rank,upper.rank, title) {
	line1 <- log10(lower.rank)
	line2 <- log10(upper.rank)
	ggplot(data.frame(nUMI=sort(log10(nUMIs), decreasing=T), rank = log10(1:length(nUMIs)))) + geom_line(aes(x=rank,y=nUMI)) + 
		geom_vline(aes(xintercept=line1),color="red",linetype="longdash") + 
		geom_text(aes(x=line1, label=paste("rank",10^line1,"\n"), y=0.15), colour="red", angle=90) +
		geom_vline(aes(xintercept=line2),color="orange",linetype="longdash") + 
		geom_text(aes(x=line2, label=paste("\nrank",10^line2), y=0.15), colour="orange", angle=90) +
		ggtitle(title) + xlab("log rank") + ylab("log nUMI") + theme_classic()
	ggsave(paste("figures/CD_martin_cb/",title,"brp.pdf", sep=""), width= 7, height= 7, units= "in")
}

# first 5
s <- lapply(dirnames[1:5], load_data)
plot_brp(s[[1]]$nCount_RNA, 1900, 10000, dirnames[1]) # 69 run
plot_brp(s[[2]]$nCount_RNA, 6000, 20000, dirnames[2]) # 68 run
plot_brp(s[[3]]$nCount_RNA, 1800, 10000, dirnames[3]) # 122 run
plot_brp(s[[4]]$nCount_RNA, 3500, 10000, dirnames[4]) # 123 run
plot_brp(s[[5]]$nCount_RNA, 5000, 20000, dirnames[5]) # 128 run
# 6-10
s <- lapply(dirnames[6:10], load_data)
plot_brp(s[[1]]$nCount_RNA, 1000, 10000, dirnames[6]) # 129 run
plot_brp(s[[2]]$nCount_RNA, 4000, 20000, dirnames[7]) # 135 run
plot_brp(s[[3]]$nCount_RNA, 3000, 20000, dirnames[8]) # 138 run v2
plot_brp(s[[4]]$nCount_RNA, 5000, 20000, dirnames[9]) # 158 run
plot_brp(s[[5]]$nCount_RNA, 5000, 30000, dirnames[10]) # 159 run
# 11-15
s <- lapply(dirnames[11:15], load_data)
plot_brp(s[[1]]$nCount_RNA, 2000, 15000, dirnames[11]) # 180 run v2
plot_brp(s[[2]]$nCount_RNA, 7000, 20000, dirnames[12]) # 181 run
plot_brp(s[[3]]$nCount_RNA, 3000, 15000, dirnames[13]) # 186 run v2
plot_brp(s[[4]]$nCount_RNA, 2000, 7000, dirnames[14]) # 187 run v2 (test convergence?)
plot_brp(s[[5]]$nCount_RNA, 4000, 30000, dirnames[15]) # 189 run v2
# 16-22
s <- lapply(dirnames[16:22], load_data)
plot_brp(s[[1]]$nCount_RNA, 1000, 10000, dirnames[16]) # 190 run
plot_brp(s[[2]]$nCount_RNA, 5000, 30000, dirnames[17]) # 192 run v2
plot_brp(s[[3]]$nCount_RNA, 5000, 30000, dirnames[18]) # 193 run v2
plot_brp(s[[4]]$nCount_RNA, 4000, 15000, dirnames[19]) # 195 run
plot_brp(s[[5]]$nCount_RNA, 4000, 50000, dirnames[20]) # 196 run
plot_brp(s[[6]]$nCount_RNA, 800, 7000, dirnames[21]) # 208 run
plot_brp(s[[7]]$nCount_RNA, 900, 10000, dirnames[22]) # 209 run


#### The following code is in bash script ####
##################################
#### Cell Bender Installation ####
##################################
conda create -n cellbender2 python=3.7
source activate cellbender2
conda install -c anaconda pytables
conda install pytorch torchvision torchaudio cudatoolkit=10.1 -c pytorch
conda install pytorch torchvision cudatoolkit=10.1 -c pytorch
git clone https://github.com/broadinstitute/CellBender.git
pip install -e CellBender

##########################
# Cell Bender Sample 69  #
##########################
CUDA_VISIBLE_DEVICES=7 cellbender remove-background \
                 --input /home/esk17/1.projects/colitis/data/CD_martin/GSM3972009_69/ \
                 --output /home/esk17/1.projects/colitis/data/CD_martin_cellbender/GSM3972009_69v2cb.h5 \
                 --cuda \
                 --expected-cells 1900 \
                 --total-droplets-included 10000 \
                 --fpr 0.01 \
                 --epochs 150

##########################
# Cell Bender Sample 68  #
##########################
CUDA_VISIBLE_DEVICES=6 cellbender remove-background \
                 --input /home/esk17/1.projects/colitis/data/CD_martin/GSM3972010_68/ \
                 --output /home/esk17/1.projects/colitis/data/CD_martin_cellbender/GSM3972010_68cb.h5 \
                 --cuda \
                 --expected-cells 6000 \
                 --total-droplets-included 20000 \
                 --fpr 0.01 \
                 --epochs 150

##########################
# Cell Bender Sample 122 #
##########################
CUDA_VISIBLE_DEVICES=5 cellbender remove-background \
                 --input /home/esk17/1.projects/colitis/data/CD_martin/GSM3972011_122/ \
                 --output /home/esk17/1.projects/colitis/data/CD_martin_cellbender/GSM3972011_122cb.h5 \
                 --cuda \
                 --low-count-threshold 5\
                 --expected-cells 1800 \
                 --total-droplets-included 10000 \
                 --fpr 0.01 \
                 --epochs 150

##########################
# Cell Bender Sample 123 #
##########################
CUDA_VISIBLE_DEVICES=7 cellbender remove-background \
                 --input /home/esk17/1.projects/colitis/data/CD_martin/GSM3972012_123/ \
                 --output /home/esk17/1.projects/colitis/data/CD_martin_cellbender/GSM3972012_123cb.h5 \
                 --cuda \
                 --expected-cells 3500 \
                 --total-droplets-included 10000 \
                 --fpr 0.01 \
                 --epochs 150

##########################
# Cell Bender Sample 128 #
##########################
CUDA_VISIBLE_DEVICES=2 cellbender remove-background \
                 --input /home/esk17/1.projects/colitis/data/CD_martin/GSM3972013_128/ \
                 --output /home/esk17/1.projects/colitis/data/CD_martin_cellbender/GSM3972013_128cb.h5 \
                 --cuda \
                 --low-count-threshold 5\
                 --expected-cells 5000 \
                 --total-droplets-included 20000 \
                 --fpr 0.01 \
                 --epochs 150

##########################
# Cell Bender Sample 129 #
##########################
CUDA_VISIBLE_DEVICES=5 cellbender remove-background \
                 --input /home/esk17/1.projects/colitis/data/CD_martin/GSM3972014_129/ \
                 --output /home/esk17/1.projects/colitis/data/CD_martin_cellbender/GSM3972014_129cb.h5 \
                 --cuda \
                 --expected-cells 1000 \
                 --total-droplets-included 10000 \
                 --fpr 0.01 \
                 --epochs 150

##########################
# Cell Bender Sample 135 #
##########################
CUDA_VISIBLE_DEVICES=3 cellbender remove-background \
                 --input /home/esk17/1.projects/colitis/data/CD_martin/GSM3972015_135/ \
                 --output /home/esk17/1.projects/colitis/data/CD_martin_cellbender/GSM3972015_135cb.h5 \
                 --cuda \
                 --low-count-threshold 10\
                 --expected-cells 4000 \
                 --total-droplets-included 20000 \
                 --fpr 0.01 \
                 --epochs 150

##########################
# Cell Bender Sample 138 #
##########################
CUDA_VISIBLE_DEVICES=7 cellbender remove-background \
                 --input /home/esk17/1.projects/colitis/data/CD_martin/GSM3972016_138/ \
                 --output /home/esk17/1.projects/colitis/data/CD_martin_cellbender/GSM3972016_138cb.h5 \
                 --cuda \
                 --expected-cells 1000 \
                 --total-droplets-included 7000 \
                 --fpr 0.01 \
                 --epochs 150

##########################
# Cell Bender Sample 138v2 #
##########################
CUDA_VISIBLE_DEVICES=3 cellbender remove-background \
                 --input /home/esk17/1.projects/colitis/data/CD_martin/GSM3972016_138/ \
                 --output /home/esk17/1.projects/colitis/data/CD_martin_cellbender/GSM3972016_138v2cb.h5 \
                 --cuda \
                 --low-count-threshold 10\
                 --expected-cells 3000 \
                 --total-droplets-included 20000 \
                 --fpr 0.01 \
                 --epochs 150

##########################
# Cell Bender Sample 158 #
##########################
CUDA_VISIBLE_DEVICES=7 cellbender remove-background \
                 --input /home/esk17/1.projects/colitis/data/CD_martin/GSM3972017_158/ \
                 --output /home/esk17/1.projects/colitis/data/CD_martin_cellbender/GSM3972017_158cb.h5 \
                 --cuda \
                 --low-count-threshold 10\
                 --expected-cells 5000 \
                 --total-droplets-included 20000 \
                 --fpr 0.01 \
                 --epochs 150

##########################
# Cell Bender Sample 159 #
##########################
CUDA_VISIBLE_DEVICES=6 cellbender remove-background \
                 --input /home/esk17/1.projects/colitis/data/CD_martin/GSM3972018_159/ \
                 --output /home/esk17/1.projects/colitis/data/CD_martin_cellbender/GSM3972018_159cb.h5 \
                 --cuda \
                 --low-count-threshold 10\
                 --expected-cells 5000 \
                 --total-droplets-included 30000 \
                 --fpr 0.01 \
                 --epochs 150

##########################
# Cell Bender Sample 180 #
##########################
CUDA_VISIBLE_DEVICES=4 cellbender remove-background \
                 --input /home/esk17/1.projects/colitis/data/CD_martin/GSM3972019_180/ \
                 --output /home/esk17/1.projects/colitis/data/CD_martin_cellbender/GSM3972019_180cb.h5 \
                 --cuda \
                 --expected-cells 1000 \
                 --total-droplets-included 5000 \
                 --fpr 0.01 \
                 --epochs 150

##########################
# Cell Bender Sample 180v2 #
##########################
CUDA_VISIBLE_DEVICES=7 cellbender remove-background \
                 --input /home/esk17/1.projects/colitis/data/CD_martin/GSM3972019_180/ \
                 --output /home/esk17/1.projects/colitis/data/CD_martin_cellbender/GSM3972019_180v2cb.h5 \
                 --cuda \
                 --low-count-threshold 10\
                 --expected-cells 2000 \
                 --total-droplets-included 15000 \
                 --fpr 0.01 \
                 --epochs 150

##########################
# Cell Bender Sample 181 #
##########################
CUDA_VISIBLE_DEVICES=1 cellbender remove-background \
                 --input /home/esk17/1.projects/colitis/data/CD_martin/GSM3972020_181/ \
                 --output /home/esk17/1.projects/colitis/data/CD_martin_cellbender/GSM3972020_181cb.h5 \
                 --cuda \
                 --low-count-threshold 10\
                 --expected-cells 7000 \
                 --total-droplets-included 20000 \
                 --fpr 0.01 \
                 --epochs 150

##########################
# Cell Bender Sample 186 #
##########################
CUDA_VISIBLE_DEVICES=2 cellbender remove-background \
                 --input /home/esk17/1.projects/colitis/data/CD_martin/GSM3972021_186/ \
                 --output /home/esk17/1.projects/colitis/data/CD_martin_cellbender/GSM3972021_186cb.h5 \
                 --cuda \
                 --expected-cells 1000 \
                 --total-droplets-included 5000 \
                 --fpr 0.01 \
                 --epochs 150

############################
# Cell Bender Sample 186v2 #
############################
CUDA_VISIBLE_DEVICES=1 cellbender remove-background \
                 --input /home/esk17/1.projects/colitis/data/CD_martin/GSM3972021_186/ \
                 --output /home/esk17/1.projects/colitis/data/CD_martin_cellbender/GSM3972021_186v2cb.h5 \
                 --cuda \
                 --low-count-threshold 5\
                 --expected-cells 3000 \
                 --total-droplets-included 15000 \
                 --fpr 0.01 \
                 --epochs 150

##########################
# Cell Bender Sample 187 #
##########################
CUDA_VISIBLE_DEVICES=3 cellbender remove-background \
                 --input /home/esk17/1.projects/colitis/data/CD_martin/GSM3972022_187/ \
                 --output /home/esk17/1.projects/colitis/data/CD_martin_cellbender/GSM3972022_187cb.h5 \
                 --cuda \
                 --expected-cells 1000 \
                 --total-droplets-included 5000 \
                 --fpr 0.01 \
                 --epochs 150

############################
# Cell Bender Sample 187v2 #
############################
CUDA_VISIBLE_DEVICES=3 cellbender remove-background \
                 --input /home/esk17/1.projects/colitis/data/CD_martin/GSM3972022_187/ \
                 --output /home/esk17/1.projects/colitis/data/CD_martin_cellbender/GSM3972022_187v2cb.h5 \
                 --cuda \
                 --low-count-threshold 5\
                 --expected-cells 2000 \
                 --total-droplets-included 7000 \
                 --fpr 0.01 \
                 --epochs 150

##########################
# Cell Bender Sample 189 #
##########################
CUDA_VISIBLE_DEVICES=1 cellbender remove-background \
                 --input /home/esk17/1.projects/colitis/data/CD_martin/GSM3972023_189/ \
                 --output /home/esk17/1.projects/colitis/data/CD_martin_cellbender/GSM3972023_189cb.h5 \
                 --cuda \
                 --expected-cells 1000 \
                 --total-droplets-included 7000 \
                 --fpr 0.01 \
                 --epochs 150

##########################
# Cell Bender Sample 189v2 #
##########################
CUDA_VISIBLE_DEVICES=6 cellbender remove-background \
                 --input /home/esk17/1.projects/colitis/data/CD_martin/GSM3972023_189/ \
                 --output /home/esk17/1.projects/colitis/data/CD_martin_cellbender/GSM3972023_189v2cb.h5 \
                 --cuda \
                 --low-count-threshold 10\
                 --expected-cells 4000 \
                 --total-droplets-included 30000 \
                 --fpr 0.01 \
                 --epochs 150

##########################
# Cell Bender Sample 190 #
##########################
CUDA_VISIBLE_DEVICES=2 cellbender remove-background \
                 --input /home/esk17/1.projects/colitis/data/CD_martin/GSM3972024_190/ \
                 --output /home/esk17/1.projects/colitis/data/CD_martin_cellbender/GSM3972024_190cb.h5 \
                 --cuda \
                 --expected-cells 5000 \
                 --total-droplets-included 10000 \
                 --fpr 0.01 \
                 --epochs 150

##########################
# Cell Bender Sample 192 #
##########################
CUDA_VISIBLE_DEVICES=5 cellbender remove-background \
                 --input /home/esk17/1.projects/colitis/data/CD_martin/GSM3972025_192/ \
                 --output /home/esk17/1.projects/colitis/data/CD_martin_cellbender/GSM3972025_192cb.h5 \
                 --cuda \
                 --expected-cells 5000 \
                 --total-droplets-included 15000 \
                 --fpr 0.01 \
                 --epochs 150

############################
# Cell Bender Sample 192v2 #
############################
CUDA_VISIBLE_DEVICES=5 cellbender remove-background \
                 --input /home/esk17/1.projects/colitis/data/CD_martin/GSM3972025_192/ \
                 --output /home/esk17/1.projects/colitis/data/CD_martin_cellbender/GSM3972025_192v2cb.h5 \
                 --cuda \
                 --low-count-threshold 10\
                 --expected-cells 5000 \
                 --total-droplets-included 30000 \
                 --fpr 0.01 \
                 --epochs 150

##########################
# Cell Bender Sample 193 #
##########################
CUDA_VISIBLE_DEVICES=7 cellbender remove-background \
                 --input /home/esk17/1.projects/colitis/data/CD_martin/GSM3972026_193/ \
                 --output /home/esk17/1.projects/colitis/data/CD_martin_cellbender/GSM3972026_193cb.h5 \
                 --cuda \
                 --expected-cells 5000 \
                 --total-droplets-included 15000 \
                 --fpr 0.01 \
                 --epochs 150

############################
# Cell Bender Sample 193v2 #
############################
CUDA_VISIBLE_DEVICES=5 cellbender remove-background \
                 --input /home/esk17/1.projects/colitis/data/CD_martin/GSM3972026_193/ \
                 --output /home/esk17/1.projects/colitis/data/CD_martin_cellbender/GSM3972026_193v2cb.h5 \
                 --cuda \
                 --low-count-threshold 10\
                 --expected-cells 5000 \
                 --total-droplets-included 30000 \
                 --fpr 0.01 \
                 --epochs 150

##########################
# Cell Bender Sample 195 #
##########################
CUDA_VISIBLE_DEVICES=6 cellbender remove-background \
                 --input /home/esk17/1.projects/colitis/data/CD_martin/GSM3972027_195/ \
                 --output /home/esk17/1.projects/colitis/data/CD_martin_cellbender/GSM3972027_195cb.h5 \
                 --cuda \
                 --expected-cells 4000 \
                 --total-droplets-included 20000 \
                 --fpr 0.01 \
                 --epochs 150

##########################
# Cell Bender Sample 196 #
##########################
CUDA_VISIBLE_DEVICES=4 cellbender remove-background \
                 --input /home/esk17/1.projects/colitis/data/CD_martin/GSM3972028_196/ \
                 --output /home/esk17/1.projects/colitis/data/CD_martin_cellbender/GSM3972028_196cb.h5 \
                 --cuda \
                 --expected-cells 4000 \
                 --total-droplets-included 20000 \
                 --fpr 0.01 \
                 --epochs 150

##########################
# Cell Bender Sample 208 #
##########################
CUDA_VISIBLE_DEVICES=3 cellbender remove-background \
                 --input /home/esk17/1.projects/colitis/data/CD_martin/GSM3972029_208/ \
                 --output /home/esk17/1.projects/colitis/data/CD_martin_cellbender/GSM3972029_208cb.h5 \
                 --cuda \
                 --expected-cells 800 \
                 --total-droplets-included 7000 \
                 --fpr 0.01 \
                 --epochs 150

##########################
# Cell Bender Sample 209 #
##########################
CUDA_VISIBLE_DEVICES=6 cellbender remove-background \
                 --input /home/esk17/1.projects/colitis/data/CD_martin/GSM3972030_209/ \
                 --output /home/esk17/1.projects/colitis/data/CD_martin_cellbender/GSM3972030_209cb.h5 \
                 --cuda \
                 --expected-cells 900 \
                 --total-droplets-included 10000 \
                 --fpr 0.01 \
                 --epochs 150

#### The following code is in R script ####
################################################
# 6- Recluster martin from scratch (with new filter method)
################################################
dirnames <- dir("data/CD_martin",pattern="GSM")
date <- "121320" # remember to make directory in "saved_objects" and "figures"
qc_CD_postcb(dirnames[1], date)
save_figures_CD(dirnames[1], date)

for(i in 2:length(input.dirs)) {
	qc_CD(input.dirs[i], date)
}
for(i in 2:length(input.dirs)) {
	save_figures_CD(input.dirs[i], date)
}

y <- sapply(input.dirs,function(x) get_info_CD(x,date))
y.print <- as.data.frame(t(y))

library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")
xtable(y.print, display=c("s","d","e","d","d","f","d","d","d","d","d"))





