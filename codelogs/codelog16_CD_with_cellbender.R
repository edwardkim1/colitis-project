############################################################
## Computational Integration Analysis
## edward_kim@college.harvard.edu - December 2020
#############################################################

##############################
# 0 - Load librairies
##############################
library(dplyr)
library(stringr)
library(ggplot2)
library(Seurat)
library(future)
library(future.apply)
library(hdf5r)
library(parallel)
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
plot_brp(s[[1]]$nCount_RNA, 1900, 10000, dirnames[1]) # 69 run v2
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
# 6- Preprocess (with new filter method)
################################################
dirnames <- dir("data/CD_martin",pattern="GSM")
dirnames[c(1,8,11,13,14,15,17,18)] <- str_glue("{dirnames[c(1,8,11,13,14,15,17,18)]}v2")
date <- "121320" # remember to make directory in "saved_objects" and "figures"

qc_CD_postcb(dirnames[1], date)
save_figures_CD_postcb(dirnames[1], date)


for(i in 2:length(dirnames)) {
	qc_CD_postcb(dirnames[i], date)
}
for(i in 2:length(dirnames)) {
	save_figures_CD_postcb(dirnames[i], date)
}



y <- sapply(dirnames,function(x) get_info_CD_postcb(x,date))
y.print <- as.data.frame(t(y))

library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")
xtable(y.print, display=c("s","d","d","f","d","d","d","d","d"))

################################################
# 6- cleaning clusters
################################################
dirnames <- dir("data/CD_martin",pattern="GSM")
dirnames[c(1,8,11,13,14,15,17,18)] <- str_glue("{dirnames[c(1,8,11,13,14,15,17,18)]}v2")
date <- "121320" # remember to make directory in "saved_objects" and "figures"
cells <- mclapply(paste("saved_objects/CD_martin_qc_", date, "/", dirnames, "_3.RDS", sep=""),readRDS)
#### Removing bad clusters from each sample
cells[[1]] <- subset(cells[[1]], idents= c("5","6"), invert=TRUE) # 69 run v2
cells[[2]] <- subset(cells[[2]], idents= c("10","12","13"), invert=TRUE)# 68 run
cells[[3]] <- subset(cells[[3]], idents= c("5","7"), invert=TRUE)# 122 run
cells[[4]] <- subset(cells[[4]], idents= c("3"), invert=TRUE) # 123 run
cells[[5]] <- subset(cells[[5]], idents= c("7","11"), invert=TRUE) # 128 run
cells[[6]] <- subset(cells[[6]], idents= c("8","9"), invert=TRUE) # 129 run
cells[[7]] <- subset(cells[[7]], idents= c("8","9"), invert=TRUE) # 135 run
cells[[8]] <- subset(cells[[8]], idents= c("9"), invert=TRUE) # 138 run v2
#cells[[9]] <- subset(cells[[9]], idents= c(...), invert=TRUE) # 158 run
cells[[10]] <- subset(cells[[10]], idents= c("4"), invert=TRUE) # 159 run
cells[[11]] <- subset(cells[[11]], idents= c("3"), invert=TRUE) # 180 run v2
#cells[[12]] <- subset(cells[[12]], idents= c(...), invert=TRUE) # 181 run
#cells[[13]] <- subset(cells[[13]], idents= c(...), invert=TRUE) # 186 run v2
cells[[14]] <- subset(cells[[14]], idents= c("9"), invert=TRUE) # 187 run v2 
cells[[15]] <- subset(cells[[15]], idents= c("12","13","14"), invert=TRUE) # 189 run v2
cells[[16]] <- subset(cells[[16]], idents= c("12"), invert=TRUE) # 190 run
cells[[17]] <- subset(cells[[17]], idents= c("10","12","13"), invert=TRUE) # 192 run v2
cells[[18]] <- subset(cells[[18]], idents= c("13"), invert=TRUE) # 193 run v2
cells[[19]] <- subset(cells[[19]], idents= c("12"), invert=TRUE) # 195 run
cells[[20]] <- subset(cells[[20]], idents= c("7","8","9","10","11","12"), invert=TRUE) # 196 run
cells[[21]] <- subset(cells[[21]], idents= c("12","13","14"), invert=TRUE) # 208 run
cells[[22]] <- subset(cells[[22]], idents= c("5","6","10","11","12","13"), invert=TRUE) # 209 run

####################################
# Naive integration with cells from all patients
####################################
s <- merge(cells[[1]], c(cells[[2]],cells[[3]],cells[[4]], cells[[5]], cells[[6]], cells[[7]], cells[[8]], cells[[9]], cells[[10]], cells[[11]], cells[[12]], cells[[13]], cells[[14]], cells[[15]], cells[[16]], cells[[17]], cells[[18]], cells[[19]], cells[[20]], cells[[21]],cells[[22]]))

########################################################
# 3 - Add patient annotations and inflammed annotations
########################################################
library(stringr)
# rename sampletype to sample.id and merge with orig.ident material
s$sample.id <- s$orig.ident
s$orig.ident.backup <- s$orig.ident
#### mutate orig.ident and edit sample.id
s$orig.ident <- ifelse(grepl(s$orig.ident,pattern="^GSM*"), "CD_martin", s$orig.ident)
s$sample.id <- ifelse(str_detect(s$sample.id,"v2"), str_sub(s$sample.id,1,-3), s$sample.id) %>%
        word(2, sep="_") %>% 
        paste("CD",.,sep="_")
#### add patient.id
s <- add_patient_id_all(s)
s$patient.id <- str_replace(s$patient.id,"Patient ","P")
s$patient.id <- factor(s$patient.id, levels = c("P5","P6","P7","P8","P10","P11","P12","P13","P14","P15","P16"))
# add inflamed id for Crohns
inflamed <-c("CD_69$","CD_122$","CD_128$","CD_138$","CD_158$","CD_181$","CD_187$","CD_193$","CD_190$","CD_196$","CD_209$")
toMatch <- paste(inflamed,collapse="|")
status <- ifelse(grepl(s$sample.id, pattern=toMatch),"CD_inflamed","CD_uninflamed")
s$colitis <- ifelse(s$orig.ident == "CD_martin", status ,s$colitis)
s$colitis2 <- ifelse(s$orig.ident == "CD_martin", status ,s$colitis2)
## Convert colitis2 and sample.id to factors
s$colitis2 <- factor(s$colitis2, levels = c("CD_inflamed","CD_uninflamed"))
s$sample.id <- factor(s$sample.id)
# removing unapplicable meta.data columns
m.names <- colnames(s@meta.data)
columns.to.remove <- m.names[grepl("^pANN*|^DF.c*|seurat_clusters|^RNA*|cells.to.remove|Single*|percent.mito",m.names)]
for(i in columns.to.remove) {
  s[[i]] <- NULL
}
# remove genes with less than 10
x <- s[["RNA"]]@counts
nnz_row <- tabulate(x@i + 1)
keep <- rownames(x)[nnz_row>10]
s <- subset(s, features = keep)

#############################
# 4 - perform clustering
#############################
s <- PercentageFeatureSet(s, pattern = "^MT-", col.name = "percent.mt")

s <- NormalizeData(s, normalization.method = "LogNormalize", scale.factor = 10000) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000) ## Highly variable genes
hvf <- VariableFeatures(s)
'%ni%' = Negate('%in%')
all.genes <- rownames(s)
trvgenes <- all.genes[grepl(x=all.genes, pattern = "^TRAV|^TRBV|^TRGV|^TRDV")]
igvgenes <- all.genes[grepl(x=all.genes, pattern = "^IG.V")]
VariableFeatures(s) <- hvf[hvf %ni% c(trvgenes,igvgenes)]
# Do the rest
s <- ScaleData(s, features = all.genes, vars.to.regress=c("percent.mt")) %>% RunPCA(features=VariableFeatures(s))
s <-  cluster_umap(s, npc=20, k.param=120, resolution=0.5, min.dist=0.3)
DimPlot(s, label=T)
ggsave("figures/CD_martin_121320_integration/martin_naive_npc20_k120_umap.pdf")
DimPlot(s, group.by="sample.id" , label=F)
ggsave("figures/CD_martin_121320_integration/martin_naive_npc20_k120_persample.pdf")
# distribution by sample
data.frame(cluster = s$seurat_clusters, patient= s$sample.id) %>%
        ggplot() + geom_bar(
                mapping = aes(x= patient, fill= cluster),
                position = "fill"
        ) + labs(y="proportion")+ theme_classic() + coord_flip()

# make directory: CD_martin_121320_integration
ggsave("figures/CD_martin_121320_integration/martin_naive_npc20_k120_distpersample.pdf")
## Replicating supplemental figure for QC validation
library(ggplot2)
df <- as.data.frame(table(s$patient.id,s$colitis))
colnames(df) <- c("patient","sampletype", "no.of.cells")
ggplot(df, aes(patient, no.of.cells, fill = sampletype)) + 
  geom_bar(stat="identity",position="dodge") + 
  scale_fill_brewer(palette = "Set1") +
  theme_classic()
ggsave("figures/CD_martin_121320_integration/martin_naive_npc20_k120_counts_per_patient.pdf")

####################################
# 5 - Cluster Annotations
####################################
load('data/ssuo_CD3/seurat.object.RData')
Tcell <- process_ssuo_Tcells(Tcell)
martin <- read.csv(file="data/martin_categories.csv", header=TRUE,sep=",") %>% process_martin()

require(SingleR)
require(BiocParallel)
pred.CPI <- SingleR(method='cluster',
        test=as.matrix(s[['RNA']]@data), 
        clusters= s$seurat_clusters,
        ref=as.matrix(Tcell[['RNA']]@data), 
        labels=Tcell$names,
        de.method="wilcox",
        BPPARAM=MulticoreParam(10))
s$SingleR.Luoma <- s$seurat_clusters
levels(s$SingleR.Luoma) <- pred.CPI$labels

pred.Martin <- SingleR(method='cluster',  
        test=s[['RNA']]@data, 
        clusters= s$seurat_clusters,
        ref=martin, 
        labels=colnames(martin),
        BPPARAM=MulticoreParam(10))
s$SingleR.Martin <- s$seurat_clusters
levels(s$SingleR.Martin) <- pred.Martin$labels

ref <- MonacoImmuneData()
pred.Monaco <- SingleR(method="cluster",
        test=s[['RNA']]@data,
        clusters= s$seurat_clusters, 
        ref=ref, 
        labels=ref$label.fine, 
        BPPARAM=MulticoreParam(10))
s$SingleR.Monaco <- s$seurat_clusters
levels(s$SingleR.Monaco) <- pred.Monaco$labels


require(ggpubr)
p1 <- DimPlot(s, reduction = "umap", label=T, group.by="SingleR.Martin") + labs(title="Martin reference")
p2 <- DimPlot(s, reduction = "umap", label=T, group.by="SingleR.Luoma") + labs(title="Luoma-Suo reference")
p3 <- DimPlot(s, reduction = "umap", label=T, group.by="SingleR.Monaco") + labs(title="Monaco reference")
ggarrange(p1,p2,p3, ncol=3, nrow=1)
ggsave("figures/CD_martin_121320_integration/martin_naive_npc20_k120_singleR.pdf", width= 24, height= 6, units= "in")
######################################################
# 3 - CD3E and CD4/CD8 and other QC metrics per cluster
######################################################
# violin plot of CD3E per cluster
p <- VlnPlot(s,features = "CD3E", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 90),legend.position = "none", axis.title.x=element_blank())
ggsave("figures/CD_martin_121320_integration/martin_naive20120_VlnPlot_CD3E_by_cluster.pdf", width= 7, height= 7, units= "in")
# feature plots the CD4 and CD8 split
FeaturePlot(s, features= c('CD4','CD8A','CD8B','CD3E'), reduction="umap")
ggsave("figures/CD_martin_121320_integration/martin_naive20120_CD4xCD8_featureplot.pdf", width= 12, height= 12, units= "in")
# violin plot of %mito per cluster
p <- VlnPlot(s,features = "percent.mt", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 90),legend.position = "none", axis.title.x=element_blank())
ggsave("figures/CD_martin_121320_integration/martin_naive20120_VlnPlot_percentmito_by_cluster.pdf", width= 7, height= 7, units= "in")
# violin plot of nFeature per cluster
p <- VlnPlot(s,features = "nFeature_RNA", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 90),legend.position = "none", axis.title.x=element_blank())
ggsave("figures/CD_martin_121320_integration/martin_naive20120_VlnPlot_nFeatureRNA_by_cluster.pdf", width= 7, height= 7, units= "in")
# feature plots the IGKC and JCHAIN
FeaturePlot(s, features= c('IGKC','JCHAIN'), reduction="umap")
ggsave("figures/CD_martin_121320_integration/martin_naive20120_IGKC-JCHAIN_featureplot.pdf", width= 12, height= 7, units= "in")

# violin plot of IGKC per cluster
p <- VlnPlot(s,features = "IGKC", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 90),legend.position = "none", axis.title.x=element_blank())
ggsave(paste("figures/CD_martin_121320_integration/martin_naive20120_VlnPlot_IGKC_by_cluster.pdf", sep=""), width= 7, height= 7, units= "in")

# save RDS
saveRDS(s, paste("saved_objects/CD_martin_qc_", date, "/martin_20120.rds", sep=""))

