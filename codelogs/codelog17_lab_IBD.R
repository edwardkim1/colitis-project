############################################################
## Mapping lab IBD on the CPI dataset
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

####################################
# cluster IBD (negative)
###################################
# mkdir CD_luoma_011621
date <- "011621"
dirname <- "p1089neg"
input.directory <- "data/CD_luoma/p1089neg-GEX-Pool7/outs/filtered_feature_bc_matrix"
qc_CD(dirname, date, "luoma", input.directory)
save_figures_CD(dirname, date, "luoma", input.directory)





