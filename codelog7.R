######################################################################
## 
## Analysis of all T cells with PD1
## 
## edward_kim@college.harvard.edu - August 2020
######################################################################

##############################
# 0 - Load librairies
##############################
library(dplyr)
library(ggplot2)
library(Seurat)
library(scran)
library(pheatmap)

##############################
# 1 - Definitions & Settings
##############################
'%ni%' = Negate('%in%')
all.genes <- rownames(s2_postQC) # do after importing counts matrix
trvgenes <- all.genes[grepl(x=all.genes, pattern = "^TRAV|^TRBV|^TRGV|^TRDV")]
igvgenes <- all.genes[grepl(x=all.genes, pattern = "^IG.V")]
options(future.globals.maxSize = 11000 * 1024^2)
setwd("~/1.projects/PD1")

############################## 
# 2 - Source file 
##############################
source("scripts/Stat111functions.R")
#source("scripts/new_clustering_functions.R")

#####################
## Basic Inference ##
#####################
# Metadata Key: 
# sampletype = patient ID
# colitis = Colitis | Control | No-Colitis
# colitis2 = Colitis(aPD1) | Colitis (combo) | Control | 
# 			 No-Colitis (aPD1) | No-colitis (combo)
# condition = Batch1 | Batch2 | Batch3 | Batch4 | BatchPD1
# orig.ident = CD3_Tcell
# others: Phase, S.Score, G2M.Score
# will create: colitis2
# 
# PD-1 info:
# colitis+: "p101519_CD3","p020620_CD45"
# no colitis: "DFCI-1294_CD3","DFCI-1545_CD3","DFCI-1618_CD3"
s <- readRDS("saved_objects/onlyTcells.rds")

# Viz 1: UMAP - group by colitis
p.status <- DimPlot(s, reduction = "umap", label=F, group.by="colitis")
ggsave("figures/onlyTcells_group_colitis.pdf")
# Viz 2: UMAP - split by colitis
p.status.split <- DimPlot(s, reduction = "umap", label=F, split.by="colitis", ncol=3)
ggsave("figures/onlyTcells_split_colitis.pdf", width= 12, height= 4, units= "in")
# Viz 3: UMAP - group by colitis2
p.colitis2 <- DimPlot(s, reduction = "umap", label=F, group.by="colitis2")
ggsave("figures/onlyTcells_group_colitis2.pdf")
# Viz 4: cluster proportions
p.prop <-


DimPlot(seurat.object, reduction = "umap", label=F, split.by="status")
ggsave("figures/onlyTcells_split_status.pdf")

