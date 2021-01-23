############################################################
## Integration Analysis
## edward_kim@college.harvard.edu - October 2020
#############################################################

##############################
# 0 - Load librairies
##############################
library(dplyr)
library(ggplot2)
library(Seurat)
library(future)
library(future.apply)
library(fgsea)
library(harmony)
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
# 2 - Remove cluster 8; negative filter of B cells
##############################
s.a <- readRDS("saved_objects/cdcpi3_anchor.rds")
s.a$



############################################
# 3 - DA CPI colitis vs CPI (no colitis and control)
############################################
s.h <- readRDS("saved_objects/cdcpi3_harmony.rds")
s.a <- readRDS("saved_objects/cdcpi3_anchor.rds")

# Anchor:
library(edgeR)
x <- subset(s.a, orig.ident == "CD3_Tcell") # extract CPI dataset
x$CPI_colitis <- x$colitis == "CPI_colitis"
design <- function(y) model.matrix(~factor(CPI_colitis),y)
out <- DA_analysis(x, design, title = "cdcpi3_anchor_CPI_colitis_vs_otherCPI")
summary(out$y.ab$common.dispersion)
summary(out$fit.ab$var.prior)
summary(out$fit.ab$df.prior)
summary(decideTests(out$res))
topTags(out$res)

# Harmony:
library(edgeR)
x <- subset(s.h, orig.ident == "CD3_Tcell") # extract CPI dataset
x$CPI_colitis <- x$colitis == "CPI_colitis"
design <- function(y) model.matrix(~factor(CPI_colitis),y)
out <- DA_analysis(x, design, title = "cdcpi3h_CPI_colitis_vs_otherCPI")
summary(out$y.ab$common.dispersion)
summary(out$fit.ab$var.prior)
summary(out$fit.ab$df.prior)
summary(decideTests(out$res))
topTags(out$res)

############################################
# 3 - DA with interactive terms (WIP)
############################################
require(edgeR)
x <- s.a
abundances <- table(x$seurat_clusters, x$sample.id) %>% unclass() %>% t() %>% as.data.frame()


# output
return(list(y.ab = y.ab, fit.ab = fit.ab,res = res))

# Anchor:
library(edgeR)
x <- subset(s.a, orig.ident == "CD3_Tcell") # extract CPI dataset
x$CPI_colitis <- x$colitis == "CPI_colitis"
design <- function(y) model.matrix(~factor(CPI_colitis),y)
out <- DA_analysis(x, design, title = "cdcpi3_anchor_CPI_colitis_vs_otherCPI")
summary(out$y.ab$common.dispersion)
summary(out$fit.ab$var.prior)
summary(out$fit.ab$df.prior)
summary(decideTests(out$res))
topTags(out$res)


############################################
# 4 - Other Figures
############################################
# proportion of each dataset by cluster (harmony)
data.frame(dataset = s.h$orig.ident, cluster= s.h$seurat_clusters) %>%
	ggplot() + geom_bar(
		mapping = aes(x= cluster, fill= dataset),
		position = "fill"
	) + labs(y="proportion")+ theme_classic()
ggsave("figures/cdcpi3_harmony/cdcpi3h_cluster_dataset_per_cluster.pdf",width= 7, height= 7, units= "in")
# proportion of each dataset by cluster (anchor)
data.frame(dataset = s.a$orig.ident, cluster= s.a$seurat_clusters) %>%
	ggplot() + geom_bar(
		mapping = aes(x= cluster, fill= dataset),
		position = "fill"
	) + labs(y="proportion")+ theme_classic()
ggsave("figures/cdcpi3/cdcpi3_anchor_dataset_proportions.pdf",width= 7, height= 7, units= "in")

# differential gene analysis by cluster split by dataset
x <- subset(s.a, orig.ident == "CD3_Tcell") # extract CPI dataset
DE_heatmap(x,"figures/cdcpi3_compare/cdcpi3_CPI_DE_avgexp.pdf")
x <- subset(s.a, orig.ident == "CD_martin") # extract CD dataset
DE_heatmap(x,"figures/cdcpi3_compare/cdcpi3_CD_DE_avgexp.pdf")
x <- subset(s.h, orig.ident == "CD3_Tcell") # extract CPI dataset
DE_heatmap(x,"figures/cdcpi3_compare/cdcpi3h_CPI_DE_avgexp.pdf")
x <- subset(s.h, orig.ident == "CD_martin") # extract CD dataset
DE_heatmap(x,"figures/cdcpi3_compare/cdcpi3h_CD_DE_avgexp.pdf")

# differential gene analysis by colitis condtion in Treg cluster: anchor
FCcutoff <- log(1.5, base=2)

x <- subset(s.a, orig.ident == "CD3_Tcell") # extract CPI dataset
s.markers <- FindMarkers(x, ident.1 ="CPI_colitis" ,group.by= "colitis", subset.ident="6", min.pct = 0.25, logfc.threshold = 0.01, verbose=F)
library(EnhancedVolcano)
head(s.markers, n=20) %>% rownames() -> labels
up.labels <- rownames(s.markers[s.markers$avg_logFC > FCcutoff,])
down.labels <- rownames(s.markers[s.markers$avg_logFC < -FCcutoff,])
labels <- unique(c(labels, up.labels,down.labels))
EnhancedVolcano(s.markers,
    lab = rownames(s.markers),
    x = 'avg_logFC',
    y = 'p_val_adj',
  #  xlim = c(-20, 20),
    title = 'Treg comparison (CPI_colitis vs. CPI_controls)',
    pCutoff = 0.05,
    FCcutoff = FCcutoff,
    pointSize = 3.0,
    labSize = 3.0,
   selectLab = labels,
    drawConnectors = T,
    widthConnectors = 0.2,
    colConnectors = 'grey30')
ggsave("figures/cdcpi3_compare/cdcpi3_DE_avgexp_treg_CPI_colitis.pdf", width= 10, height= 8, units= "in")

x <- subset(s.a, orig.ident == "CD_martin") # extract CD dataset
s.markers <- FindMarkers(x, ident.1 ="CD_inflamed" ,group.by= "colitis", subset.ident="6", min.pct = 0.25, logfc.threshold = 0.01, verbose=F)
head(s.markers, n=20) %>% rownames() -> labels
up.labels <- rownames(s.markers[s.markers$avg_logFC > FCcutoff,])
down.labels <- rownames(s.markers[s.markers$avg_logFC < -FCcutoff,])
labels <- unique(c(labels, up.labels,down.labels))
EnhancedVolcano(s.markers,
    lab = rownames(s.markers),
    x = 'avg_logFC',
    y = 'p_val_adj',
  #  xlim = c(-20, 20),
    title = 'Treg comparison (CD_inflamed vs. CD_uninflamed)',
    pCutoff = 0.05,
    FCcutoff = FCcutoff,
    pointSize = 3.0,
    labSize = 3.0,
   selectLab = labels,
    drawConnectors = T,
    widthConnectors = 0.2,
    colConnectors = 'grey30')
ggsave("figures/cdcpi3_compare/cdcpi3_DE_avgexp_treg_CD_inflamed.pdf", width= 10, height= 8, units= "in")

# differential gene analysis by colitis condtion in Treg cluster: harmony
x <- subset(s.h, orig.ident == "CD3_Tcell") # extract CPI dataset
s.markers <- FindMarkers(x, ident.1 ="CPI_colitis" ,group.by= "colitis", subset.ident="6", min.pct = 0.25, logfc.threshold = 0.01, verbose=F)
library(EnhancedVolcano)
head(s.markers, n=20) %>% rownames() -> labels
EnhancedVolcano(s.markers,
    lab = rownames(s.markers),
    x = 'avg_logFC',
    y = 'p_val_adj',
  #  xlim = c(-20, 20),
    title = 'Treg comparison (CPI_colitis vs. CPI_controls)',
    pCutoff = 0.05,
    FCcutoff = 0.8,
    pointSize = 3.0,
    labSize = 3.0,
   selectLab = labels,
    drawConnectors = F,
    widthConnectors = 0.2,
    colConnectors = 'grey30')
ggsave("figures/cdcpi3_compare/cdcpi3h_DE_avgexp_treg_CPI_colitis.pdf", width= 10, height= 8, units= "in")

x <- subset(s.h, orig.ident == "CD_martin") # extract CD dataset
s.markers <- FindMarkers(x, ident.1 ="CD_inflamed" ,group.by= "colitis", subset.ident="6", min.pct = 0.25, logfc.threshold = 0.01, verbose=F)
head(s.markers, n=20) %>% rownames() -> labels
EnhancedVolcano(s.markers,
    lab = rownames(s.markers),
    x = 'avg_logFC',
    y = 'p_val_adj',
  #  xlim = c(-20, 20),
    title = 'Treg comparison (CD_inflamed vs. CD_uninflamed)',
    pCutoff = 0.05,
    FCcutoff = 0.8,
    pointSize = 3.0,
    labSize = 3.0,
   selectLab = labels,
    drawConnectors = F,
    widthConnectors = 0.2,
    colConnectors = 'grey30')
ggsave("figures/cdcpi3_compare/cdcpi3h_DE_avgexp_treg_CD_inflamed.pdf", width= 10, height= 8, units= "in")


# differential gene analysis by colitis condtion in CD4 Trm cluster: anchor
x <- subset(s.a, orig.ident == "CD3_Tcell") # extract CPI dataset
s.markers <- FindMarkers(x, ident.1 ="CPI_colitis" ,group.by= "colitis", subset.ident= "0", min.pct = 0.25, logfc.threshold = 0.01, verbose=F)
head(s.markers, n=20) %>% rownames() -> labels
up.labels <- rownames(s.markers[s.markers$avg_logFC > FCcutoff,])
down.labels <- rownames(s.markers[s.markers$avg_logFC < -FCcutoff,])
labels <- unique(c(labels, up.labels,down.labels))
EnhancedVolcano(s.markers,
    lab = rownames(s.markers),
    x = 'avg_logFC',
    y = 'p_val_adj',
  #  xlim = c(-20, 20),
    title = 'CD4 Trm comparison (CPI_colitis vs. CPI_controls)',
    pCutoff = 0.05,
    FCcutoff = FCcutoff,
    pointSize = 3.0,
    labSize = 3.0,
   selectLab = labels,
    drawConnectors = T,
    widthConnectors = 0.2,
    colConnectors = 'grey30')
ggsave("figures/cdcpi3_compare/cdcpi3_DE_CD4Trm_CPI_colitis.pdf", width= 10, height= 8, units= "in")

x1 <- subset(s.a, orig.ident == "CD_martin") # extract CD dataset
s.markers <- FindMarkers(x1, ident.1 ="CD_inflamed" ,group.by= "colitis", subset.ident="0", min.pct = 0.25, logfc.threshold = 0.01, verbose=F)
head(s.markers, n=20) %>% rownames() -> labels
up.labels <- rownames(s.markers[s.markers$avg_logFC > FCcutoff,])
down.labels <- rownames(s.markers[s.markers$avg_logFC < -FCcutoff,])
labels <- unique(c(labels, up.labels,down.labels))
EnhancedVolcano(s.markers,
    lab = rownames(s.markers),
    x = 'avg_logFC',
    y = 'p_val_adj',
  #  xlim = c(-20, 20),
    title = 'CD4 Trm comparison (CD_inflamed vs. CD_uninflamed)',
    pCutoff = 0.05,
    FCcutoff = FCcutoff,
    pointSize = 3.0,
    labSize = 3.0,
   selectLab = labels,
    drawConnectors = F,
    widthConnectors = 0.2,
    colConnectors = 'grey30')
ggsave("figures/cdcpi3_compare/cdcpi3_DE_CD4Trm_CD_inflamed.pdf", width= 10, height= 8, units= "in")

# differential gene analysis by colitis condtion in CD8 Trm cluster: anchor
s.markers <- FindMarkers(x, ident.1 ="CPI_colitis" ,group.by= "colitis", subset.ident= "1", min.pct = 0.25, logfc.threshold = 0.01, verbose=F)
head(s.markers, n=20) %>% rownames() -> labels
up.labels <- rownames(s.markers[s.markers$avg_logFC > FCcutoff,])
down.labels <- rownames(s.markers[s.markers$avg_logFC < -FCcutoff,])
labels <- unique(c(labels, up.labels,down.labels))
EnhancedVolcano(s.markers,
    lab = rownames(s.markers),
    x = 'avg_logFC',
    y = 'p_val_adj',
  #  xlim = c(-20, 20),
    title = 'CD8 Trm comparison (CPI_colitis vs. CPI_controls)',
    pCutoff = 0.05,
    FCcutoff = FCcutoff,
    pointSize = 3.0,
    labSize = 3.0,
   selectLab = labels,
    drawConnectors = T,
    widthConnectors = 0.2,
    colConnectors = 'grey30')
ggsave("figures/cdcpi3_compare/cdcpi3_DE_CD8Trm_CPI_colitis.pdf", width= 10, height= 8, units= "in")


s.markers <- FindMarkers(x1, ident.1 ="CD_inflamed" ,group.by= "colitis", subset.ident="1", min.pct = 0.25, logfc.threshold = 0.01, verbose=F)
head(s.markers, n=20) %>% rownames() -> labels
up.labels <- rownames(s.markers[s.markers$avg_logFC > FCcutoff,])
down.labels <- rownames(s.markers[s.markers$avg_logFC < -FCcutoff,])
labels <- unique(c(labels, up.labels,down.labels))
EnhancedVolcano(s.markers,
    lab = rownames(s.markers),
    x = 'avg_logFC',
    y = 'p_val_adj',
  #  xlim = c(-20, 20),
    title = 'CD8 Trm comparison (CD_inflamed vs. CD_uninflamed)',
    pCutoff = 0.05,
    FCcutoff = FCcutoff,
    pointSize = 3.0,
    labSize = 3.0,
   selectLab = labels,
    drawConnectors = F,
    widthConnectors = 0.2,
    colConnectors = 'grey30')
ggsave("figures/cdcpi3_compare/cdcpi3_DE_CD8Trm_CD_inflamed.pdf", width= 10, height= 8, units= "in")

# differential gene analysis by colitis condtion in Th1 effector cluster: anchor
s.markers <- FindMarkers(x, ident.1 ="CPI_colitis" ,group.by= "colitis", subset.ident= "7", min.pct = 0.25, logfc.threshold = 0.01, verbose=F)
FCcutoff <- log(1.5, base=2)
head(s.markers, n=20) %>% rownames() -> labels
up.labels <- rownames(s.markers[s.markers$avg_logFC > FCcutoff,])
down.labels <- rownames(s.markers[s.markers$avg_logFC < -FCcutoff,])
labels <- unique(c(labels, up.labels,down.labels))
EnhancedVolcano(s.markers,
    lab = rownames(s.markers),
    x = 'avg_logFC',
    y = 'p_val_adj',
  #  xlim = c(-20, 20),
    title = 'Th1 effector comparison (CPI_colitis vs. CPI_controls)',
    pCutoff = 0.05,
    FCcutoff = FCcutoff,
    pointSize = 3.0,
    labSize = 3.0,
   selectLab = labels,
    drawConnectors = T,
    widthConnectors = 0.2,
    colConnectors = 'grey30')
ggsave("figures/cdcpi3_compare/cdcpi3_DE_Th1effector_CPI_colitis.pdf", width= 10, height= 8, units= "in")


s.markers <- FindMarkers(x1, ident.1 ="CD_inflamed" ,group.by= "colitis", subset.ident="7", min.pct = 0.25, logfc.threshold = 0.01, verbose=F)
FCcutoff <- log(1.5, base=2)
head(s.markers, n=20) %>% rownames() -> labels
up.labels <- rownames(s.markers[s.markers$avg_logFC > FCcutoff,])
down.labels <- rownames(s.markers[s.markers$avg_logFC < -FCcutoff,])
labels <- unique(c(labels, up.labels,down.labels))
EnhancedVolcano(s.markers,
    lab = rownames(s.markers),
    x = 'avg_logFC',
    y = 'p_val_adj',
  #  xlim = c(-20, 20),
    title = 'Th1 effector comparison (CD_inflamed vs. CD_uninflamed)',
    pCutoff = 0.05,
    FCcutoff = FCcutoff,
    pointSize = 3.0,
    labSize = 3.0,
   selectLab = labels,
    drawConnectors = F,
    widthConnectors = 0.2,
    colConnectors = 'grey30')
ggsave("figures/cdcpi3_compare/cdcpi3_DE_Th1effector_CD_inflamed.pdf", width= 10, height= 8, units= "in")

# differential gene analysis by colitis condtion in cytotoxic cluster: anchor
s.markers <- FindMarkers(x, ident.1 ="CPI_colitis" ,group.by= "colitis", subset.ident= c("3","4"), min.pct = 0.25, logfc.threshold = 0.01, verbose=F)
head(s.markers, n=20) %>% rownames() -> labels
up.labels <- rownames(s.markers[s.markers$avg_logFC > FCcutoff,])
down.labels <- rownames(s.markers[s.markers$avg_logFC < -FCcutoff,])
labels <- unique(c(labels, up.labels,down.labels))
EnhancedVolcano(s.markers,
    lab = rownames(s.markers),
    x = 'avg_logFC',
    y = 'p_val_adj',
  #  xlim = c(-20, 20),
    title = 'Cytotoxic (3+4) (CPI_colitis vs. CPI_controls)',
    pCutoff = 0.05,
    FCcutoff = FCcutoff,
    pointSize = 3.0,
    labSize = 3.0,
   selectLab = labels,
    drawConnectors = T,
    widthConnectors = 0.2,
    colConnectors = 'grey30')
ggsave("figures/cdcpi3_compare/cdcpi3_DE_cytotoxic34_CPI_colitis.pdf", width= 10, height= 8, units= "in")


s.markers <- FindMarkers(x1, ident.1 ="CD_inflamed" ,group.by= "colitis", subset.ident=c("3","4"), min.pct = 0.25, logfc.threshold = 0.05, verbose=F)
head(s.markers, n=20) %>% rownames() -> labels
up.labels <- rownames(s.markers[s.markers$avg_logFC > FCcutoff,])
down.labels <- rownames(s.markers[s.markers$avg_logFC < -FCcutoff,])
labels <- unique(c(labels, up.labels,down.labels))
EnhancedVolcano(s.markers,
    lab = rownames(s.markers),
    x = 'avg_logFC',
    y = 'p_val_adj',
  #  xlim = c(-20, 20),
    title = 'Cytotoxic (3+4) (CD_inflamed vs. CD_uninflamed)',
    pCutoff = 0.05,
    FCcutoff = FCcutoff,
    pointSize = 3.0,
    labSize = 3.0,
   selectLab = labels,
    drawConnectors = T,
    widthConnectors = 0.2,
    colConnectors = 'grey30')
ggsave("figures/cdcpi3_compare/cdcpi3_DE_cytotoxic34_CD_inflamed.pdf", width= 10, height= 8, units= "in")

# differential gene analysis by colitis condtion in naive cluster: anchor
s.markers <- FindMarkers(x, ident.1 ="CPI_colitis" ,group.by= "colitis", subset.ident= "2", min.pct = 0.25, logfc.threshold = 0.01, verbose=F)
head(s.markers, n=20) %>% rownames() -> labels
up.labels <- rownames(s.markers[s.markers$avg_logFC > FCcutoff,])
down.labels <- rownames(s.markers[s.markers$avg_logFC < -FCcutoff,])
labels <- unique(c(labels, up.labels,down.labels))
EnhancedVolcano(s.markers,
    lab = rownames(s.markers),
    x = 'avg_logFC',
    y = 'p_val_adj',
  #  xlim = c(-20, 20),
    title = 'Naive (CCR7+) CPI_colitis vs. CPI_controls',
    pCutoff = 0.05,
    FCcutoff = FCcutoff,
    pointSize = 3.0,
    labSize = 3.0,
   selectLab = labels,
    drawConnectors = T,
    widthConnectors = 0.2,
    colConnectors = 'grey30')
ggsave("figures/cdcpi3_compare/cdcpi3_DE_naive_CPI_colitis.pdf", width= 10, height= 8, units= "in")


s.markers <- FindMarkers(x1, ident.1 ="CD_inflamed" ,group.by= "colitis", subset.ident="2", min.pct = 0.25, logfc.threshold = 0.05, verbose=F)
head(s.markers, n=20) %>% rownames() -> labels
up.labels <- rownames(s.markers[s.markers$avg_logFC > FCcutoff,])
down.labels <- rownames(s.markers[s.markers$avg_logFC < -FCcutoff,])
labels <- unique(c(labels, up.labels,down.labels))
EnhancedVolcano(s.markers,
    lab = rownames(s.markers),
    x = 'avg_logFC',
    y = 'p_val_adj',
  #  xlim = c(-20, 20),
    title = 'Naive (CCR7+) CD_inflamed vs. CD_uninflamed',
    pCutoff = 0.05,
    FCcutoff = FCcutoff,
    pointSize = 3.0,
    labSize = 3.0,
   selectLab = labels,
    drawConnectors = T,
    widthConnectors = 0.2,
    colConnectors = 'grey30')
ggsave("figures/cdcpi3_compare/cdcpi3_DE_naive_CD_inflamed.pdf", width= 10, height= 8, units= "in")

# differential gene analysis by colitis condtion in MAIT cluster: anchor
s.markers <- FindMarkers(x, ident.1 ="CPI_colitis" ,group.by= "colitis", subset.ident= "11", min.pct = 0.25, logfc.threshold = 0.01, verbose=F)
head(s.markers, n=20) %>% rownames() -> labels
up.labels <- rownames(s.markers[s.markers$avg_logFC > FCcutoff,])
down.labels <- rownames(s.markers[s.markers$avg_logFC < -FCcutoff,])
labels <- unique(c(labels, up.labels,down.labels))
EnhancedVolcano(s.markers,
    lab = rownames(s.markers),
    x = 'avg_logFC',
    y = 'p_val_adj',
  #  xlim = c(-20, 20),
    title = 'MAIT (CPI_colitis vs. CPI_controls)',
    pCutoff = 0.05,
    FCcutoff = FCcutoff,
    pointSize = 3.0,
    labSize = 3.0,
   selectLab = labels,
    drawConnectors = T,
    widthConnectors = 0.2,
    colConnectors = 'grey30')
ggsave("figures/cdcpi3_compare/cdcpi3_DE_MAIT_CPI_colitis.pdf", width= 10, height= 8, units= "in")


s.markers <- FindMarkers(x1, ident.1 ="CD_inflamed" ,group.by= "colitis", subset.ident="11", min.pct = 0.25, logfc.threshold = 0.05, verbose=F)
head(s.markers, n=20) %>% rownames() -> labels
up.labels <- rownames(s.markers[s.markers$avg_logFC > FCcutoff,])
down.labels <- rownames(s.markers[s.markers$avg_logFC < -FCcutoff,])
labels <- unique(c(labels, up.labels,down.labels))
EnhancedVolcano(s.markers,
    lab = rownames(s.markers),
    x = 'avg_logFC',
    y = 'p_val_adj',
  #  xlim = c(-20, 20),
    title = 'MAIT (CD_inflamed vs. CD_uninflamed)',
    pCutoff = 0.05,
    FCcutoff = FCcutoff,
    pointSize = 3.0,
    labSize = 3.0,
   selectLab = labels,
    drawConnectors = T,
    widthConnectors = 0.2,
    colConnectors = 'grey30')
ggsave("figures/cdcpi3_compare/cdcpi3_DE_MAIT_CD_inflamed.pdf", width= 10, height= 8, units= "in")

###########################
## GSEA
###########################

x <- subset(s.a, orig.ident == "CD3_Tcell") # extract CPI dataset
s.markers <- FindMarkers(x, ident.1 ="CPI_colitis" ,group.by= "colitis", subset.ident="6", min.pct = 0.25, logfc.threshold = 0.01, verbose=F)

x1 <- subset(s.a, orig.ident == "CD_martin") # extract CD dataset
s.markers1 <- FindMarkers(x1, ident.1 ="CD_inflamed" ,group.by= "colitis", subset.ident="6", min.pct = 0.25, logfc.threshold = 0.01, verbose=F)

gene_list <- s.markers$avg_logFC[s.markers$p_val_adj <0.05]
names(gene_list) <- rownames(s.markers[s.markers$p_val_adj <0.05,])
gene_list <- sort(gene_list, decreasing = TRUE)


GO_file <- "/home/esk17/1.projects/colitis/data/GO_annotations/c7.all.v7.2.symbols.gmt"
res <- GSEA(gene_list, GO_file, pval= 0.05)
res$Plot
ggsave("figures/cdcpi3_compare/cdcpi3_treg_CPI_GSEA.pdf", width= 20, height= 10, units= "in")

######################
## QC comparison
######################
x <- subset(s.a, orig.ident == "CD3_Tcell") # extract CPI dataset
x1 <- subset(s.a, orig.ident == "CD_martin") # extract CD dataset

# nFeature_RNA split by dataset
p <- VlnPlot(s.a,features = "nFeature_RNA", pt.size= 0, group.by="seurat_clusters", split.by = "orig.ident") + theme(axis.text.x = element_text(angle = 90),legend.position = "right", axis.title.x=element_blank())
ggsave("figures/cdcpi3_compare/CPIvsCD_nFeatureRNA_by_cluster.pdf", width= 10, height= 7, units= "in")

# nUMI split by dataset
p <- VlnPlot(s.a,features = "nCount_RNA", pt.size= 0, group.by="seurat_clusters", split.by = "orig.ident") + theme(axis.text.x = element_text(angle = 90),legend.position = "right", axis.title.x=element_blank())
ggsave("figures/cdcpi3_compare/CPIvsCD_nUMI_by_cluster.pdf", width= 10, height= 7, units= "in")

# correlation between nUMI and nFeature_RNA
FeatureScatter(s.a,feature1="nFeature_RNA", feature2="nCount_RNA", group.by="orig.ident")
ggsave("figures/cdcpi3_compare/CPIvsCD_nUMIvsnFeature.pdf", width= 7, height= 7, units= "in")

##########################################
## Bar plot of colitis vs ctrl split by dataset
##########################################
x <- subset(s.a, orig.ident == "CD3_Tcell") # extract CPI dataset
x1 <- subset(s.a, orig.ident == "CD_martin") # extract CD dataset
x$colitis <- factor(x$colitis, levels=c("CPI_colitis","CPI_no-colitis","Control"))
# proportion of each dataset by cluster (harmony)
data.frame(group = x$colitis, cluster= x$seurat_clusters) %>%
	ggplot() + geom_bar(
		mapping = aes(x= cluster, fill= group),
		position = "fill"
	) + labs(y="proportion")+ theme_classic()
ggsave("figures/cdcpi3_compare/CPI_colitis_by_cluster.pdf", width= 10, height= 7, units= "in")
# proportion of each dataset by cluster (anchor)
data.frame(group = x1$colitis, cluster= x1$seurat_clusters) %>%
	ggplot() + geom_bar(
		mapping = aes(x= cluster, fill= group),
		position = "fill"
	) + labs(y="proportion")+ theme_classic()
ggsave("figures/cdcpi3_compare/CD_inflamed_by_cluster.pdf", width= 10, height= 7, units= "in")

#####################################################
## (11/30/20) Removing B cells
#####################################################
s <- readRDS("saved_objects/CD_martin_qc_100720/martin_naive_CD3_15100.rds")
date <- "113020"

# figure showing which clusters have the B cells
# violin plot of IGKC per cluster
p <- VlnPlot(s,features = "IGKC", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 90),legend.position = "none", axis.title.x=element_blank())
ggsave(paste("figures/CD_martin_", date, "/martin_naive_CD3_15100_VlnPlot_IGKC_by_cluster.pdf", sep=""), width= 7, height= 7, units= "in")
# violin plot of JCHAIN per cluster
p <- VlnPlot(s,features = "JCHAIN", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 90),legend.position = "none", axis.title.x=element_blank())
ggsave(paste("figures/CD_martin_", date, "/martin_naive_CD3_15100_VlnPlot_JCHAIN_by_cluster.pdf", sep=""), width= 7, height= 7, units= "in")
# feature plots the IGKC and JCHAIN split
FeaturePlot(s, features= c('IGKC','JCHAIN'), reduction="umap")
ggsave(paste("figures/CD_martin_", date, "/martin_naive_CD3_15100_Featureplot_IGKC_JCHAIN.pdf", sep=""), width= 14, height= 7, units= "in")
# Correlation between IGKC expression and CD3E
FeatureScatter(s,feature1="IGKC", feature2="CD3E")
ggsave(paste("figures/CD_martin_", date, "/martin_naive_CD3_15100_Scatter_IGKC_CD3E.pdf", sep=""), width= 7, height= 7, units= "in")
# feature plots the IGKC and JCHAIN split in cdcpi3 anchor
FeaturePlot(s.a, features= c('rna_IGKC','rna_JCHAIN'), reduction="umap")
ggsave(paste("figures/CD_martin_", date, "/cdcpi3a_Featureplot_IGKC_JCHAIN.pdf", sep=""), width= 14, height= 7, units= "in")
# Correlation between IGKC expression and CD3E in cdcpi3 anchor
FeatureScatter(s.a,feature1="rna_IGKC", feature2="rna_CD3E", group.by="orig.ident")
ggsave(paste("figures/CD_martin_", date, "/cdcpi3a_Scatter_IGKC_CD3E.pdf", sep=""), width= 7, height= 7, units= "in")
# Correlation between IGKC expression and JCHAIN in both
FeatureScatter(s.a,feature1="rna_IGKC", feature2="rna_JCHAIN", group.by="orig.ident")
ggsave(paste("figures/CD_martin_", date, "/cdcpi3a_vs_CD_Scatter_IGKC_JCHAIN.pdf", sep=""), width= 7, height= 7, units= "in")
# Select cells that are in the thresholds
potential.Bcells <- WhichCells(s.a, expression=rna_IGKC > 1 | rna_JCHAIN > 1)
#length(potential.Bcells)
# length(potential.Bcells)/dim(s.a)[2]
s.a$potential.Bcells <- colnames(s.a) %in% potential.Bcells
FeatureScatter(s.a,feature1="rna_IGKC", feature2="rna_JCHAIN", group.by = "potential.Bcells", cols=c("green","grey"))
ggsave(paste("figures/CD_martin_", date, "/cdcpi3a_vs_CD_Scatter_IGKC_JCHAIN_highlightBcells.pdf", sep=""), width= 7, height= 7, units= "in")

### Distributional stats of cells that are in the thresholds
# distribution by sample
data.frame(potential.Bcells = s.a$potential.Bcells, patient= s.a$sample.id) %>%
	ggplot() + geom_bar(
		mapping = aes(x= patient, fill= potential.Bcells),
		position = "fill") + labs(y="proportion")+ theme_classic() + coord_flip() +
	scale_fill_manual(values=c("grey","#F8766D"))
ggsave(paste("figures/CD_martin_", date, "/cdcpi3a_highlightBcells_per_patient.pdf", sep=""), width= 7, height= 7, units= "in")
# distribution by condition
data.frame(potential.Bcells = s.a$potential.Bcells, condition= s.a$colitis) %>%
	ggplot() + geom_bar(
		mapping = aes(x= condition, fill= potential.Bcells),
		position = "fill"
	) + labs(y="proportion")+ theme_classic() + scale_fill_manual(values=c("grey","#F8766D"))
ggsave(paste("figures/CD_martin_", date, "/cdcpi3a_highlightBcells_per_condition.pdf", sep=""), width= 7, height= 7, units= "in")

# apply negative filter to remove potential B cells
s.a <- s.a[,!s.a$potential.Bcells]
# save as new object
saveRDS(s.a, paste("saved_objects/CD_martin_qc_", date, "/cdcpi3a_Bminus.rds", sep=""))
# total number of cells per patient
# distribution by sample
data.frame(patient= s.a$sample.id) %>%
	ggplot() + geom_bar(mapping = aes(x= patient)) + labs(y="total no. of cells")+ theme_classic() + coord_flip()
ggsave(paste("figures/CD_martin_", date, "/cdcpi3aB_totalcount_per_patient.pdf", sep=""), width= 7, height= 7, units= "in")








#saveRDS(s, paste("saved_objects/CD_martin_qc_", date, "/martin_naive_CD3_15100.rds", sep=""))





