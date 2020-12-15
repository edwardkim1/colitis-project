##############################
# 0 - Load librairies
##############################
library(dplyr)
library(ggplot2)
library(Seurat)
library(future)
library(future.apply)
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

######################################
# 3 - Preparing Martin for integration
######################################
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
s <- readRDS("saved_objects/martinTcells.rds")
s$sampletype <- s$orig.ident
s1 <- readRDS("saved_objects/onlyTcells1550_labeled.rds")

####################################
# 4 - Merging data
####################################
s <- merge(s, s1)
# add inflamed id for Crohns
inflammed <-c("_1$","_3$","_6$","_7$","_10$","_12$","_14$","_16$","_18$")
toMatch <- paste(inflammed,collapse="|")
status <- ifelse(grepl(colnames(s), pattern=toMatch),"CD_inflamed","CD_uninflamed")
s$colitis <- ifelse(is.na(s$colitis), status ,s$colitis)
s$colitis <- ifelse(s$colitis =="Colitis", "CPI_colitis", s$colitis)
s$colitis <- ifelse(s$colitis =="No-Colitis", "CPI_no-colitis", s$colitis)
s$colitis2 <- ifelse(is.na(s$colitis2), status ,s$colitis2)
table(s$colitis2)
# rename patient.id to CD.patient.id
s$CD.patient.id <- s$patient.id
s$patient.id <- NULL
# rename sampletype to sample.id and merge with orig.ident material
s$sample.id <- s$sampletype
s$sampletype <- NULL
s$sample.id <- ifelse(is.na(s$sample.id), s$orig.ident ,s$sample.id)
table(s$sample.id)
# removing unapplicable meta.data columns
m.names <- colnames(s@meta.data)
columns.to.remove <- m.names[grepl("^pANN*|^DF.c*|seurat_clusters|^RNA*|cells.to.remove|Single*|percent.mito",m.names)]
for(i in columns.to.remove) {
  s[[i]] <- NULL
}
#24385x119063 to 21237x119063
x <- s[["RNA"]]@counts
nnz_row <- tabulate(x@i + 1)
keep <- rownames(x)[nnz_row>10]
s <- subset(s, features = keep)

####################################
# 5 - Clustering
####################################
s <- PercentageFeatureSet(s, pattern = "^MT-", col.name = "percent.mt")
#stats <- filter_stats(s, save=T, filename="saved_objects/allTcells_filterstats.RDS")
#s$cells.to.remove <- stats$cells.to.remove
#saveRDS(s,"saved_objects/allTcells_1.RDS")
#s <- s[,!s$cells.to.remove]
s <- NormalizeData(s, normalization.method = "LogNormalize", scale.factor = 10000) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000) ## Highly variable genes
hvf <- VariableFeatures(s)
'%ni%' = Negate('%in%')
all.genes <- rownames(s)
trvgenes <- all.genes[grepl(x=all.genes, pattern = "^TRAV|^TRBV|^TRGV|^TRDV")]
igvgenes <- all.genes[grepl(x=all.genes, pattern = "^IG.V")]
VariableFeatures(s) <- hvf[hvf %ni% c(trvgenes,igvgenes)]
# Do the rest
s <- ScaleData(s, features = all.genes, vars.to.regress=c("percent.mt")) %>% RunPCA(features=VariableFeatures(s))
s <-  cluster_umap(s, npc=15, k.param=50, resolution=0.5, min.dist=0.3)

##############################
# 6 - Diagnostic Plots
##############################
## Figures from QC
stats <- readRDS("saved_objects/cdcpi_filterstats.RDS")

# Mito% Pre/Post Violin Plot
p <- ggplot(data=data.frame(x= c(rep("pre mito%", length(stats$mt.pre)),rep("post mito%", length(stats$mt.post))), percent.mito = c(stats$mt.pre,stats$mt.post))) + geom_violin(aes(x= x, y=percent.mito))
ggsave("figures/cdcpi_MITOviolin.pdf")

# Log(No. of Features) Pre/Post Violin Plot
p <- ggplot(data=data.frame(x= c(rep("pre lnf", length(stats$lnf.pre)),rep("post lnf", length(stats$lnf.post))), log.num.feats = c(stats$lnf.pre,stats$lnf.post))) + geom_violin(aes(x= x, y=log.num.feats))
ggsave("figures/cdcpi_LNFviolin.pdf")

# Log.no.features vs. mito% Figure
s$cells.to.remove <-  s$nFeature_RNA > exp(median(stats$lnf.pre)+3*mad(stats$lnf.pre)) #greater than 2177.043 features
p <- FeatureScatter(s,feature1="nFeature_RNA", feature2="percent.mt", group.by="cells.to.remove")
ggsave("figures/cdcpi_LNFxMITO.pdf")

# Mito% Model Fitting Figure
p <- plot_kde(stats$mt.pre, kernel = "gaussian", bw=0.5, lab.x = "mitochondrial gene percentage (per cell)", overlay = T, model="Norm", dmodel = function(x) {dnorm(x,stats$mt.median,stats$mt.mad)}, 
	support=seq(0,100,0.1))
p + geom_vline(xintercept= stats$mt.lim, color = "red") +
	geom_vline(xintercept= stats$mt.median-3*stats$mt.mad, linetype = "dotted") + 
	geom_vline(xintercept= stats$mt.median+3*stats$mt.mad, linetype = "dotted")
ggsave("figures/cdcpi_MITO-KDE-Normal.pdf")

# Log.no.features Model Fitting Figure
p <- plot_kde(stats$lnf.pre, kernel = "gaussian", bw=0.01, lab.x = "Log(No. of Features)", overlay = T, model="Norm", dmodel = function(x) {dnorm(x,stats$lnf.median,stats$lnf.mad)})
p + geom_vline(xintercept= stats$lnf.lim, color = "red") +
	geom_vline(xintercept= stats$lnf.median-3*stats$lnf.mad, linetype = "dotted") + 
	geom_vline(xintercept= stats$lnf.median+3*stats$lnf.mad, linetype = "dotted")
ggsave("figures/cdcpi_LNF-KDE-Normal.pdf")


## Figures from Object 2
s <- readRDS("saved_objects/cdcpi_2.RDS")

# Variable Features Save & Figure
VariableFeatures(s)%>%saveRDS(file="saved_objects/cdcpi_hvf.RDS")
p <- VariableFeaturePlot.Tcells(s) %>% LabelPoints(points=head(VariableFeatures(s),10), repel= TRUE)
ggsave("figures/cdcpi_VariableFeatures.pdf")

# PCA Elbow Plot Figure
p <- ElbowPlot(s, ndims= 30)
ggsave("figures/cdcpi_elbowplot.pdf")

# Evaling cluster stability: Bootstrap
nPCs<- 20
myClusterFUN <- function(seurat.object) {
    g <- seurat.object %>% FindNeighbors(dims = 1:nPCs, verbose=F) %>% FindClusters(resolution = 0.5, verbose=F)
    g$seurat_clusters
}
originals <- s$seurat_clusters
set.seed(111)
coassign <- bootstrapCluster(s, FUN=myClusterFUN, clusters=originals, iterations=50)
pdf("figures/cdcpi_ClusterStability_50iter.pdf")
pheatmap(coassign, cluster_row=FALSE, cluster_col=FALSE, color=rev(viridis::magma(100)))
dev.off()

# Clustering Figure
p <- DimPlot(s, label=T)
ggsave("figures/cdcpi_UMAP.pdf")

#####################################
## Basic Inference: cdcpi
#####################################
## Viz 1: UMAP - group by colitis and colitis2
p <- DimPlot(s, reduction = "umap", label=F, group.by="colitis")
q <- DimPlot(s, reduction = "umap", label=F, group.by="colitis2")
p+q
ggsave("figures/cdcpi_group_colitis.pdf",width= 16, height= 6, units= "in")




