############################################################
## Harmony Integration CD (15100) + CPI -> cdcpi3_harmony
## edward_kim@college.harvard.edu - June 2020
#############################################################

##############################
# 0 - Load librairies
##############################
library(dplyr)
library(ggplot2)
library(Seurat)
library(future)
library(future.apply)
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

####################################
# 3 - Calculate Harmony Embeddings
####################################
library(harmony)
s1 <- readRDS("saved_objects/CD_martin_qc_100720/martin_naive_CD3_15100.rds")
s2 <- readRDS("saved_objects/onlyTcells1550_labeled.rds")
s <- merge(s1,s2)
# removing low frequency features
x <- s[["RNA"]]@counts
nnz_row <- tabulate(x@i + 1)
keep <- rownames(x)[nnz_row>10]
s <- subset(s, features = keep)
#### Adding/Editing meta data
library(stringr)
# removing unapplicable meta.data columns
m.names <- colnames(s@meta.data)
columns.to.remove <- m.names[grepl("^pANN*|^DF.c*|seurat_clusters|^RNA*|cells.to.remove|Single*|percent.mito",m.names)]
for(i in columns.to.remove) {
  s[[i]] <- NULL
}
# rename sampletype to sample.id and merge with orig.ident material
s$sample.id <- s$sampletype
s$sampletype <- NULL
s$sample.id <- ifelse(is.na(s$sample.id), s$orig.ident ,s$sample.id)
#### mutate orig.ident and edit sample.id
s$orig.ident <- ifelse(grepl(s$orig.ident,pattern="^GSM*"), "CD_martin", s$orig.ident)
s$sample.id <- ifelse(s$orig.ident=="CD_martin", str_glue("CD_{str_sub(s$sample.id,-3)}"), s$sample.id)
s$sample.id <- ifelse(s$orig.ident=="CD3_Tcell", str_glue("CPI_{s$sample.id}"), s$sample.id)
#### add and merge patient.id
s <- add_patient_id_all(s)
s$patient.id <- str_replace(s$patient.id,"Patient ","CD_P")
s$patient.id <- ifelse(s$orig.ident == "CD3_Tcell",str_glue("CPI_{s$sample.id}"),s$patient.id)
# add inflamed id for Crohns
inflamed <-c("CD__69$","CD_122$","CD_128$","CD_138$","CD_158$","CD_181$","CD_187$","CD_193$","CD_190$","CD_196$","CD_209$")
toMatch <- paste(inflammed,collapse="|")
status <- ifelse(grepl(s$sample.id, pattern=toMatch),"CD_inflamed","CD_uninflamed")
s$colitis <- ifelse(s$orig.ident == "CD_martin", status ,s$colitis)
s$colitis <- ifelse(s$colitis =="Colitis", "CPI_colitis", s$colitis)
s$colitis <- ifelse(s$colitis =="No-Colitis", "CPI_no-colitis", s$colitis)
s$colitis2 <- ifelse(s$orig.ident == "CD_martin", status ,s$colitis2)
## Convert colitis2 and sample.id to factors
s$colitis2 <- factor(s$colitis2, levels = c("CD_inflamed","CD_uninflamed", "Colitis (aPD1)" ,"No-Colitis (aPD1)" , "Colitis (combo)" ,  "No-Colitis (combo)", "Control"))
conditions <- levels(s$colitis2)
order <- vector(mode="character", length=0)
for(i in 1:length(conditions)) {
    ids <- as.character(s$sample.id)
    order <- c(order, rownames(table(ids[s$colitis2 == conditions[i]])))
}
s$sample.id <- factor(s$sample.id, levels = order)

#### Clustering
s <- NormalizeData(s, normalization.method = "LogNormalize", scale.factor = 10000) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000) ## Highly variable genes
hvf <- VariableFeatures(s)
'%ni%' = Negate('%in%')
all.genes <- rownames(s)
trvgenes <- all.genes[grepl(x=all.genes, pattern = "^TRAV|^TRBV|^TRGV|^TRDV")]
igvgenes <- all.genes[grepl(x=all.genes, pattern = "^IG.V")]
VariableFeatures(s) <- hvf[hvf %ni% c(trvgenes,igvgenes)]
# Do the rest
s <- ScaleData(s, features = all.genes) %>% RunPCA(features=VariableFeatures(s))

# Harmony part
options(repr.plot.height = 2.5, repr.plot.width = 6)
s <- s %>% RunHarmony("orig.ident", plot_convergence = TRUE)

#save Seurat object with harmony embeddings
saveRDS(s, "saved_objects/cdcpi3_harmony.rds")

#################################################################
# 4 - Basic Plots
#################################################################
# Plot Harmony embedding
options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = s, reduction = "harmony", pt.size = .1, group.by = "orig.ident")
p2 <- VlnPlot(object = s, features = "harmony_1", group.by = "orig.ident", pt.size = .1)
p1 + p2
ggsave("figures/cdcpi3_harmony/cdcpi3_plot.pdf", width= 14, height= 7, units= "in")


# Downstream
s <- s %>% 
    RunUMAP(reduction = "harmony", dims = 1:15) %>% 
    FindNeighbors(reduction = "harmony", dims = 1:15) %>% 
    FindClusters(resolution = 0.5) %>% 
    identity()

#plot UMAP split by dataset
options(repr.plot.height = 4, repr.plot.width = 10)
DimPlot(s, reduction = "umap", group.by = "seurat_clusters", pt.size = .1, split.by = 'orig.ident')
ggsave("figures/cdcpi3_harmony/cdcpi3_UMAP_splitby_dataset.pdf", width= 14, height= 7, units= "in")
#plot UMAP
DimPlot(s, reduction = "umap", label = TRUE, pt.size = .1)
ggsave("figures/cdcpi3_harmony/cdcpi3_UMAP.pdf", width= 7, height= 7, units= "in")


#################################################################
# 10 - Clustering check with comparison to within-batch clusters
#################################################################
library(pheatmap)
s<- readRDS("saved_objects/cdcpi3_harmony.rds")
s1 <- readRDS("saved_objects/CD_martin_qc_100720/martin_naive_CD3_15100.rds")
s2 <- readRDS("saved_objects/onlyTcells1550_labeled.rds")
# For the first batch (adding +10 for a smoother color transition
# from zero to non-zero counts for any given matrix entry).
tab <- table(paste("after", s$seurat_clusters[s$orig.ident=="CD_martin"]),
    paste("before", s1$seurat_clusters))
heatCD <- pheatmap(matrix.sort.no.diag(log10(tab+10)), cluster_row=FALSE, cluster_col=FALSE,
    main="CD comparison", silent=TRUE)

# For the second batch.
tab <- table(paste("after", s$seurat_clusters[s$orig.ident=="CD3_Tcell"]),
    paste("before", s2$seurat_clusters))
heatCPI <- pheatmap(matrix.sort.no.diag.T(log10(tab+10)), cluster_row=FALSE, cluster_col=FALSE,
    main="CPI comparison", silent=TRUE)

pdf("figures/cdcpi3_harmony/cdcpi3h_cluster_comparison.pdf")
gridExtra::grid.arrange(heatCD[[4]], heatCPI[[4]])
dev.off()

####################################
# 5 - Basic Vizualizations
####################################
# violin plot of %mito per cluster
p <- VlnPlot(s,features = "percent.mt", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 90),legend.position = "none", axis.title.x=element_blank())
ggsave("figures/cdcpi3_harmony/cdcpi3h_VlnPlot_percentmito_by_cluster.pdf", width= 7, height= 7, units= "in")
# violin plot of nFeature per cluster
p <- VlnPlot(s,features = "nFeature_RNA", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 90),legend.position = "none", axis.title.x=element_blank())
ggsave("figures/cdcpi3_harmony/cdcpi3h_VlnPlot_nFeatureRNA_by_cluster.pdf", width= 7, height= 7, units= "in")
# feature plots the CD4 and CD8 split
FeaturePlot(s, features= c('CD4','CD8A','CD8B','CD3E'), reduction="umap")
ggsave("figures/cdcpi3_harmony/cdcpi3h_CD4xCD8_featureplot.pdf", width= 12, height= 12, units= "in")
# violin plot of CD4 per cluster
p <- VlnPlot(s,features = "CD4", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 90),legend.position = "none", axis.title.x=element_blank())
ggsave("figures/cdcpi3_harmony/cdcpi3h_VlnPlot_CD4_by_cluster.pdf", width= 7, height= 7, units= "in")
# violin plot of CD8 per cluster
p <- VlnPlot(s,features = "CD8A", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 90),legend.position = "none", axis.title.x=element_blank())
ggsave("figures/cdcpi3_harmony/cdcpi3h_VlnPlot_CD8_by_cluster.pdf", width= 7, height= 7, units= "in")
# violin plot of CD3E per cluster
p <- VlnPlot(s,features = "CD3E", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 90),legend.position = "none", axis.title.x=element_blank())
ggsave("figures/cdcpi3_harmony/cdcpi3h_VlnPlot_CD3E_by_cluster.pdf", width= 7, height= 7, units= "in")
# distribution by sample
data.frame(cluster = s$seurat_clusters, patient= s$sample.id) %>%
	ggplot() + geom_bar(
		mapping = aes(x= patient, fill= cluster),
		position = "fill"
	) + labs(y="proportion")+ theme_classic() + coord_flip()
ggsave("figures/cdcpi3_harmony/cdcpi3h_cluster_prop_per_sample.pdf",width= 7, height= 7, units= "in")
# differential gene analysis
DE_heatmap(s,"figures/cdcpi3_harmony/cdcpi3h_DE_avgexp.pdf")

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

#saveRDS(list(pred.CPI,pred.Martin,pred.Monaco),"saved_objects/label_info_for_cdcpi21550.rds")


require(ggpubr)
p1 <- DimPlot(s, reduction = "umap", label=T, group.by="SingleR.Martin") + labs(title="Martin reference")
p2 <- DimPlot(s, reduction = "umap", label=T, group.by="SingleR.Luoma") + labs(title="Luoma-Suo reference")
p3 <- DimPlot(s, reduction = "umap", label=T, group.by="SingleR.Monaco") + labs(title="Monaco reference")
ggarrange(p1,p2,p3, ncol=3, nrow=1)
ggsave("figures/cdcpi3_harmony/cdcpi3h_labeling.pdf", width= 24, height= 6, units= "in")




