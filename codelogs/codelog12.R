############################################################
## Reinvestigation of the CD T cells
## edward_kim@college.harvard.edu - June 2020
#############################################################
​
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

######################################################
# 3 - after -> before label transfer
######################################################
s <- readRDS("saved_objects/martinTcells15200.rds")
s1 <- readRDS("saved_objects/cdcpi2_anchor.rds")
## Get annotations from Luoma-Suo cells
require(stringr)
# Finding the correct mapping: turns out it is the "-1_1, -2_1, etc." endings
# that match the Suo Tcells
colnames(s) %>% str_sub(17) %>% table()
colnames(s1) %>% str_sub(17) %>% table()

s$cdcpi2.annotations <- s1$seurat_clusters
# graph showing clusters: cdcpi2 -> martin15200
p <- DimPlot(s,reduction="umap",group.by="cdcpi2.annotations")
ggsave("figures/apply_cdcpi2_labels_to_martin15200.pdf")
p <- DimPlot(s,reduction="umap",group.by="cdcpi2.annotations", cells=!is.na(s$cdcpi2.annotations))
ggsave("figures/apply_cdcpi2_labels_to_martin15200_noNA.pdf")

######################################################
# 3 - CD3E and CD4/CD8 and other QC metrics per cluster
######################################################
# violin plot of CD3E per cluster
p <- VlnPlot(s,features = "CD3E", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 90),legend.position = "none", axis.title.x=element_blank())
ggsave("figures/martin15200_VlnPlot_CD3E_by_cluster.pdf", width= 7, height= 7, units= "in")
# feature plots the CD4 and CD8 split
FeaturePlot(s, features= c('CD4','CD8A','CD8B','CD3E'), reduction="umap")
ggsave("figures/martin15200_CD4xCD8_featureplot.pdf", width= 12, height= 12, units= "in")
# violin plot of %mito per cluster
p <- VlnPlot(s,features = "percent.mt", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 90),legend.position = "none", axis.title.x=element_blank())
ggsave("figures/martin15200_VlnPlot_percentmito_by_cluster.pdf", width= 7, height= 7, units= "in")
# violin plot of nFeature per cluster
p <- VlnPlot(s,features = "nFeature_RNA", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 90),legend.position = "none", axis.title.x=element_blank())
ggsave("figures/martin15200_VlnPlot_nFeatureRNA_by_cluster.pdf", width= 7, height= 7, units= "in")

################################################
# 4- DA analysis: CD inflamed vs CD non-inflamed
################################################
library(edgeR)
x <- s
x$CD_inf <- x$status == "inflamed"
x$sample.id <- x$orig.ident
design <- function(y) model.matrix(~factor(patient.id)+factor(CD_inf),y)
out <- DA_analysis(x, design, title = "martin15200_CD_inf_vs_uninf")
summary(out$y.ab$common.dispersion)
summary(out$fit.ab$var.prior)
summary(out$fit.ab$df.prior)
summary(decideTests(out$res))
topTags(out$res)

################################################
# 5- DE of martin15200
################################################
DE_heatmap(s,"figures/martin15200_CD3_DE.pdf")

################################################
# 6- Recluster martin from scratch
################################################
input.dirs <- dir(path="data/CD_martin",pattern="_")
for(i in 1:length(input.dirs)) {
	qc_CD(input.dirs[i])
}
################################################
# 7- Diagnostic plots
################################################
input.dirs <- dir(path="data/CD_martin",pattern="_")
save_figures_CD(input.dirs[1])
for(i in 2:length(input.dirs)) {
	save_figures_CD(input.dirs[i])
}

################################################
# 7- Table with total cells removed
################################################
input.dirs <- dir(path="data/CD_martin",pattern="_")
y <- sapply(input.dirs,function(x) get_info_CD(x,"092120"))
y.print <- as.data.frame(t(y))

library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")
xtable(y.print, display=c("s","d","e","d","d","f","d","d","d","d","d"))

################################################
# 6- Recluster martin from scratch (with new filter method)
################################################
input.dirs <- dir(path="data/CD_martin",pattern="_")
qc_CD(input.dirs[1], "100720")
save_figures_CD(input.dirs[1], "100720")

for(i in 2:length(input.dirs)) {
	qc_CD(input.dirs[i], "100720")
}
for(i in 2:length(input.dirs)) {
	save_figures_CD(input.dirs[i], "100720")
}

y <- sapply(input.dirs,function(x) get_info_CD(x,"100720"))
y.print <- as.data.frame(t(y))

library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")
xtable(y.print, display=c("s","d","e","d","d","f","d","d","d","d","d"))

################################################
# 6- cleaning clusters
################################################
input.dirs <- dir(path="data/CD_martin",pattern="_")
date <- "100720"
cells <- mclapply(paste("saved_objects/CD_martin_qc_", date, "/", input.dirs, "_3.RDS", sep=""),readRDS)
#### Removing bad clusters from each sample
cells[[1]] <- subset(cells[[1]], idents= c("5","6","7","8","9","10"), invert=TRUE)
cells[[2]] <- subset(cells[[2]], idents= c("11","12","6","7","5","15"), invert=TRUE)
cells[[4]] <- subset(cells[[4]], idents= c("3"), invert=TRUE)
cells[[5]] <- subset(cells[[5]], idents= c("12","10","11"), invert=TRUE)
cells[[6]] <- subset(cells[[6]], idents= c("5","6","7","8","10"), invert=TRUE)
cells[[7]] <- subset(cells[[7]], idents= c("6","7","9","10"), invert=TRUE)
cells[[8]] <- subset(cells[[8]], idents= c("6","7","8","9"), invert=TRUE)
cells[[9]] <- subset(cells[[9]], idents= c("10","11"), invert=TRUE)
cells[[10]] <- subset(cells[[10]], idents= c("8","5","13","11","7"), invert=TRUE)
cells[[11]] <- subset(cells[[11]], idents= c("4","13","14"), invert=TRUE)
cells[[12]] <- subset(cells[[12]], idents= c("4","5","7","11","12","13"), invert=TRUE)
cells[[13]] <- subset(cells[[13]], idents= c("8","11","13"), invert=TRUE)
cells[[14]] <- subset(cells[[14]], idents= c("5","10","11"), invert=TRUE)
cells[[15]] <- subset(cells[[15]], idents= c("3","4","5","10","11","12","13","14"), invert=TRUE)
cells[[16]] <- subset(cells[[16]], idents= c("9","10","11"), invert=TRUE)
cells[[17]] <- subset(cells[[17]], idents= c("7","9","12"), invert=TRUE) # for 17, I am concerned about the quality of cluster 0
cells[[18]] <- subset(cells[[18]], idents= c("5","8","10","11","12"), invert=TRUE)
cells[[19]] <- subset(cells[[19]], idents= c("1","11","12"), invert=TRUE)
cells[[20]] <- subset(cells[[20]], idents= c("7","9","10","11"), invert=TRUE)
cells[[21]] <- subset(cells[[21]], idents= c("7","8","9","12","13","14","15"), invert=TRUE)
cells[[22]] <- subset(cells[[22]], idents= c("5","6","10","11","12","13"), invert=TRUE)
#### Getting final cell counts:
#sapply(cells,dim)[2,]

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
#### mutate orig.ident and edit sample.id
s$orig.ident <- ifelse(grepl(s$orig.ident,pattern="^GSM*"), "CD_martin", s$orig.ident)
s$sample.id <- ifelse(s$orig.ident=="CD_martin", str_glue("CD_{str_sub(s$sample.id,-3)}"), s$sample.id)
#### add patient.id
s <- add_patient_id_all(s)
s$patient.id <- str_replace(s$patient.id,"Patient ","P")
s$patient.id <- factor(s$patient.id, levels = c("P5","P6","P7","P8","P10","P11","P12","P13","P14","P15","P16"))
# add inflamed id for Crohns
inflamed <-c("CD__69$","CD_122$","CD_128$","CD_138$","CD_158$","CD_181$","CD_187$","CD_193$","CD_190$","CD_196$","CD_209$")
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
#33694x86308 to 20450x86308 
x <- s[["RNA"]]@counts
nnz_row <- tabulate(x@i + 1)
keep <- rownames(x)[nnz_row>10]
s <- subset(s, features = keep)

#############################
# 4 - perform clustering
#############################
s <- PercentageFeatureSet(s, pattern = "^MT-", col.name = "percent.mt")
date <- "100720"
stats <- filter_stats(s, save=T, filename=paste("saved_objects/CD_martin_qc_", date, "/martin_naive_intg_filterstats.RDS", sep=""))

#custom removal based on LNFxMITO plot
#s$cells.to.remove <-  s$nFeature_RNA > exp(median(stats$lnf.pre)+3*mad(stats$lnf.pre))

# 869 cells removed; nFeature_RNA > 2177
#s <- s[,!s$cells.to.remove]
s <- NormalizeData(s, normalization.method = "LogNormalize", scale.factor = 10000) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000) ## Highly variable genes
hvf <- VariableFeatures(s)
'%ni%' = Negate('%in%')
all.genes <- rownames(s)
trvgenes <- all.genes[grepl(x=all.genes, pattern = "^TRAV|^TRBV|^TRGV|^TRDV")]
igvgenes <- all.genes[grepl(x=all.genes, pattern = "^IG.V")]
VariableFeatures(s) <- hvf[hvf %ni% c(trvgenes,igvgenes)]
# Do the rest
s <- ScaleData(s, features = all.genes, vars.to.regress=c("percent.mt","nFeature_RNA")) %>% RunPCA(features=VariableFeatures(s))
s <-  cluster_umap(s, npc=20, k.param=120, resolution=0.5, min.dist=0.3)
DimPlot(s, label=T)
ggsave("figures/CD_martin_100720_integration/martin_naive_npc20_k120_umap.pdf")
DimPlot(s, group.by="orig.ident" , label=F)
ggsave("figures/CD_martin_100720_integration/martin_naive_npc20_k120_persample.pdf")
# distribution by sample
data.frame(cluster = s$seurat_clusters, patient= s$orig.ident) %>%
	ggplot() + geom_bar(
		mapping = aes(x= patient, fill= cluster),
		position = "fill"
	) + labs(y="proportion")+ theme_classic() + coord_flip()
ggsave("figures/CD_martin_100720_integration/martin_naive_npc20_k120_distpersample.pdf")
## Replicating supplemental figure for QC validation
library(ggplot2)
df <- as.data.frame(table(s$patient.id,s$colitis))
colnames(df) <- c("patient","sampletype", "no.of.cells")
ggplot(df, aes(patient, no.of.cells, fill = sampletype)) + 
  geom_bar(stat="identity",position="dodge") + 
  scale_fill_brewer(palette = "Set1") +
  theme_classic()
ggsave("figures/CD_martin_100720_integration/martin_naive_npc20_k120_counts_per_patient.pdf")


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
ggsave("figures/CD_martin_100720_integration/martin_naive_npc20_k120_singleR.pdf", width= 24, height= 6, units= "in")
######################################################
# 3 - CD3E and CD4/CD8 and other QC metrics per cluster
######################################################
# violin plot of CD3E per cluster
p <- VlnPlot(s,features = "CD3E", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 90),legend.position = "none", axis.title.x=element_blank())
ggsave("figures/CD_martin_100720_integration/martin_naive20120_VlnPlot_CD3E_by_cluster.pdf", width= 7, height= 7, units= "in")
# feature plots the CD4 and CD8 split
FeaturePlot(s, features= c('CD4','CD8A','CD8B','CD3E'), reduction="umap")
ggsave("figures/CD_martin_100720_integration/martin_naive20120_CD4xCD8_featureplot.pdf", width= 12, height= 12, units= "in")
# violin plot of %mito per cluster
p <- VlnPlot(s,features = "percent.mt", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 90),legend.position = "none", axis.title.x=element_blank())
ggsave("figures/CD_martin_100720_integration/martin_naive20120_VlnPlot_percentmito_by_cluster.pdf", width= 7, height= 7, units= "in")
# violin plot of nFeature per cluster
p <- VlnPlot(s,features = "nFeature_RNA", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 90),legend.position = "none", axis.title.x=element_blank())
ggsave("figures/CD_martin_100720_integration/martin_naive20120_VlnPlot_nFeatureRNA_by_cluster.pdf", width= 7, height= 7, units= "in")

date <- "100720"
saveRDS(s, paste("saved_objects/CD_martin_qc_", date, "/martin_naive_clustering_20120.rds", sep=""))

######################################################
# 3 - isolate CD3 T cells
######################################################
# violin plot of CD3E per cluster
s1 <- s
s <- subset(s1, idents=c("0","1","2","17"))
s <- NormalizeData(s, normalization.method = "LogNormalize", scale.factor = 10000) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000) 
## Highly variable genes
hvf <- VariableFeatures(s)
'%ni%' = Negate('%in%')
all.genes <- rownames(s)
trvgenes <- all.genes[grepl(x=all.genes, pattern = "^TRAV|^TRBV|^TRGV|^TRDV")]
igvgenes <- all.genes[grepl(x=all.genes, pattern = "^IG.V")]
VariableFeatures(s) <- hvf[hvf %ni% c(trvgenes,igvgenes)]
# Do the rest
nPCs <- 20
s <- ScaleData(s, features = all.genes, vars.to.regress=c("percent.mt","nFeature_RNA")) %>% RunPCA(features=VariableFeatures(s))
s <-  cluster_umap(s, npc=20, k.param=100, resolution=0.5, min.dist=0.3)
DimPlot(s, label=T)
ggsave("figures/CD_martin_100720_integration/martin_naive_CD3_20100_umap.pdf")
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
ggsave("figures/CD_martin_100720_integration/martin_naive_CD3_20100_singleR.pdf", width= 24, height= 6, units= "in")


######################################################
# 3 - CD3E and CD4/CD8 and other QC metrics per cluster
######################################################
# violin plot of CD3E per cluster
p <- VlnPlot(s,features = "CD3E", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 90),legend.position = "none", axis.title.x=element_blank())
ggsave("figures/CD_martin_100720_integration/martin_naive_CD3_20100_VlnPlot_CD3E_by_cluster.pdf", width= 7, height= 7, units= "in")
# feature plots the CD4 and CD8 split
FeaturePlot(s, features= c('CD4','CD8A','CD8B','CD3E'), reduction="umap")
ggsave("figures/CD_martin_100720_integration/martin_naive_CD3_20100_CD4xCD8_featureplot.pdf", width= 12, height= 12, units= "in")
# violin plot of %mito per cluster
p <- VlnPlot(s,features = "percent.mt", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 90),legend.position = "none", axis.title.x=element_blank())
ggsave("figures/CD_martin_100720_integration/martin_naive_CD3_20100_VlnPlot_percentmito_by_cluster.pdf", width= 7, height= 7, units= "in")
# violin plot of nFeature per cluster
p <- VlnPlot(s,features = "nFeature_RNA", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 90),legend.position = "none", axis.title.x=element_blank())
ggsave("figures/CD_martin_100720_integration/martin_naive_CD3_20100_VlnPlot_nFeatureRNA_by_cluster.pdf", width= 7, height= 7, units= "in")


#########################
## Repeat clustering
#########################
# Do the rest
s <-  cluster_umap(s, npc=15, k.param=100, resolution=0.5, min.dist=0.3)
DimPlot(s, label=T)
ggsave("figures/CD_martin_100720_integration/martin_naive_CD3_15100_umap.pdf")

# distribution by sample
data.frame(cluster = s1$seurat_clusters, patient= s1$orig.ident) %>%
	ggplot() + geom_bar(
		mapping = aes(x= patient, fill= cluster),
		position = "fill"
	) + labs(y="proportion")+ theme_classic() + coord_flip()
ggsave("figures/CD_martin_100720_integration/martin_naive_CD3_15100_distpersample.pdf")
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
ggsave("figures/CD_martin_100720_integration/martin_naive_CD3_15100_singleR.pdf", width= 24, height= 6, units= "in")
######################################################
# 3 - CD3E and CD4/CD8 and other QC metrics per cluster
######################################################
# violin plot of CD3E per cluster
p <- VlnPlot(s,features = "CD3E", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 90),legend.position = "none", axis.title.x=element_blank())
ggsave("figures/CD_martin_100720_integration/martin_naive_CD3_15100_VlnPlot_CD3E_by_cluster.pdf", width= 7, height= 7, units= "in")
# feature plots the CD4 and CD8 split
FeaturePlot(s, features= c('CD4','CD8A','CD8B','CD3E'), reduction="umap")
ggsave("figures/CD_martin_100720_integration/martin_naive_CD3_15100_CD4xCD8_featureplot.pdf", width= 12, height= 12, units= "in")
# violin plot of %mito per cluster
p <- VlnPlot(s,features = "percent.mt", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 90),legend.position = "none", axis.title.x=element_blank())
ggsave("figures/CD_martin_100720_integration/martin_naive_CD3_15100_VlnPlot_percentmito_by_cluster.pdf", width= 7, height= 7, units= "in")
# violin plot of nFeature per cluster
p <- VlnPlot(s,features = "nFeature_RNA", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 90),legend.position = "none", axis.title.x=element_blank())
ggsave("figures/CD_martin_100720_integration/martin_naive_CD3_15100_VlnPlot_nFeatureRNA_by_cluster.pdf", width= 7, height= 7, units= "in")
# differential gene analysis
DE_heatmap(s1,"figures/CD_martin_100720_integration/martin_naive_15100_DE_avgexp.pdf")


date <- "100720"
saveRDS(s, paste("saved_objects/CD_martin_qc_", date, "/martin_naive_CD3_15100.rds", sep=""))

#####################################################
## Replicating supplemental figure for QC validation (T cell only)
#####################################################
library(ggplot2)
df <- as.data.frame(table(s$patient.id,s$colitis))
colnames(df) <- c("patient","sampletype", "no.of.cells")
ggplot(df, aes(patient, no.of.cells, fill = sampletype)) + 
  geom_bar(stat="identity",position="dodge") + 
  scale_fill_brewer(palette = "Set1") +
  theme_classic()
ggsave("figures/CD_martin_100720_integration/martin_naive_15100_counts_per_patient.pdf")






