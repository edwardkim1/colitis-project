############################################################
## Anchor redo with the CD + CPI -> cdcpi2
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
# 3 - Integration
##############################
s1 <- readRDS("saved_objects/martinTcells15200.rds")
s2 <- readRDS("saved_objects/onlyTcells1550_labeled.rds")
s.list <- list(CPI=s2,CD=s1)
# variable features
s.features <- SelectIntegrationFeatures(object.list = s.list, nfeatures = 2000)
'%ni%' = Negate('%in%')
trvgenes <- s.features[grepl(x=s.features, pattern = "^TRAV|^TRBV|^TRGV|^TRDV")]
igvgenes <- s.features[grepl(x=s.features, pattern = "^IG.V")]
s.features <- s.features[s.features %ni% c(trvgenes,igvgenes)]
# Finding the integration anchors and integrating
reference.number <- 1
s.anchors <- FindIntegrationAnchors(object.list = s.list, 
	reference = reference.number, 
	normalization.method= "LogNormalize",
	anchor.features= s.features,
	dims = 1:20)
s <- IntegrateData(anchorset = s.anchors)

#### Observation
# length(unique(c(rownames(s.list$CPI), rownames(s.list$CD))))
# > 24385
# is different from
# 26385 x 119063 integrated matrix (26385 # of features)

####################################
# 4 - Adding/Editing meta data
####################################
library(stringr)
# rename sampletype to sample.id and merge with orig.ident material
s$sample.id <- s$sampletype
s$sampletype <- NULL
s$sample.id <- ifelse(is.na(s$sample.id), s$orig.ident ,s$sample.id)
#### mutate orig.ident and edit sample.id
s$orig.ident <- ifelse(grepl(s$orig.ident,pattern="^GSM*"), "CD_martin", s$orig.ident)
s$sample.id <- ifelse(s$orig.ident=="CD_martin", str_glue("CD_{str_sub(s$sample.id,-3)}"), s$sample.id)
s$sample.id <- ifelse(s$orig.ident=="CD3_Tcell", str_glue("CPI_{s$sample.id}"), s$sample.id)
#### merge patient.id
s$patient.id <- str_replace(s$patient.id,"Patient ","CD_P")
s$patient.id <- ifelse(s$orig.ident == "CD3_Tcell",str_glue("CPI_{s$sample.id}"),s$patient.id)
# add inflamed id for Crohns
inflammed <-c("_1$","_3$","_6$","_7$","_10$","_12$","_14$","_16$","_18$")
toMatch <- paste(inflammed,collapse="|")
status <- ifelse(grepl(colnames(s), pattern=toMatch),"CD_inflamed","CD_uninflamed")
s$colitis <- ifelse(is.na(s$colitis), status ,s$colitis)
s$colitis <- ifelse(s$colitis =="Colitis", "CPI_colitis", s$colitis)
s$colitis <- ifelse(s$colitis =="No-Colitis", "CPI_no-colitis", s$colitis)
s$colitis2 <- ifelse(is.na(s$colitis2), status ,s$colitis2)
## Convert colitis2 and sample.id to factors
s$colitis2 <- factor(s$colitis2, levels = c("CD_inflamed","CD_uninflamed", "Colitis (aPD1)" ,"No-Colitis (aPD1)" , "Colitis (combo)" ,  "No-Colitis (combo)", "Control"))
conditions <- levels(s$colitis2)
order <- vector(mode="character", length=0)
for(i in 1:length(conditions)) {
	order <- c(order, rownames(table(s$sample.id[s$colitis2 == conditions[i]])))
}
s$sample.id <- factor(s$sample.id, levels = order)
# removing unapplicable meta.data columns
m.names <- colnames(s@meta.data)
columns.to.remove <- m.names[grepl("^pANN*|^DF.c*|seurat_clusters|^RNA*|cells.to.remove|Single*|percent.mito",m.names)]
for(i in columns.to.remove) {
  s[[i]] <- NULL
}
#26385x119063 to 23202x119063
x <- s[["RNA"]]@counts
nnz_row <- tabulate(x@i + 1)
keep <- rownames(x)[nnz_row>10]
s <- subset(s, features = keep)

####################################
# 5 - Clustering
####################################
DefaultAssay(s) <- "integrated"
# Do the rest; no need to regress out percent.mt because those genes are not in the list
s <- ScaleData(s) %>% RunPCA()
# PCA Elbow Plot Figure
p <- ElbowPlot(s, ndims= 30)
ggsave("figures/cdcpi2_anchor_elbowplot.pdf")

s <- FindNeighbors(s, dims = 1:15, k.param = 50) 
# names(s@graphs)
s <- FindClusters(s,resolution = 0.5, graph.name ="integrated_snn") 
s <- RunUMAP(s,min.dist = 0.3, graph="integrated_snn")

#visualize
p <- DimPlot(s, label=T)
ggsave("figures/cdcpi2_anchor_PCA.pdf")

p <- DimPlot(s, reduction="umap", label=T)
ggsave("figures/cdcpi2_anchor_UMAP.pdf")

## Viz 1: UMAP - group by colitis and colitis2
p <- DimPlot(s, reduction = "umap", label=F, group.by="colitis")
q <- DimPlot(s, reduction = "umap", label=F, group.by="colitis2")
p+q
ggsave("figures/cdcpi2_anchor_group_colitis.pdf",width= 16, height= 6, units= "in")

## Viz 4a: cluster proportions
table(s$seurat_clusters,s$colitis2) %>% unclass() %>% prop.table()
data.frame(cluster =s$seurat_clusters, condition= s$colitis) %>%
	ggplot() + geom_bar(
		mapping = aes(x= cluster, fill= condition),
		position = "fill"
	) + labs(y="proportion")+ theme_classic()
ggsave("figures/cdcpi2_anchor_cluster_prop.pdf")

## Viz 4c: cluster proportions by condition colitis 2
reordered.conditions <- factor(s$colitis2, levels = c("Colitis (aPD1)" ,"No-Colitis (aPD1)" , "Colitis (combo)" ,  "No-Colitis (combo)", "Control", "CD_inflamed","CD_uninflamed"))
data.frame(cluster =s$seurat_clusters, condition= reordered.conditions) %>%
	ggplot() + geom_bar(
		mapping = aes(x= condition, fill= cluster),
		position = "fill"
	) + labs(y="proportion")+ theme_classic()
ggsave("figures/cdcpi2_anchor_cluster_prop_colitis2_by_condition.pdf")

# save
saveRDS(s,"saved_objects/cdcpi2_anchor.rds")
##################################################
# 10 - Clustering check with reverse annotating
##################################################
## Get annotations from Luoma-Suo cells
require(stringr)
require(future.apply)
# Finding the correct mapping: turns out it is the "-1_1, -2_1, etc." endings
# that match the Suo Tcells
colnames(s) %>% str_sub(17) %>% table()
colnames(Tcell) %>% str_sub(17) %>% table()
colnames(s)[as.logical(cells)] %>% str_sub(17) %>% table()
Tcell$names %>% class()
# executing the addtion
x <- Tcell$names
names(x) <- str_glue("{names(x)}_1")
s$old.Tcell.annotations <- x
# graph showing clusters: onlyTcells (CPI)
p <- DimPlot(s,reduction="umap",group.by="old.Tcell.annotations")
ggsave("figures/apply_old_CPI_labels_to_cdcpi21550.pdf")
p <- DimPlot(s,reduction="umap",group.by="old.Tcell.annotations", cells=!is.na(s$old.Tcell.annotations))
ggsave("figures/apply_old_CPI_labels_to_cdcpi21550_noNA.pdf")

## Check for version 1
s4 <- readRDS("saved_objects/cdcpi_anchor.rds") 
s4$old.Tcell.annotations <- x
# graph showing clusters: onlyTcells (CPI)
p <- DimPlot(s4,reduction="umap",group.by="old.Tcell.annotations")
ggsave("figures/apply_old_CPI_labels_to_cdcpi1550.pdf")
p <- DimPlot(s4,reduction="umap",group.by="old.Tcell.annotations", cells=!is.na(s4$old.Tcell.annotations))
ggsave("figures/apply_old_CPI_labels_to_cdcpi1550_noNA.pdf")

#################################################################
# 10 - Clustering check with comparison to within-batch clusters
#################################################################
library(pheatmap)
s<- readRDS("saved_objects/cdcpi2_anchor.rds")
s1 <- readRDS("saved_objects/martinTcells15200.rds")
s2 <- readRDS("saved_objects/onlyTcells1550_labeled.rds")
# For the first batch (adding +10 for a smoother color transition
# from zero to non-zero counts for any given matrix entry).
tab <- table(paste("after", s$seurat_clusters[s$orig.ident=="CD_martin"]),
    paste("before", s1$seurat_clusters))
heatCD <- pheatmap(log10(tab+10), cluster_row=FALSE, cluster_col=FALSE,
    main="CD comparison", silent=TRUE)
# For the second batch.
tab <- table(paste("after", s$seurat_clusters[s$orig.ident=="CD3_Tcell"]),
    paste("before", s2$seurat_clusters))
heatCPI <- pheatmap(log10(tab+10), cluster_row=FALSE, cluster_col=FALSE,
    main="CPI comparison", silent=TRUE)

pdf("figures/cdcpi2_cluster_comparison.pdf")
gridExtra::grid.arrange(heatCD[[4]], heatCPI[[4]])
dev.off()

####################################
# 5 - Basic Vizualizations
####################################
# feature plots the CD4 and CD8 split
FeaturePlot(s, features= c('CD4','CD8A','CD8B','CD3E'), reduction="umap")
ggsave("figures/cdcpi2_CD4xCD8_featureplot.pdf", width= 12, height= 12, units= "in")
# violin plot of CD4 per cluster
p <- VlnPlot(s,features = "CD4", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 90),legend.position = "none", axis.title.x=element_blank())
ggsave("figures/cdcpi2_VlnPlot_CD4_by_cluster.pdf", width= 7, height= 7, units= "in")
# violin plot of CD8 per cluster
p <- VlnPlot(s,features = "CD8A", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 90),legend.position = "none", axis.title.x=element_blank())
ggsave("figures/cdcpi2_VlnPlot_CD8_by_cluster.pdf", width= 7, height= 7, units= "in")
# violin plot of CD3E per cluster
p <- VlnPlot(s,features = "CD3E", pt.size= 0, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 90),legend.position = "none", axis.title.x=element_blank())
ggsave("figures/cdcpi2_VlnPlot_CD3E_by_cluster.pdf", width= 7, height= 7, units= "in")
# distribution by sample
data.frame(cluster = s$seurat_clusters, patient= s$sample.id) %>%
	ggplot() + geom_bar(
		mapping = aes(x= patient, fill= cluster),
		position = "fill"
	) + labs(y="proportion")+ theme_classic() + coord_flip()
ggsave("figures/cdcpi2_cluster_prop_per_sample.pdf",width= 7, height= 7, units= "in")
# differential gene analysis
DE_heatmap(s,"figures/cdcpi2_DE_avgexp.pdf")
DE_heatmap(s4,"figures/cdcpi_DE_avgexp.pdf")

####################################
# DA analysis: CD inflamed vs CD non-inflamed
####################################
library(edgeR)
x <- subset(s, orig.ident == "CD_martin")
x$sample.id <- droplevels(x$sample.id)
x$CD_inf <- x$colitis == "CD_inflamed"
design <- function(y) model.matrix(~factor(CD_inf),y)
out <- DA_analysis(x, design, title = "cdcpi2_CD_inf_vs_uninf")
summary(out$y.ab$common.dispersion)
summary(out$fit.ab$var.prior)
summary(out$fit.ab$df.prior)
summary(decideTests(out$res))
topTags(out$res)

####################################
# DA analysis: CPI colitis vs CPI no-colitis
####################################
x <- subset(s, colitis == "CPI_colitis" | colitis == "CPI_no-colitis")
x$sample.id <- droplevels(x$sample.id)
x$CPI_colitis <- x$colitis == "CPI_colitis"
design <- function(y) model.matrix(~factor(CPI_colitis),y)
out <- DA_analysis(x, design, title = "cdcpi2_CPI_colitis_vs_healthy")
summary(out$y.ab$common.dispersion)
summary(out$fit.ab$var.prior)
summary(out$fit.ab$df.prior)
summary(decideTests(out$res))
topTags(out$res)

####################################
# DA analysis: CPI no-colitis vs healthy
####################################
x <- subset(s, colitis == "Control" | colitis == "CPI_no-colitis")
x$sample.id <- droplevels(x$sample.id)
x$CPI <- x$colitis == "CPI_no-colitis"
design <- function(y) model.matrix(~factor(CPI),y)
out <- DA_analysis(x, design, title = "cdcpi2_CPI-no-colitis_vs_healthy")
summary(out$y.ab$common.dispersion)
summary(out$fit.ab$var.prior)
summary(out$fit.ab$df.prior)
summary(decideTests(out$res))
topTags(out$res)

####################################
# DA analysis: CD inf vs healthy
####################################
x <- subset(s, colitis == "Control" | colitis == "CD_inflamed")
x$sample.id <- droplevels(x$sample.id)
x$CD_inflamed_over_Control <- x$colitis == "CD_inflamed"
design <- function(y) model.matrix(~factor(CD_inflamed_over_Control),y)
out <- DA_analysis(x, design, title = "cdcpi2_CD_inflamed_vs_healthy")
summary(out$y.ab$common.dispersion)
summary(out$fit.ab$var.prior)
summary(out$fit.ab$df.prior)
summary(decideTests(out$res))
topTags(out$res)

########################################
# DA analysis: CD samples vs CPI samples
########################################
x <- s
x$CDallvsCPIall <- x$orig.ident == "CD_martin"
design <- function(y) model.matrix(~factor(CDallvsCPIall),y)
out <- DA_analysis(x, design, title = "cdcpi2_martin_vs_luoma-suo")
summary(out$y.ab$common.dispersion)
summary(out$fit.ab$var.prior)
summary(out$fit.ab$df.prior)
summary(decideTests(out$res))
topTags(out$res)

x <- subset(s, colitis!= "Control")
x$sample.id <- droplevels(x$sample.id)
x$CDvsCPInocontrol <- x$orig.ident == "CD_martin"
design <- function(y) model.matrix(~factor(CDvsCPInocontrol),y)
out <- DA_analysis(x, design, title = "cdcpi2_martin_vs_luoma-suo_nocontrol")
summary(out$y.ab$common.dispersion)
summary(out$fit.ab$var.prior)
summary(out$fit.ab$df.prior)
summary(decideTests(out$res))
topTags(out$res)

########################################
# DA analysis: CD vs Control
########################################
x <- subset(s, colitis == "Control"| orig.ident=="CD_martin")
x$sample.id <- droplevels(x$sample.id)
x$sample.id <- x$patient.id
x$CDvsControl <- x$orig.ident == "CD_martin"
design <- function(y) model.matrix(~factor(CDvsControl),y)
out <- DA_analysis(x, design, title = "cdcpi2_CD_vs_Control")
summary(out$y.ab$common.dispersion)
summary(out$fit.ab$var.prior)
summary(out$fit.ab$df.prior)
summary(decideTests(out$res))
topTags(out$res)


########################################
# DE analysis: CD vs CPI: Treg comparison
########################################
library(SingleCellExperiment)
library(scater)
x <- s
x$CDallvsCPIall <- x$orig.ident == "CD_martin"
DefaultAssay(x) <- "RNA"
merged <- as.SingleCellExperiment(x)
summed <- aggregateAcrossCells(merged, 
    id=colData(merged)[,c("seurat_clusters", "sample.id")])
summed
library(edgeR)
label <- 7
current <- summed[, summed$seurat_clusters==label]
y <- DGEList(counts(current), samples = colData(current))
# remove label-sample combinations less than 20
discarded <- current$ncells < 20
y <- y[,!discarded]
summary(discarded)
# filter the genes
keep <- filterByExpr(y, group=current$CDallvsCPIall)
y <- y[keep,]
summary(keep)
# add normalizing factors
y <- calcNormFactors(y)
design <- model.matrix(~factor(CDallvsCPIall), y$samples)
design
y <- estimateDisp(y, design)
summary(y$trended.dispersion)
fit <- glmQLFit(y, design, robust=TRUE)
summary(fit$var.prior)
res <- glmQLFTest(fit, coef=ncol(design))
summary(decideTests(res))
topTags(res)
# figures
	pdf("figures/cdcpi2_DE_CDallvsCPIall_BVC.pdf")
	plotBCV(y)
	dev.off()
	pdf("figures/cdcpi2_DE_CDallvsCPIall_QLDisp.pdf")
	plotQLDisp(fit)
	dev.off()

library(EnhancedVolcano)

DE_table <- topTags(res, n=6074)
topTags(res, n=50) %>% rownames() -> labels
EnhancedVolcano(DE_table$table,
    lab = rownames(DE_table$table),
    x = 'logFC',
    y = 'FDR',
    xlim = c(-20, 20),
    title = 'CD versus CPI: Treg comparison',
    pCutoff = 0.05,
    FCcutoff = 5,
    pointSize = 3.0,
    labSize = 3.0,
 #   selectLab = labels,
    drawConnectors = F,
    widthConnectors = 0.2,
    colConnectors = 'grey30')
ggsave("figures/cdcpi2_DE_CDallvsCPIall_volcano.pdf", width= 10, height= 8, units= "in")

########################################
# DE analysis: CD vs CPI
########################################
library(scran)
x <- s
x$CDallvsCPIall <- x$orig.ident == "CD_martin"
DefaultAssay(x) <- "RNA"
merged <- as.SingleCellExperiment(x)
summed <- aggregateAcrossCells(merged, 
    id=colData(merged)[,c("seurat_clusters", "sample.id")])
summed.filt$sample.id <- as.factor(summed.filt$sample.id)
levels(summed.filt$sample.id) <- levels(merged$sample.id)
summed.filt <- summed[,summed$ncells >= 20]
# Pulling out a sample-level 'targets' data.frame:
targets <- colData(merged)[!duplicated(merged$sample.id),]
# Constructing the design matrix:
design <-  model.matrix(~factor(CDallvsCPIall), data=targets)
rownames(design) <- targets$sample.id
factor(summed.filt$sample.id , levels(merged$sample.id))
de.results <- pseudoBulkDGE(summed.filt, 
    sample=summed.filt$sample.id,
    label=summed.filt$seurat_clusters,
    design=design,
    coef=ncol(design),
    # 'condition' sets the group size for filterByExpr(),
    # to perfectly mimic our previous manual analysis.
    condition=targets$CDallvsCPIall 
)
is.de <- decideTestsPerLabel(de.results, threshold=0.05)
summarizeTestsPerLabel(is.de)
# Upregulated across most cell types.
up.de <- is.de > 0 & !is.na(is.de)
head(sort(rowMeans(up.de), decreasing=TRUE), 50)
# Downregulated across cell types.
down.de <- is.de < 0 & !is.na(is.de)
head(sort(rowMeans(down.de), decreasing=TRUE), 50)
#label specific DE
remotely.de <- decideTestsPerLabel(de.results, threshold=0.5)
not.de <- remotely.de==0 | is.na(remotely.de)

other.labels <- setdiff(colnames(not.de), "7")
unique.degs <- is.de[,"7"]!=0 & rowMeans(not.de[,other.labels])==1
unique.degs <- names(which(unique.degs))
unique.degs

# Choosing the top-ranked gene for inspection:
de.treg <- de.results$"7"
de.treg <- de.treg[order(de.treg$PValue),]
de.treg <- de.treg[rownames(de.treg) %in% unique.degs,]

sizeFactors(summed.filt) <- NULL
plotExpression(logNormCounts(summed.filt), 
    features=rownames(de.treg)[1],
    x="CDallvsCPIall", colour_by="CDallvsCPIall", 
    other_fields="seurat_clusters") + 
    facet_wrap(~seurat_clusters)
ggsave("figures/cdcpi2_DE_CDallvsCPIall_plotexpression_treg.pdf", width= 10, height= 10, units= "in")

# Choosing the top-ranked gene for volcano plot:
library(EnhancedVolcano)
de.treg <- de.results$"7"
de.treg <- de.treg[order(de.treg$PValue),]
EnhancedVolcano(de.treg,
    lab = rownames(de.treg),
    x = 'logFC',
    y = 'FDR',
    xlim = c(-20, 20),
    title = 'CD versus CPI: Treg comparison',
    pCutoff = 0.05,
    FCcutoff = 5,
    pointSize = 3.0,
    labSize = 3.0,
    selectLab = unique.degs,
    drawConnectors = T,
    widthConnectors = 0.2,
    colConnectors = 'grey30')
ggsave("figures/cdcpi2_DE_CDallvsCPIall_volcano_tregs_2.pdf", width= 10, height= 8, units= "in")
####################################
# 5 - Cluster Annotations
####################################
load('data/ssuo_CD3/seurat.object.RData')
Tcell <- process_ssuo_Tcells(Tcell)
martin <- read.csv(file="data/martin_categories.csv", header=TRUE,sep=",") %>% process_martin()

require(SingleR)
require(BiocParallel)
pred.CPI <- SingleR(method='cluster',
	test=as.matrix(s[['integrated']]@data), 
	clusters= s$seurat_clusters,
	ref=as.matrix(Tcell[['RNA']]@data), 
	labels=Tcell$names,
	de.method="wilcox",
	BPPARAM=MulticoreParam(10))
s$SingleR.Luoma <- s$seurat_clusters
levels(s$SingleR.Luoma) <- pred.CPI$labels

pred.Martin <- SingleR(method='cluster',  
	test=s[['integrated']]@data, 
	clusters= s$seurat_clusters,
	ref=martin, 
	labels=colnames(martin),
	BPPARAM=MulticoreParam(10))
s$SingleR.Martin <- s$seurat_clusters
levels(s$SingleR.Martin) <- pred.Martin$labels

ref <- MonacoImmuneData()
pred.Monaco <- SingleR(method="cluster",
	test=s[['integrated']]@data,
	clusters= s$seurat_clusters, 
	ref=ref, 
	labels=ref$label.fine, 
	BPPARAM=MulticoreParam(10))
s$SingleR.Monaco <- s$seurat_clusters
levels(s$SingleR.Monaco) <- pred.Monaco$labels

saveRDS(list(pred.CPI,pred.Martin,pred.Monaco),"saved_objects/label_info_for_cdcpi21550.rds")


require(ggpubr)
p1 <- DimPlot(s, reduction = "umap", label=T, group.by="SingleR.Martin") + labs(title="Martin reference")
p2 <- DimPlot(s, reduction = "umap", label=T, group.by="SingleR.Luoma") + labs(title="Luoma-Suo reference")
p3 <- DimPlot(s, reduction = "umap", label=T, group.by="SingleR.Monaco") + labs(title="Monaco reference")
ggarrange(p1,p2,p3, ncol=3, nrow=1)
ggsave("figures/cdcpi2_anchor_labeling.pdf", width= 24, height= 6, units= "in")

saveRDS(s, "saved_objects/cdcpi2_anchor.rds")

