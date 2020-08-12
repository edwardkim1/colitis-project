# 6/26-27/20
library(ggplot2)
library(Seurat)
setwd("~/1.projects/PD1")
p1294 <- Read10X(data.dir = "data/DFCI-1294-CD3/outs/filtered_feature_bc_matrix")

#filter based on low library size
filter_by_libsize <- function(seurat.object, save = FALSE, filename = "") {
	lib.sizes = Matrix::colSums(seurat.object)
	lls <- log(lib.sizes)
	lls.p = pnorm(lls, mean = median(lls), sd = mad(lls), lower.tail = TRUE)
	lls.lim = max(lls[which(p.adjust(lls.p, method = "fdr") < 1e-2)])

	filtered.object <- seurat.object[,which(lls > lls.lim)]
	new.object = list(filtered.counts.mat = filtered.object, pre.lls = lls, post.lls = lls[which(lls > lls.lim)], total.removed = sum(!(lls > lls.lim)) , median = median(lls), mad =mad(lls), lls.lim = lls.lim)
	if(save==TRUE) {
		saveRDS(new.object, filename)
	}

	return(new.object)
}
q1.p1294 <- filter_by_libsize(p1294)

p <- ggplot(data=data.frame(x= c(rep("pre lls", length(q1.p1294$pre.lls)),rep("post lls", length(q1.p1294$post.lls))), log.lib.size = c(q1.p1294$pre.lls,q1.p1294$post.lls))) + geom_violin(aes(x= x, y=log.lib.size))
ggsave("figures/loglibsizeviolin.pdf")


## filter based on number of features
# need to check for seurat.object compatibility; 
# currently works on sparse matrices 
filter_by_numfeats <- function(seurat.object, save = FALSE, filename = "") {
	n.feats = diff(seurat.object@p)
	lnf <- log(n.feats)
	lnf.p = pnorm(lnf, mean = median(lnf), sd = mad(lnf), lower.tail = TRUE)
	lnf.lim = max(lnf[which(p.adjust(lnf.p, method = "fdr") < 1e-2)])

	filtered.object <- seurat.object[,which(lnf > lnf.lim)]
	new.object = list(filtered.counts.mat = filtered.object, pre.lnf = lnf, post.lnf = lnf[which(lnf > lnf.lim)], total.removed = sum(!(lnf > lnf.lim)) , median = median(lnf), mad =mad(lnf), lnf.lim = lnf.lim)
	if(save==TRUE) {
		saveRDS(new.object, filename)
	}

	return(new.object)
}
q1.p1294 <- filter_by_numfeats(p1294)
q1.p1294$lnf.lim
q1.p1294$total.removed
q1.p1294$median - 3 * q1.p1294$mad 

p <- ggplot(data=data.frame(x= c(rep("pre lnf", length(q1.p1294$pre.lnf)),rep("post lnf", length(q1.p1294$post.lnf))), log.num.feats = c(q1.p1294$pre.lnf,q1.p1294$post.lnf))) + geom_violin(aes(x= x, y=log.num.feats))
ggsave("figures/lognumfeatsviolin.pdf")

#### Filter based on mitochondrial
filter_by_mito <- function(seurat.object, save = FALSE, filename = "") {
	genename=as.matrix(rownames(seurat.object))
	rownames(genename)=rownames(seurat.object)
	colnames(genename)='Symbol'
	lib.sizes = Matrix::colSums(seurat.object)
	mito.gene=rownames(genename)[grepl(x=rownames(genename), pattern = "^MT-")]
	mt.counts = seurat.object[mito.gene,]
	mt.fraction = Matrix::colSums(mt.counts)/lib.sizes
	mt.p = pnorm(mt.fraction, mean = median(mt.fraction), sd = mad(mt.fraction), lower.tail = FALSE)
	mt.lim = min(mt.fraction[which(p.adjust(mt.p, method = "fdr") < 1e-2)])
	filtered.object <- seurat.object[,which(mt.fraction < mt.lim)]
	new.object = list(filtered.counts.mat = filtered.object, pre.mt.fraction = mt.fraction, post.mt.fraction= mt.fraction[which(mt.fraction < mt.lim)], total.removed = sum(!(mt.fraction < mt.lim)) , median = median(mt.fraction), mad =mad(mt.fraction), mt.lim = mt.lim)
	if(save==TRUE) {
		saveRDS(new.object, filename)
	}

	return(new.object)
}
q1.p1294 <-filter_by_mito(p1294)
q1.p1294 <-filter_by_mito(p1294, save=T, filename= "saved_objects/p1294_mito_filtering.RDS")
### Figures to evaluate how good the mito filtering is performing
p <- ggplot(data=data.frame(x= c(rep("pre mito%", length(q1.p1294$pre.mt.fraction)),rep("post mito%", length(q1.p1294$post.mt.fraction))), mt.fraction = c(q1.p1294$pre.mt.fraction,q1.p1294$post.mt.fraction))) + geom_violin(aes(x= x, y=mt.fraction))
ggsave("figures/mitohist.pdf")

source("scripts/Stat111functions.R")
mito.data <- q1.p1294$pre.mt.fraction
p <- plot_kde(mito.data, kernel = "gaussian", bw=0.01, lab.x = "mitochondrial gene fraction (per cell)", overlay = T, model="Norm", 
	dmodel = function(x) {dnorm(x,median(mito.data),mad(mito.data))}, 
	support=seq(0,1,0.001))
p + geom_vline(xintercept= q1.p1294$mt.lim, color = "red") +
	geom_vline(xintercept= median(mito.data)+ 3*mad(mito.data), linetype = "dotted")
ggsave("figures/mito-kde-Normal.pdf")

q1.p1294$total.removed/length(q1.p1294$pre.mt.fraction) #to find percentage filtered

### Proceeding with cells only filtered by the mito filtering scheme
## scaling and log-transformation
post.qc <- readRDS("saved_objects/p1294_mito_filtering.RDS")
counts.mat <- post.qc$filtered.counts.mat
lib.sizes <- Matrix::colSums(counts.mat)
lib.sizes
head(counts.mat/lib.sizes,50)

library(future)
plan("multiprocess", workers = 10)

s1 <- CreateSeuratObject(counts = counts.mat, project = "p1294_CD3", min.cells = 0, min.features = 0)
s1 <- NormalizeData(s1, normalization.method = "LogNormalize", scale.factor = 10000)
s1 <- FindVariableFeatures(s1, selection.method = "vst", nfeatures = 3000) ## Highly variable genes
hvf <- VariableFeatures(s1)
'%ni%' = Negate('%in%')
all.genes <- rownames(s1)
trvgenes <- all.genes[grepl(x=all.genes, pattern = "^TR.V")]
igvgenes <- all.genes[grepl(x=all.genes, pattern = "^IG.V")]
# Normalized
VariableFeatures(s1) <- hvf[hvf %ni% c(trvgenes,igvgenes)]
# Center expression for each feature
s1 <- ScaleData(s1, features = all.genes)
s1 <- RunPCA(s1, features = VariableFeatures(s1)) ## PCA
pdf("figures/p1294_CD3_3000hvf.pdf", width = 10, height = 10)
ElbowPlot(s1, ndims=30)
dev.off()
s1 <- FindNeighbors(s1, dims = 1:20) ## KNN + Louvian algorithm ##
s1 <- FindClusters(s1, resolution = 0.5)
s1 <- RunUMAP(s1, dims = 1:20) ## UMAP ##

pdf("figures/p1294_CD3_post_mito_umap.pdf", width = 10, height = 10)
DimPlot(s1, reduction = "umap", label=T)
dev.off()

dirname= "p1294_CD3"
nPCs <- 25
## Detecting Doublets: Doublet Finder ##
## pK identification
library(DoubletFinder)
sweep.res.list <- paramSweep_v3(s1, PCs = 1:nPCs)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
pdf("figures/p1294_CD3_DF_pK.pdf")
s1_pK <- find.pK(sweep.stats)
dev.off()
pK_max <- s1_pK[which.max(s1_pK$BCmetric),]$pK
pK_max <- as.numeric(levels(pK_max)[as.integer(pK_max)])
## Output pK maximization using ggplot2 as pdf
s1_pK$group <- numeric(length(s1_pK$ParamID))
p <- ggplot(data=s1_pK, aes(x=pK, y=BCmetric, group=0)) +
  		geom_line()+ scale_x_discrete(breaks=seq(0,0.3,0.05)) +
  		geom_point() + ggtitle(dirname) + theme_classic() 
ggsave("figures/p1294_CD3_DF_pK_ggplot.pdf")
## Calculation
homotypic.prop <- modelHomotypic(s1@meta.data$seurat_clusters)   ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.1*nrow(s1@meta.data))   ## Assuming 10% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

s1 <- doubletFinder_v3(s1, PCs = 1:nPCs, pN = 0.25, pK= pK_max, nExp = nExp_poi, reuse.pANN = FALSE)
first_pANN <- names(s1@meta.data)[grep("pANN", names(s1@meta.data))]
s1 <- doubletFinder_v3(s1, PCs = 1:nPCs, pN = 0.25, pK = pK_max, nExp = nExp_poi.adj, reuse.pANN = first_pANN)
DFs <- names(s1@meta.data)[grep("DF.class", names(s1@meta.data))]

s1@meta.data[,"DF_hi.lo"] <- s1[[DFs[1]]]
s1@meta.data$DF_hi.lo[which(s1@meta.data$DF_hi.lo == "Doublet" & s1[[DFs[2]]] == "Singlet")] <- "Doublet_lo"
s1@meta.data$DF_hi.lo[which(s1@meta.data$DF_hi.lo == "Doublet")] <- "Doublet_hi"
# Output predoublet removal
saveRDS(s1,"saved_objects/p1294_CD3_labeled_doublets.RDS")
pdf("figures/p1294_CD3_DF_umap.pdf")
DimPlot(s1, reduction = "umap", group.by="DF_hi.lo")+ scale_colour_manual(values=c("red","yellow","gray"))
dev.off()

table(s1$DF_hi.lo)

s1 <- subset(s1, subset = DF_hi.lo != "Doublet_hi" ) #eliminating hi probability doublets
# Output post doublet removal
saveRDS(s1, paste(dirname,"postQC.rds", sep="_"))

#### Investigating the properties of Shengbao's T cell Seurat object
load("../CD3/colitis/seurat.object.RData")
Tcell
dim(Tcell)
tcell.qc <- filter_by_mito(Tcell)
tcell.qc$total.removed/length(tcell.qc$pre.mt.fraction)
tcell.qc

p <- ggplot(data=data.frame(x= c(rep("pre mito%", length(tcell.qc$pre.mt.fraction)),rep("post mito%", length(tcell.qc$post.mt.fraction))), mt.fraction = c(tcell.qc$pre.mt.fraction,tcell.qc$post.mt.fraction))) + geom_violin(aes(x= x, y=mt.fraction))
ggsave("figures/ShengbaoTcell_mitoviolin.pdf")

#separate analysis
tcell.qc <- filter_by_numfeats(Tcell[['RNA']]@counts)
tcell.qc$total.removed/length(tcell.qc$pre.lnf) # 0.00569
tcell.qc

diff(Tcell[['RNA']]@counts@p)
Tcell$nFeature_RNA # nFeature_RNA provides the same information!!

p <- ggplot(data=data.frame(x= c(rep("pre lnf", length(tcell.qc$pre.lnf)),rep("post lnf", length(tcell.qc$post.lnf))), log.num.feats = c(tcell.qc$pre.lnf,tcell.qc$post.lnf))) + geom_violin(aes(x= x, y=log.num.feats))
ggsave("figures/ShengbaoTcell_lognumfeatsviolin.pdf")


####### Completely fresh start at filtering
filter_stats <- function(seurat.object, save = FALSE, filename = "") {
	#thresholding number of unique features
	lnf <- log(seurat.object$nFeature_RNA)
	lnf.p = pnorm(lnf, mean = median(lnf), sd = mad(lnf), lower.tail = TRUE)
	lnf.lim = max(lnf[which(p.adjust(lnf.p, method = "fdr") < 1e-2)])
	#thresholding by mitochondrial fraction
	mt.fraction <- seurat.object$percent.mt
	mt.p = pnorm(mt.fraction, mean = median(mt.fraction), sd = mad(mt.fraction), lower.tail = FALSE)
	mt.lim = min(mt.fraction[which(p.adjust(mt.p, method = "fdr") < 1e-2)])
	# output
	filter.stats = list(mt.pre = mt.fraction, mt.post = mt.fraction[which(mt.fraction < mt.lim)], lnf.pre = lnf, lnf.post = lnf[which(lnf > lnf.lim)], cells.to.remove= !(lnf > lnf.lim & mt.fraction < mt.lim), mt.remove = sum(!(mt.fraction < mt.lim)) , mt.median = median(mt.fraction), mt.mad =mad(mt.fraction), mt.lim = mt.lim, lnf.remove = sum(!(lnf > lnf.lim)) , lnf.median = median(lnf), lnf.mad =mad(lnf), lnf.lim = lnf.lim, total.remove = sum(!(lnf > lnf.lim & mt.fraction < mt.lim)))
	if(save==TRUE) {
		saveRDS(filter.stats, filename)
	}

	return(filter.stats)
}

s <- CreateSeuratObject(counts = p1294, project = "p1294_CD3", min.cells = 0, min.features = 0) %>% PercentageFeatureSet(pattern = "^MT-", col.name = "percent.mt")
stats <- filter_stats(s)
stats <- filter_stats(s, save=T, filename="saved_objects/p1294_CD3_filterstats.RDS")
stats$lnf.remove
stats$mt.remove
stats$total.remove
stats$cells.to.remove
length(stats$cells.to.remove)
which(!(lnf > lnf.lim & mt.fraction < mt.lim))
s$cells.to.remove <- stats$cells.to.remove

pdf("figures/p1294_CD3_featurescatter.pdf")
FeatureScatter(s, feature1="nFeature_RNA", feature2="percent.mt", group.by="cells.to.remove")
dev.off()

#remove
s$cells.to.remove <- stats$cells.to.remove
s <- s[,!s$cells.to.remove]
sum(s$cells.to.remove) # == 0
# normalize
s <- NormalizeData(s, normalization.method = "LogNormalize", scale.factor = 10000) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 3000) ## Highly variable genes
hvf <- VariableFeatures(s)
'%ni%' = Negate('%in%')
all.genes <- rownames(s)
trvgenes <- all.genes[grepl(x=all.genes, pattern = "^TR.V")]
igvgenes <- all.genes[grepl(x=all.genes, pattern = "^IG.V")]
VariableFeatures(s) <- hvf[hvf %ni% c(trvgenes,igvgenes)]
# Center expression for each feature
s <- ScaleData(s, features = all.genes, vars.to.regress="percent.mt") %>% RunPCA(features=VariableFeatures(s)) %>% RunUMAP(dims = 1:20) %>% FindNeighbors(dims = 1:20) %>% FindClusters(resolution = 0.5) 
s <- remove_doublets(s, nPCs=20)

DimPlot(s, reduction = "umap", group.by="DF_hi.lo")+ scale_colour_manual(values=c("red","yellow","gray"))
ggsave("figures/p1294_CD3_DF_umap_2.pdf")

saveRDS(s, "saved_objects/p1294_CD3_labeled_doublets.RDS")

#### sctransform
library(sctransform)
s <- SCTransform(s, vars.to.regress = "percent.mt") 
hvf <- VariableFeatures(s)
'%ni%' = Negate('%in%')
all.genes <- rownames(s)
trvgenes <- all.genes[grepl(x=all.genes, pattern = "^TR.V")]
igvgenes <- all.genes[grepl(x=all.genes, pattern = "^IG.V")]
VariableFeatures(s) <- hvf[hvf %ni% c(trvgenes,igvgenes)]
s <- RunPCA(s) %>% FindNeighbors(dims = 1:nPCs) %>% RunUMAP(dims = 1:nPCs) %>% FindNeighbors(dims = 1:nPCs) %>% FindClusters(resolution= 0.5)

DimPlot(s)
ggsave("figures/p1294_CD3_SCT.pdf")


### Previous log-normalization method
s <- NormalizeData(s, normalization.method = "LogNormalize", scale.factor = 10000)
mean(colSums(s[['RNA']]@data)) #cannot seem to interpret



