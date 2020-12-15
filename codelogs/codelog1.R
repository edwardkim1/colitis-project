https://bioconductor.org/packages/release/bioc/vignettes/DropletUtils/inst/doc/DropletUtils.html


library(DropletUtils)
sce <- read10xCounts("../data/martin/GSM3972009_69")


Users can also load multiple samples at once by supplying a character vector to read10xCounts. This will return a single SingleCellExperiment where all of the individual matrices are combined by column. Obviously, this only makes sense when the same set of genes is being used across samples.

my.counts <- counts(sce)
br.out <- barcodeRanks(my.counts)

plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total")
o <- order(br.out$rank)
lines(br.out$rank[o], br.out$fitted[o], col="red")

abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"), 
    legend=c("knee", "inflection"))


set.seed(100)
e.out <- emptyDrops(my.counts)
e.out

is.cell <- e.out$FDR <= 0.01
sum(is.cell, na.rm=TRUE)

# Diagnostics

table(Limited=e.out$Limited, Significant=is.cell)

plot(e.out$Total, -e.out$LogProb, col=ifelse(is.cell, "red", "black"),
    xlab="Total UMI count", ylab="-Log Probability")


# swapped drops? <- applicable to multiplex sequencing


#### Now load R with Seurat ####
library(Seurat)
e.out <- readRDS("e.out.GSM3972009_69")
s1.data <- Read10X(data.dir = "../data/martin/GSM3972009_69/")

## Quality Control (remove emptydroplets)
is.cell <- e.out$FDR <= 0.01
is.cell <- ifelse(is.na(is.cell),F,is.cell)
s1.qc1 <- s1.data[,is.cell]

## Quality Control (mito hypothesis testing)
genename=as.matrix(rownames(s1.qc1))
rownames(genename)=rownames(s1.qc1)
colnames(genename)='Symbol'
lib.sizes = Matrix::colSums(s1.qc1)
mito.gene=rownames(genename)[grepl(x=rownames(genename), pattern = "^MT-")]
mt.counts = s1.qc1[mito.gene,]
mt.fraction = Matrix::colSums(mt.counts)/lib.sizes
mt.p = pnorm(mt.fraction, mean = median(mt.fraction), sd = mad(mt.fraction), lower.tail = FALSE)
mt.lim = min(mt.fraction[which(p.adjust(mt.p, method = "fdr") < 1e-2)])
s1.qc2 <- s1.qc1[,which(mt.fraction < mt.lim)]

## Some qc stats:
sum(mt.counts) #verify that counts are present

## Create Seurat Object ##
s1 <- CreateSeuratObject(counts = s1.qc2, project = "GSM3972009_69", min.cells = 0, min.features = 0)
s1[["percent.mt"]] <- PercentageFeatureSet(s1, pattern = "^MT-")
VlnPlot(s1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

s1 <- NormalizeData(s1, normalization.method = "LogNormalize", scale.factor = 10000)
## Highly variable genes
s1 <- FindVariableFeatures(s1, selection.method = "vst", nfeatures = 3000)
top10 <- head(VariableFeatures(s1), 10)
plot1 <- VariableFeaturePlot(s1)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2
## Scale to mean 0 var 1
all.genes <- rownames(s1)
s1 <- ScaleData(s1, features = all.genes)
## PCA
s1 <- RunPCA(s1, features = VariableFeatures(object = s1))
## KNN + Louvian algorithm ##
s1 <- FindNeighbors(s1, dims = 1:25)
s1 <- FindClusters(s1, resolution = 0.4)
head(Idents(s1), 5)
## UMAP ##
s1 <- RunUMAP(s1, dims = 1:25)
DimPlot(s1, reduction = "umap")
s1 <- RunTSNE(s1, dims = 1:25)
DimPlot(s1, reduction = "tsne")
## Detecting Doublets: Doublet Finder ##
## pK identification
library(DoubletFinder)
sweep.res.list <- paramSweep_v3(s1, PCs = 1:25)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
s1_pK <- find.pK(sweep.stats)
pK_max <- s1_pK[which.max(s1_pK$BCmetric),]$pK
pK_max <- as.numeric(levels(pK_max)[as.integer(pK_max)])

## Visualize pK using ggplot2
s1_pK$group <- numeric(length(s1_pK$ParamID))
ggplot(data=s1_pK, aes(x=pK, y=BCmetric, group=0)) +
  geom_line()+ scale_x_discrete(breaks=seq(0,0.3,0.05)) +
  geom_point() + theme_classic()


homotypic.prop <- modelHomotypic(s1@meta.data$seurat_clusters)   ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.1*nrow(s1@meta.data))   ## Assuming 10% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

s1 <- doubletFinder_v3(s1, PCs = 1:25, pN = 0.25, pK= pK_max, nExp = nExp_poi, reuse.pANN = FALSE)
s1 <- doubletFinder_v3(s1, PCs = 1:25, pN = 0.25, pK = pK_max, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.04_361")

s1@meta.data[,"DF_hi.lo"] <- s1@meta.data$DF.classifications_0.25_0.04_361
s1@meta.data$DF_hi.lo[which(s1@meta.data$DF_hi.lo == "Doublet" & s1@meta.data$DF.classifications_0.25_0.04_309 == "Singlet")] <- "Doublet_lo"
s1@meta.data$DF_hi.lo[which(s1@meta.data$DF_hi.lo == "Doublet")] <- "Doublet_hi"
DimPlot(s1, reduction = "umap", group.by="DF_hi.lo", plot.order=c("Singlet","Doublet_lo","Doublet_hi")) + scale_colour_manual(values=c("red","yellow","gray"))

table(s1$DF_hi.lo)
s1 <- subset(s1, subset = DF_hi.lo != "Doublet_hi" ) #eliminating hi probability doublets

s1.markers <- FindAllMarkers(s1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
library(dplyr)
top10 <- s1.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

# Save the heatmap
pdf("s1_heatmap.pdf", width = 14, height = 10)
DoHeatmap(s1, features = top10$gene) + NoLegend()
dev.off()
scp -r esk17@rcapps5.dfci.harvard.edu:/mnt/beegfs/home/esk17/live/*.pdf ~/Desktop 


#### Next steps: Integration (in the qc_files directory)
dirname <- dir("../",pattern="GSM")
cells <- mclapply(paste(dirname,"postQC.rds", sep="_"),readRDS)

## Visualizing the valid cells bargraph
cellbg <- numeric(22)
for(i in 1:22) {
	cellbg[i] <- ncol(cells[[i]])
}

library(ggplot2)
df.bg <- data.frame(patient = paste("rp",rep(c("05","06","07","08",10,11,12,13,14,15,16),each=2)),
	number = cellbg[c(2,1,4,3,6,5,7,8,10,9,11,12,13,14,15,16,17,18,19,20,21,22)],
	categories = rep(c("uninflammed","inflammed"),11))
ggplot(df.bg, aes(patient, number, fill = categories)) + 
  geom_bar(stat="identity",position="dodge") + 
  scale_fill_brewer(palette = "Set1")


x <- c(2,1,4,3,6,5,7,8,10,9,11,12,13,14,15,16,17,18,19,20,21,22)
names(x) <- df.bg$categories
x[names(x)=="uninflammed"]


#### Generating the list below excluding patient 5 and 16
input <- "cells[[4]]"
for(i in 5:20) {
	input <- paste(input,", cells[[",sep="")
	input <- paste(input,i,sep="")
	input <- paste(input,"]]",sep="")
}


s2.wo5 <- merge(cells[[3]], c(cells[[4]], cells[[5]], cells[[6]], cells[[7]], cells[[8]], cells[[9]], cells[[10]], cells[[11]], cells[[12]], cells[[13]], cells[[14]], cells[[15]], cells[[16]], cells[[17]], cells[[18]], cells[[19]], cells[[20]]))
#FeatureScatter(s2.wo5, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
s2 <- s2.wo5[['RNA']]@counts
# Filter based on mitochondrial reads
genename=as.matrix(rownames(s2))
rownames(genename)=rownames(s2)
colnames(genename)='Symbol'
lib.sizes = Matrix::colSums(s2)
mito.gene=rownames(genename)[grepl(x=rownames(genename), pattern = "^MT-")]
mt.counts = s2[mito.gene,]
mt.fraction = Matrix::colSums(mt.counts)/lib.sizes
mt.p = pnorm(mt.fraction, mean = median(mt.fraction), sd = mad(mt.fraction), lower.tail = FALSE)
mt.lim = min(mt.fraction[which(p.adjust(mt.p, method = "fdr") < 1e-2)])
s2.qc1 <- s2[,which(mt.fraction < mt.lim)]
#s2[["percent.mt"]] <- PercentageFeatureSet(s2, pattern = "^MT-")
#VlnPlot(s2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
saveRDS(c(median(mt.fraction),mad(mt.fraction),mt.lim), "s2_wo5_mito_stats.rds")
# Filter the genes the appear infrequently
s2.qc2 <- CreateSeuratObject(s2.qc1,min.cells=10)
s2.qc2$orig.ident <- s2.wo5$orig.ident
# Filter cells that have low feature count
s2_postQC <- subset(s2.qc2, subset = nFeature_RNA > 200)

### filtering Genes from count matrix (WIP: doesn't work!)
#filterGenes <- function(object, min.cells) {
#	genes.use <- rownames(object[['RNA']]@counts)
#	if (min.cells > 0) {
#		num.cells <- Matrix::rowSums(object[['RNA']]@counts > 0)
#		genes.use <- names(num.cells[which(num.cells >= min.cells)])
#		object[['RNA']]@counts <- object[['RNA']]@counts[genes.use,]
#	}
#}
# Define Functions, Lists, and Operators
'%ni%' = Negate('%in%')
all.genes <- rownames(s2_postQC)
trvgenes <- all.genes[grepl(x=all.genes, pattern = "^TR.V")]
igvgenes <- all.genes[grepl(x=all.genes, pattern = "^IG.V")]
# Normalized
s2_postQC <- NormalizeData(s2_postQC, normalization.method = "LogNormalize", scale.factor = 10000)
## Highly variable genes
s2_postQC <- FindVariableFeatures(s2_postQC, selection.method = "vst", nfeatures = 2000)
hvf <- VariableFeatures(s2_postQC)
VariableFeatures(s2_postQC) <- hvf[hvf %ni% c(trvgenes,igvgenes)] 
top10 <- head(VariableFeatures(s2_postQC), 10)
plot1 <- VariableFeaturePlot(s2_postQC)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2
## Scale to mean 0 var 1
s2_postQC <- ScaleData(s2_postQC, features = all.genes)
## PCA
s2_postQC <- RunPCA(s2_postQC, features = VariableFeatures(s2_postQC))
## KNN + Louvian algorithm ##
s2_postQC <- FindNeighbors(s2_postQC, dims = 1:25)
s2_postQC <- FindClusters(s2_postQC, resolution = 0.5)
head(Idents(s2_postQC), 5)
## UMAP ##
s2_postQC <- RunUMAP(s2_postQC, dims = 1:25)
DimPlot(s2_postQC, reduction = "umap", label=T)
## Differential Gene Analysis with parallelization ##
library(future)
plan("multiprocess", workers = 40)
s2_postQC.markers <- FindAllMarkers(s2_postQC, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose=F)
top10 <- s2_postQC.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

#### Subset for T cells  ####
s3 <- subset(s2_postQC,idents= c(0,1,3,7,8,10,19))
s3[["percent.cd3"]] <- PercentageFeatureSet(s3, pattern = "^CD3")
s3

## Investigation
s3[['RNA']]@counts[c('CD3E','CD3G','CD3D','CD4','CD8A','CD8B','CD28','ICOS','CD2','CD7'),]
cd.counts <- Matrix::colSums(s3[['RNA']]@counts[c('CD3E','CD3G','CD3D','CD4','CD8A','CD8B','CD28','ICOS','CD2','CD7'),])
s3[['RNA']]@counts[c('CD3E'),]
s3[['RNA']]@counts[c('CD3G'),]
s3[['RNA']]@counts[c('CD3D'),]


#### s3 renormalization 1: without additional filtering
s3 <- NormalizeData(s3, normalization.method = "LogNormalize", scale.factor = 10000)
## Highly variable genes
s3 <- FindVariableFeatures(s3, selection.method = "vst", nfeatures = 2000)
hvf <- VariableFeatures(s3)
VariableFeatures(s3) <- hvf[hvf %ni% c(trvgenes,igvgenes)] 
top10 <- head(VariableFeatures(s3), 10)
plot1 <- VariableFeaturePlot(s3)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2
## Scale to mean 0 var 1
all.genes <- rownames(s3)
s3 <- ScaleData(s3, features = all.genes)
## PCA
s3 <- RunPCA(s3, features = VariableFeatures(s3))
## KNN + Louvian algorithm ##
s3 <- FindNeighbors(s3, dims = 1:10)
s3 <- FindClusters(s3, resolution = 0.5)
head(Idents(s3), 5)
## UMAP ##
s3 <- RunUMAP(s3, dims = 1:10)
DimPlot(s3, reduction = "umap", label=T)

## s3 renormalization 2: with addititional filtering
s3_2 <- s3[,cd.counts>0]
s3_2 <- NormalizeData(s3_2, normalization.method = "LogNormalize", scale.factor = 10000)
## Highly variable genes
s3_2 <- FindVariableFeatures(s3_2, selection.method = "vst", nfeatures = 2000)
hvf <- VariableFeatures(s3_2)
VariableFeatures(s3_2) <- hvf[hvf %ni% c(trvgenes,igvgenes)] 
top10 <- head(VariableFeatures(s3_2), 10)
plot1 <- VariableFeaturePlot(s3_2)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2
## Scale to mean 0 var 1
all.genes <- rownames(s3_2)
s3_2 <- ScaleData(s3_2, features = all.genes)
## PCA
s3_2 <- RunPCA(s3_2, features = VariableFeatures(s3_2))
## KNN + Louvian algorithm ##
s3_2 <- FindNeighbors(s3_2, dims = 1:10)
s3_2 <- FindClusters(s3_2, resolution = 0.5)
head(Idents(s3_2), 5)
## UMAP ##
s3_2 <- RunUMAP(s3_2, dims = 1:10)
DimPlot(s3_2, reduction = "umap", label=T)


#### Subset T cells (take 2) ####
s3 <- subset(s2_postQC,idents= c(0,1,3,10,17,19))
cd.counts <- Matrix::colSums(s3[['RNA']]@counts[c('CD3E','CD3G','CD3D','CD4','CD8A','CD8B','CD28','ICOS','CD2','CD7'),])
s3_2 <- s3[,cd.counts>0]
s3_2 <- NormalizeData(s3_2, normalization.method = "LogNormalize", scale.factor = 10000)
## Highly variable genes
s3_2 <- FindVariableFeatures(s3_2, selection.method = "vst", nfeatures = 2000)
hvf <- VariableFeatures(s3_2)
'%ni%' = Negate('%in%')
all.genes <- rownames(s3_2)
trvgenes <- all.genes[grepl(x=all.genes, pattern = "^TR.V")]
igvgenes <- all.genes[grepl(x=all.genes, pattern = "^IG.V")]
# Normalized
VariableFeatures(s3_2) <- hvf[hvf %ni% c(trvgenes,igvgenes)] 
top10 <- head(VariableFeatures(s3_2), 10)
# Adding patient ID
p.id <- ifelse(grepl(colnames(s3_2),pattern="_1$|_2$"),"Patient 6","")
p.id <- ifelse(grepl(colnames(s3_2),pattern="_3$|_4$"),"Patient 7",p.id)
p.id <- ifelse(grepl(colnames(s3_2),pattern="_5$|_6$"),"Patient 8",p.id)
p.id <- ifelse(grepl(colnames(s3_2),pattern="_7$|_8$"),"Patient 10",p.id)
p.id <- ifelse(grepl(colnames(s3_2),pattern="_9$|_10$"),"Patient 11",p.id)
p.id <- ifelse(grepl(colnames(s3_2),pattern="_11$|_12$"),"Patient 12",p.id)
p.id <- ifelse(grepl(colnames(s3_2),pattern="_13$|_14$"),"Patient 13",p.id)
p.id <- ifelse(grepl(colnames(s3_2),pattern="_15$|_16$"),"Patient 14",p.id)
p.id <- ifelse(grepl(colnames(s3_2),pattern="_17$|_18$"),"Patient 15",p.id)
s3_2$patient.id <- p.id
s3_2[["percent.mt"]] <- PercentageFeatureSet(s3_2, pattern = "^MT-")
## Scale to mean 0 var 1 and regress out patient variance
s3_2 <- ScaleData(s3_2, features = all.genes, vars.to.regress=c("patient.id","percent.mt"))
s3_2 <- RunPCA(s3_2, features = VariableFeatures(s3_2))
npc <- 30 # 20 PCS
s3_2 <- FindNeighbors(s3_2, dims = 1:npc, k.param = 50)
s3_2 <- FindClusters(s3_2, resolution = 0.5, graph.name ="RNA_snn")
## UMAP ##
s3_2 <- RunUMAP(s3_2, min.dist =0.30, graph="RNA_snn")
DimPlot(s3_2, reduction = "umap", label=T)

### Label Inflammed and Uninflammed ###
inflammed <-c("_1$","_3$","_6$","_7$","_10$","_12$","_14$","_16$","_18$")
toMatch <- paste(inflammed,collapse="|")
s3_2$status <- ifelse(grepl(colnames(s3_2), pattern=toMatch),"inflamed","not inflamed")

## Differential Gene Analysis with parallelization ##
library(future)
plan("multiprocess", workers = 40)
s3_2_regress.markers <- FindAllMarkers(s3_2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose=F)
top10 <- s3_2_regress.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

#### Single R stuff ####
## Prepare reference matrix
martin <- read.csv(file="martin_categories.csv", header=TRUE,sep=",")
df.martin <- as.data.frame(martin)
filtered.martin <- df.martin[which(df.martin$X != "28-Feb" & df.martin$X != "1-Mar"),] # need to remove the rows with "28-Feb" "1-Mar"
rownames(filtered.martin) <- filtered.martin$X
martin <- as.matrix(filtered.martin)
martin <- martin[, -1]
class(martin) <- "numeric"

## martin reference
sr.out <- SingleR(method='single',
		sc_data=s3_2[['RNA']]@counts, 
		ref_data=martin, 
		types=colnames(martin),
		numCores=40)
s3_2$SingleR <- sr.out$labels

pdf("s3v2_singler_labels.pdf", width = 10, height = 10)
DimPlot(s3_2, reduction = "umap", label=F, group.by="SingleR")
dev.off()

# @data, cluster
sr.out <- SingleR(method='cluster',
		sc_data=s3_2[['RNA']]@data, 
		ref_data=martin, 
		types=colnames(martin),
		numCores=40)

sr.out.s2.cluster <- SingleR(method='cluster',
		sc_data=s2_postQC[['RNA']]@counts, 
		clusters=s2_postQC$seurat_clusters,
		ref_data=martin, 
		types=colnames(martin),
		numCores=40)
levels(s2_postQC$seurat_clusters) <- sr.out.s2.cluster$labels

sr.out.cluster <- SingleR(method='cluster',
		sc_data=s3_2[['RNA']]@counts, 
		clusters= s3_2$seurat_clusters,
		ref_data=martin, 
		types=colnames(martin),
		numCores=40)
levels(s3_2$seurat_clusters) <- sr.out.cluster$labels

# SingleR with Luoma reference
load('../../CD3/colitis/seurat.object.RData')
load('../../CD3/colitis/seurat.marker.list.TRV.RData')

# scale data reference
sr.out.Luoma.cluster <- SingleR(method='cluster',
		sc_data=s3_2[['RNA']]@counts, 
		clusters= s3_2$seurat_clusters,
		ref_data=Tcell[['RNA']]@scale.data, 
		types=Tcell$seurat_clusters,
		numCores=40)
levels(s3_4$seurat_clusters) <- sr.out.Luoma.cluster$labels

sr.out.Luoma.cluster.2 <- SingleR(method='cluster',
		sc_data=s3_2[['RNA']]@counts, 
		clusters= s3_2$seurat_clusters,
		ref_data=as.matrix(Tcell[['RNA']]@counts), 
		types=Tcell$seurat_clusters,
		numCores=40)
levels(s3_3$seurat_clusters) <- sr.out.Luoma.cluster.2$labels

# trying with @data for both input and reference
# all cells
sr.out.s2.cluster <- SingleR(method='cluster',
		sc_data=s2_postQC[['RNA']]@data, 
		clusters=s2_postQC$seurat_clusters,
		ref_data=martin, 
		types=colnames(martin),
		numCores=40)
s2_postQC$SingleR <- s2_postQC$seurat_clusters
levels(s2_postQC$SingleR) <- sr.out.s2.cluster$labels
# T cells only
sr.out.Luoma.cluster.2 <- SingleR(method='cluster',
		sc_data=s3_2[['RNA']]@data, 
		clusters= s3_2$seurat_clusters,
		ref_data=as.matrix(Tcell[['RNA']]@data), 
		types=Tcell$seurat_clusters,
		numCores=40)
s3_2$SingleR <- s3_2$seurat_clusters
levels(s3_2$SingleR) <- sr.out.Luoma.cluster.2$labels


#### Bar graph ####
# Get the stacked barplot by proportion (Figure 1)
viz.data <- prop.table(table(Idents(s3_2),s3_2$status), margin = 2)
p <- Seurat::DimPlot(s3_2, reduction = "umap", do.return = TRUE) # Generate the tSNE plot, but save it as an object
pbuild <- ggplot2::ggplot_build(p) # Use ggpdlot_build to deconstruct the ggplot object
pdata <- pbuild$data[[1]] # Pull the data used for the plot
pdata <-  pdata[order(pdata$group), ] # Order the plot data by group
ucols <- unique(pdata$colour) # Get a vector of unique colors
ucols <- unname(ucols)
pdf("s3v2_barplot_gbStatus.pdf", width = 10, height = 10)
par(xpd=T, mar=c(5, 4, 4, 10) + 0.1)
barplot(viz.data, col=ucols, space=0.04, font.axis=2, 
        xlab="Status", ylab = "Proportion",legend.text = T,
        args.legend=list(x = ncol(viz.data)+2.5, y= max(colSums(viz.data))-0.35,bty = "n"))
dev.off()

# Get the dodged bar plot for each cluster with error bars (Figure 2)
#+++++++++++++++++++++++++
# Function to calculate the mean and the standard deviation
  # for each group
#+++++++++++++++++++++++++
# data : a data frame
# varname : the name of a column containing the variable
  #to be summarized
# groupnames : vector of column names to be used as
  # grouping variables
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE), sem = sd(x[[col]], na.rm=TRUE)/sqrt(length(x[[col]])))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
 return(data_sum)
}
wilcoxon.plot <- function(s3_3) {
	viz.data <- prop.table(table(Idents(s3_3),s3_3$orig.ident), margin = 2)

	inflamed <- c(1,3,6,7,10,12,14,16,18)
	status <- rep("not inflamed", ncol(viz.data))
	status[inflamed] <- "inflamed"

	df.prop <- data.frame(cluster= rep(rownames(viz.data),ncol(viz.data)), proportion= as.vector(viz.data), status= rep(status, each=nrow(viz.data)))
	df.prop.summary<- data_summary(df.prop, varname="proportion", groupnames=c("cluster", "status"))
	p <- ggplot(df.prop.summary, aes(x=cluster, y=proportion, fill=status)) + 
	   geom_bar(stat="identity", position=position_dodge()) +
	   geom_errorbar(aes(ymin=proportion-sem, ymax=proportion+sem), width=.2, position=position_dodge(.9)) +
	   labs(y = "Fraction of total cells") + scale_fill_brewer(palette="Paired") + theme_minimal()

	pdf("s3v2_barplot_sem_error.pdf", width = 10, height = 10)
	p
	dev.off()
}

#### Wilcoxon Test (Paired) ####
#res <- wilcox.test(proportion ~ status, data = subset(df.prop, cluster==0), paired = TRUE)
wilcoxon.p.values <- function(clusters) {
	p.values <- unlist(lapply(c(0,1,2,3,4,5,8,9),function(x) {wilcox.test(proportion ~ status, data = subset(df.prop, cluster==x), paired = TRUE)$p.value}))
	names(p.values) <- unique(df.prop$cluster)
	p.values
}


#### DGA between clusters (w/o ILC1 and ILC3 clusters) ####
library(future)
plan("multiprocess", workers = 40)
s3_3.markers <- FindAllMarkers(s3_3, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose=F)
top10 <- s3_3.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

#### DGA between inflamed and not inflamed clusters ####
library(future)
plan("multiprocess", workers = 40)
s3_3.inflamed.markers <- FindMarkers(s3_3, group.by="status", ident.1 = "inflamed", ident.2="not inflamed", min.pct = 0.25, logfc.threshold = 0.25, verbose=F)
top10 <- s3_2_regress.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)





#### Saving PDFs ####
pdf("s3v3_feature_CD4CD8_9.pdf", width = 20, height = 10)
FeaturePlot(s3_3, features= c("CD4","CD8A"))
dev.off()

pdf("s3v3_features_9.pdf", width = 14, height = 21)
FeaturePlot(s3_3, features= c("FOXP3","TNF", "GZMB", "IFNG", "CXCR3", "CXCR6"), ncol=2)
dev.off()

pdf("s3v3_PCelbow_3000hvf.pdf", width = 10, height = 10)
ElbowPlot(s3_3, ndims=40)
dev.off()

pdf("s3v3_heatmap.pdf", width = 38, height = 30)
DoHeatmap(s3_3, features = top10$gene) + NoLegend()
dev.off()

pdf("s3v2_feature.pdf", width = 20, height = 10)
FeaturePlot(s3_2, features= c("CD4","CD8A"))
dev.off()

pdf("s3v2_feature.pdf", width = 10, height = 10)
FeaturePlot(s3_2, features= "FOXP3")
dev.off()

pdf("s2v2_singler_labels_cluster_data.pdf", width = 10, height = 10)
DimPlot(s2_postQC, reduction = "umap", label=T, group.by="SingleR")
dev.off()

pdf("s3v2_singler_cluster_Luoma_data2.pdf", width = 10, height = 10)
DimPlot(s3_2, reduction = "umap", label=T, group.by="SingleR")
dev.off()

pdf("s3v2_singler_cluster_Luoma_scaled.pdf", width = 10, height = 10)
DimPlot(s3_4, reduction = "umap", label=T, group.by="seurat_clusters")
dev.off()

pdf("s2v2_postQC_UMAP_3_gbStatus.pdf", width = 10, height = 10)
DimPlot(s2_postQC, reduction = "umap", label=F, group.by="status")
dev.off()

pdf("s2v2_singler_labels_cluster.pdf", width = 10, height = 10)
DimPlot(s2_postQC, reduction = "umap", label=T, group.by="seurat_clusters")
dev.off()

pdf("s3v2_singler_labels_cluster.pdf", width = 10, height = 10)
DimPlot(s3_2, reduction = "umap", label=T)
dev.off()

pdf("s3v2_singler_labels_cluster.pdf", width = 10, height = 10)
DimPlot(s3_2, reduction = "umap", label=T)
dev.off()

pdf("s3v2_test.pdf", width = 10, height = 10)
DimPlot(s3_2, reduction = "umap", label=F, group.by="seurat_clusters")
dev.off()

pdf("s3v2_heatmap.pdf", width = 38, height = 30)
DoHeatmap(s3_2, features = top10$gene) + NoLegend()
dev.off()

# Save the heatmap
pdf("s2v2_heatmap.pdf", width = 48, height = 40)
DoHeatmap(s2_postQC, features = top10$gene) + NoLegend()
dev.off()
scp -r esk17@rcapps5.dfci.harvard.edu:/mnt/beegfs/home/esk17/live/*.pdf ~/Desktop/qc_files 
scp -r esk17@rcapps5.dfci.harvard.edu:/mnt/beegfs/home/esk17/data/martin/qc_files/s2_heatmap.pdf ~/Desktop/qc_files 
scp -r esk19@o2.hms.harvard.edu:/n/scratch2/ssuo/CD3/colitis/ ~/data/
scp -r esk19@o2.hms.harvard.edu:/n/scratch2/ssuo/CD3 ~/data/


rsync -r -z -v -e ssh esk19@o2.hms.harvard.edu:/n/scratch2/ssuo/CD3 ~/data/
