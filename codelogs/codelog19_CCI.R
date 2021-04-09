##############################
# 0 - Load librairies
##############################
library(dplyr)
library(stringr)
library(Seurat)
library(future)
library(ggplot2)
#library(gridExtra)
#library(grid)
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
options(ggrepel.max.overlaps = Inf)

############################## 
# 2 - Source file 
##############################
source("scripts/Stat111functions.R")
source("scripts/new_clustering_functions.R")


############################
##  CellphoneDB attempt
############################
# pip3 install cellphonedb numpy==1.19.5
conda create -n cpdb python=3.7
## for debugging:
#lsof +D ./ | awk '{print $2}' | tail -n +2 | xargs -r kill -9


#export PATH="/home/sv467/anaconda3/bin:$PATH"
export PATH="/home/sv467/anaconda3/env/cpdb/bin:$PATH"


# CellphoneDB trial 1
data_object <- s
# Take raw data and normalize it.
count_norm <- data_object[["RNA"]]@data
write.table(count_norm,"cpdb/trial1/cellphonedb_count.txt", sep="\t", quote=F)
# Generating metadatafile.
meta_data <- cbind(rownames(data_object@meta.data), data_object@meta.data[,"seurat_clusters", drop=F])
colnames(meta_data) <- c("cell","cell_type")
# Cluster is the user’s corresponding cluster column.
write.table(meta_data,"cpdb/trial1/cellphonedb_meta.txt",  sep="\t",  quote=F, row.names=F)

#83,000 cells - estimated time is 1.5hr * 8 = 12 hours on 4 e.t.a. 1am-1pm
# cellbender

cellphonedb method analysis test_meta.txt test_counts.txt


############################
##  CellChat attempt
############################
library(CellChat)
library(patchwork)
#load all cells
s <- readRDS("saved_objects/UC_chang_021321/uccpirec-allCD45.rds")
# split
s.list <- SplitObject(s[,s$colitis == "UC" | s$colitis == "Colitis"], "orig.ident")

# extract the data needed for CCI
data.input <- lapply(s.list, function(x) GetAssayData(x, assay = "RNA", slot = "data")) # normalized data matrix
labels <- lapply(s.list, function(x) x$SingleR.Luoma2)
meta <- lapply(labels, function(x) data.frame(group = x, row.names = names(x))) # create a dataframe of the cell labels
cellchat.list <- lapply(1:2, function(x) createCellChat(object = data.input[[x]], meta = meta[[x]], group.by = "group"))
names(cellchat.list) <- c("CPI-colitis", "UC")


future::plan("multiprocess", workers = 20) # do parallel

apply.DB <- function(cellchat, db.subset = NULL) {
	CellChatDB <- CellChatDB.human
	if(is.null(db.subset)) {
		CellChatDB.use <- CellChatDB
	} else {
		CellChatDB.use <- subsetDB(CellChatDB, search = db.subset)
	}	
	cellchat@DB <- CellChatDB.use
	cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
	cellchat <- identifyOverExpressedGenes(cellchat)
	cellchat <- identifyOverExpressedInteractions(cellchat)
	cellchat <- projectData(cellchat, PPI.human)
	cellchat
}


# "Secreted Signaling"
cellchat.list <- lapply(cellchat.list, apply.DB)


apply.computeProb <- function(cellchat, raw.use=T) {
	cellchat <- computeCommunProb(cellchat, raw.use = raw.use)
	# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
	cellchat <- filterCommunication(cellchat, min.cells = 10)
	cellchat <- computeCommunProbPathway(cellchat)
	cellchat <- aggregateNet(cellchat)
	cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
	cellchat
}

cellchat.list <- lapply(cellchat.list, apply.computeProb)


#cellchat.list2 <- lapply(cellchat.list2, function(cellchat) {
#cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
#cellchat
#	})

# merge cells
cellchat <- mergeCellChat(cellchat.list, add.names = names(cellchat.list))


# figure 1
lapply(1:2, function(j) {
	groupSize <- as.numeric(table(cellchat.list[[j]]@idents))
	mat <- cellchat.list[[j]]@net$weight
	pdf(paste("figures/UCCPI_CellChat_030921/cci_fig1-",j, ".pdf", sep=""), width=3*2, height=3*4)
	par(mfrow = c(4,2), xpd=TRUE)
	for (i in 1:nrow(mat)) {
	  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
	  mat2[i, ] <- mat[i, ]
	  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
	}
	dev.off()
})

# figure 2
object.list <- cellchat.list
pathways.show <- c("TNF") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
pdf(paste("figures/UCCPI_CellChat_030921/cci_fig2.pdf", sep=""), width=4*2, height=4)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()


#figure 3
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2
ggsave("figures/UCCPI_CellChat_030921/cci_fig3.pdf", width=4*2, height=4,units="in")

# redo with raw.use=F
cellchat.list2 <- lapply(1:2, function(x) createCellChat(object = data.input[[x]], meta = meta[[x]], group.by = "group"))
names(cellchat.list2) <- c("CPI-colitis", "UC")
cellchat.list2 <- lapply(cellchat.list2, apply.DB)
cellchat.list2 <- lapply(cellchat.list2, function(x) apply.computeProb(x,raw.use=F))
cellchat <- mergeCellChat(cellchat.list2, add.names = names(cellchat.list2))


# figure 1rawF
lapply(1:2, function(j) {
	groupSize <- as.numeric(table(cellchat.list2[[j]]@idents))
	mat <- cellchat.list2[[j]]@net$weight
	pdf(paste("figures/UCCPI_CellChat_030921/cci_fig1-",j, "rawF.pdf", sep=""), width=3*2, height=3*4)
	par(mfrow = c(4,2), xpd=TRUE)
	for (i in 1:nrow(mat)) {
	  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
	  mat2[i, ] <- mat[i, ]
	  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
	}
	dev.off()
})

# figure 2rawF
object.list <- cellchat.list2
pathways.show <- c("TNF") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
pdf(paste("figures/UCCPI_CellChat_030921/cci_fig2rawF.pdf", sep=""), width=4*2, height=4)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()


#figure 3rawF
gg1 <- compareInteractions(cellchat2, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat2, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2
ggsave("figures/UCCPI_CellChat_030921/cci_fig3rawF.pdf", width=4*2, height=4,units="in")

# figure 4rawF
object.list <- cellchat.list2
pathways.show <- c("CXCL") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
pdf(paste("figures/UCCPI_CellChat_030921/cci_fig4rawF.pdf", sep=""), width=4*2, height=4)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()

# figure 4rawT-CADM
object.list <- cellchat.list
pathways.show <- c("CADM") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
pdf(paste("figures/UCCPI_CellChat_030921/cci_fig4rawT-CADM.pdf", sep=""), width=4*2, height=4)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()


# figure 5rawF
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
netAnalysis_computeCentrality()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
pdf(paste("figures/UCCPI_CellChat_030921/cci_fig5rawF.pdf", sep=""), width=4*2, height=4)
patchwork::wrap_plots(plots = gg)
dev.off()


# figure 6rawF
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "functional")
#> Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
pdf(paste("figures/UCCPI_CellChat_030921/cci_fig6rawF.pdf", sep=""), width=5*1.1, height=5)
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
dev.off()
#> 2D visualization of signaling networks from datasets 1 2


# figure 7rawF
pdf(paste("figures/UCCPI_CellChat_030921/cci_fig7rawF.pdf", sep=""), width=3*1.1, height=4)
rankSimilarity(cellchat, type = "functional")
dev.off()

# figure 8rawF
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
pdf(paste("figures/UCCPI_CellChat_030921/cci_fig8rawF.pdf", sep=""), width=4*2, height=9)
gg1 + gg2
dev.off()

# figure 9rawF-outgoing
library(ComplexHeatmap)
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 20)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 20)
pdf(paste("figures/UCCPI_CellChat_030921/cci_fig9rawF-outgoing.pdf", sep=""), width=4*2, height=9)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()


# figure 9rawF-incoming
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 20, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 20, color.heatmap = "GnBu")
pdf(paste("figures/UCCPI_CellChat_030921/cci_fig9rawF-incoming.pdf", sep=""), width=4*2, height=11)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

# figure 10rawF
pdf(paste("figures/UCCPI_CellChat_030921/cci_fig10rawF.pdf", sep=""), width=4*2, height=4*5)
netVisual_bubble(cellchat, sources.use = 5, targets.use = c(1:7),  comparison = c(1, 2), angle.x = 45)
dev.off()

# figure 10rawF-diff
gg1 <- netVisual_bubble(cellchat, sources.use = 5, targets.use = c(1:7),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in UC", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = 5, targets.use = c(1:7),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in UC", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
pdf(paste("figures/UCCPI_CellChat_030921/cci_fig10rawF-diff.pdf", sep=""), width=4*4, height=4*5)
gg1 + gg2
dev.off()


# figure 11rawF
pdf(paste("figures/UCCPI_CellChat_030921/cci_fig11rawF.pdf", sep=""), width=4*2, height=4*2)
plotGeneExpression(cellchat, signaling = "CXCL", split.by = "datasets", colors.ggplot = T)
dev.off()

# figure 12rawF
pdf(paste("figures/UCCPI_CellChat_030921/cci_fig12rawF.pdf", sep=""), width=4*2, height=4*2)
plotGeneExpression(cellchat, signaling = "TNF", split.by = "datasets", colors.ggplot = T)
dev.off()
# save merged cellchat object
saveRDS(cellchat, file = "saved_objects/UCCPI_CellChat_030921/cellchat_comparisonAnalysis_UCvsCPI-colitis.rds")



# figure 12rawF - CADM
cellchat <- readRDS("saved_objects/UCCPI_CellChat_030921/cellchat_comparisonAnalysis_UCvsCPI-colitis.rds")

pdf(paste("figures/UCCPI_CellChat_030921/cci_fig12rawF-CADM.pdf", sep=""), width=4*2, height=4*2)
plotGeneExpression(cellchat, signaling = "CADM", split.by = "datasets", colors.ggplot = T)
dev.off()


############################
##  CellChat attempt comparing controls
############################
library(CellChat)
library(patchwork)
#load all cells
s <- readRDS("saved_objects/UC_chang_021321/uccpirec-allCD45.rds")
# split
s.list <- SplitObject(s[,s$colitis3 == "UC-Control" | s$colitis3 == "Control"], "orig.ident")

# extract the data needed for CCI
data.input <- lapply(s.list, function(x) GetAssayData(x, assay = "RNA", slot = "data")) # normalized data matrix
labels <- lapply(s.list, function(x) x$SingleR.Luoma2)
meta <- lapply(labels, function(x) data.frame(group = x, row.names = names(x))) # create a dataframe of the cell labels
cellchat.list2 <- lapply(1:2, function(x) createCellChat(object = data.input[[x]], meta = meta[[x]], group.by = "group"))
names(cellchat.list2) <- c("CPI-colitis", "UC")
cellchat.list2 <- lapply(cellchat.list2, apply.DB)
cellchat.list2 <- lapply(cellchat.list2, function(x) apply.computeProb(x,raw.use=T))
cellchat2 <- mergeCellChat(cellchat.list2, add.names = names(cellchat.list2))


# figure 4rawT-CADM
object.list <- cellchat.list2
pathways.show <- c("CADM") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
pdf(paste("figures/UCCPI_CellChat_030921/cci_fig4rawT-CADM-controls.pdf", sep=""), width=4*2, height=4)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 2) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()


############################
##  CellChat attempt comparing controls
############################
 library(CellChat)
library(patchwork)
