#############################################################
## Updated Clustering Functions
##
## filter_stats # has to be done before scaling/normalizing
##
## edward_kim@college.harvard.edu - Jun. 26, 2020
############################################################
require(Seurat)


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
# filtered.object <- seurat.object[,which(lnf > lnf.lim & mt.fraction < mt.lim)]
# 

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

find.pK.noPlot <- function(sweep.stats) {
    "%ni%" <- Negate("%in%")
    if ("AUC" %ni% colnames(sweep.stats) == TRUE) {
        bc.mvn <- as.data.frame(matrix(0L, nrow = length(unique(sweep.stats$pK)), ncol = 5))
        colnames(bc.mvn) <- c("ParamID", "pK", "MeanBC", "VarBC", "BCmetric")
        bc.mvn$pK <- unique(sweep.stats$pK)
        bc.mvn$ParamID <- 1:nrow(bc.mvn)
        x <- 0
        for (i in unique(bc.mvn$pK)) {
            x <- x + 1
            ind <- which(sweep.stats$pK == i)
            bc.mvn$MeanBC[x] <- mean(sweep.stats[ind, "BCreal"])
            bc.mvn$VarBC[x] <- sd(sweep.stats[ind, "BCreal"])^2
            bc.mvn$BCmetric[x] <- mean(sweep.stats[ind, "BCreal"])/(sd(sweep.stats[ind,"BCreal"])^2)
        }
        return(bc.mvn)
    }
    if ("AUC" %in% colnames(sweep.stats) == TRUE) {
        bc.mvn <- as.data.frame(matrix(0L, nrow = length(unique(sweep.stats$pK)), ncol = 6))
        colnames(bc.mvn) <- c("ParamID", "pK", "MeanAUC", "MeanBC",
            "VarBC", "BCmetric")
        bc.mvn$pK <- unique(sweep.stats$pK)
        bc.mvn$ParamID <- 1:nrow(bc.mvn)
        x <- 0
        for (i in unique(bc.mvn$pK)) {
            x <- x + 1
            ind <- which(sweep.stats$pK == i)
            bc.mvn$MeanAUC[x] <- mean(sweep.stats[ind, "AUC"])
            bc.mvn$MeanBC[x] <- mean(sweep.stats[ind, "BCreal"])
            bc.mvn$VarBC[x] <- sd(sweep.stats[ind, "BCreal"])^2
            bc.mvn$BCmetric[x] <- mean(sweep.stats[ind, "BCreal"])/(sd(sweep.stats[ind,"BCreal"])^2)
        }
        return(bc.mvn)
    }
}

label_doublets <- function(seurat.object, nPCs= 30, doublet.formation.rate= 0.1, save.pK.figure = F, filename = "") {
	## Create Seurat Object ##
	require(DoubletFinder)
	require(dplyr)
	## pK identification
	seurat.object_pK <- paramSweep_v3(s, PCs = 1:nPCs) %>% summarizeSweep(GT = FALSE) %>% find.pK.noPlot()
	pK_max <- seurat.object_pK[which.max(seurat.object_pK$BCmetric),]$pK
	pK_max <- as.numeric(levels(pK_max)[as.integer(pK_max)])
	## Output pK maximization using ggplot2 as pdf
	seurat.object_pK$group <- numeric(length(seurat.object_pK$ParamID))
	if(save.pK.figure==TRUE) {
		p <- ggplot(data=seurat.object_pK, aes(x=pK, y=BCmetric, group=0)) +
	  		geom_line()+ scale_x_discrete(breaks=seq(0,0.3,0.05)) +
	  		geom_point() + theme_classic() 
		ggsave(filename)
	}
	## Computation of parameters
	homotypic.prop <- modelHomotypic(seurat.object@meta.data$seurat_clusters)
	nExp_poi <- round(doublet.formation.rate *nrow(seurat.object@meta.data)) #tailor doublet.formation.rate to your dataset
	nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
	## Run Doublet Finder algorithm
	seurat.object <- doubletFinder_v3(seurat.object, PCs = 1:nPCs, pN = 0.25, pK= pK_max, nExp = nExp_poi, reuse.pANN = FALSE)
	first_pANN <- names(seurat.object@meta.data)[grep("pANN", names(seurat.object@meta.data))]
	seurat.object <- doubletFinder_v3(seurat.object, PCs = 1:nPCs, pN = 0.25, pK = pK_max, nExp = nExp_poi.adj, reuse.pANN = first_pANN)
	DFs <- names(seurat.object@meta.data)[grep("DF.class", names(seurat.object@meta.data))]
	## Saving DF status in the seurat.object
	seurat.object@meta.data[,"DF_hi.lo"] <- seurat.object[[DFs[1]]]
	seurat.object@meta.data$DF_hi.lo[which(seurat.object@meta.data$DF_hi.lo == "Doublet" & seurat.object[[DFs[2]]] == "Singlet")] <- "Doublet_lo"
	seurat.object@meta.data$DF_hi.lo[which(seurat.object@meta.data$DF_hi.lo == "Doublet")] <- "Doublet_hi"
	## Returning the labeled object
	return(seurat.object)
}

VariableFeaturePlot.Tcells <- function(object, cols = c("black", "red"), pt.size = 1, log = NULL, selection.method = NULL, assay = NULL) {
    if (length(x = cols) != 2) {
        stop("'cols' must be of length 2")
    }
    hvf.info <- HVFInfo(object = object, assay = assay, selection.method = selection.method, status = TRUE)
    ### The following code was added by Edward Kim 6-30-2020 ###
    # Removes TCR variable genes and IG variable genes 
    '%ni%' = Negate('%in%')
	all.genes <- rownames(hvf.info) # do after importing counts matrix
	trvgenes <- all.genes[grepl(x=all.genes, pattern = "^TRAV|^TRBV|^TRGV|^TRDV")]
	igvgenes <- all.genes[grepl(x=all.genes, pattern = "^IG.V")]
	hvf.info <- hvf.info[all.genes %ni% c(trvgenes,igvgenes),]
	### End ###
    var.status <- c("no", "yes")[unlist(x = hvf.info[, ncol(x = hvf.info)]) + 1]
    hvf.info <- hvf.info[, c(1, 3)]
    axis.labels <- switch(EXPR = colnames(x = hvf.info)[2], variance.standarized = c("Average Expression",
        "Standardized Variance"), dispersion.scaled = c("Average Expression", "Dispersion"), residual_variance = c("Geometric Mean of Expression",
        "Residual Variance"))
    log <- log %||% (any(c("variance.standardized", "residual_variance") %in% colnames(x = hvf.info)))
    plot <- SingleCorPlot(data = hvf.info, col.by = var.status, pt.size = pt.size)
    plot <- plot + labs(title = NULL, x = axis.labels[1], y = axis.labels[2]) + scale_color_manual(labels = paste(c("Non-variable", "Variable"), "count:", table(var.status)), values = cols)
	if (log) {plot <- plot + scale_x_log10()}
    return(plot)
}
environment(VariableFeaturePlot.Tcells) <- environment(ElbowPlot)

VariableFeaturePlot.Tcells.SCT <- function(object, cols = c("black", "red"), pt.size = 1, log = NULL, selection.method = NULL, assay = NULL) {
    if (length(x = cols) != 2) {
        stop("'cols' must be of length 2")
    }
    hvf.info <- HVFInfo(object = object, assay = assay, selection.method = selection.method, status = TRUE)
    ### The following code was added by Edward Kim 6-30-2020 ###
    # Removes TCR variable genes and IG variable genes 
    hvf.info$variable <- rownames(hvf.info) %in% VariableFeatures(object)
    '%ni%' = Negate('%in%')
	all.genes <- rownames(hvf.info) # do after importing counts matrix
	trvgenes <- all.genes[grepl(x=all.genes, pattern = "^TRAV|^TRBV|^TRGV|^TRDV")]
	igvgenes <- all.genes[grepl(x=all.genes, pattern = "^IG.V")]
	hvf.info <- hvf.info[all.genes %ni% c(trvgenes,igvgenes),]
	### End ###
    var.status <- c("no", "yes")[unlist(x = hvf.info[, ncol(x = hvf.info)]) + 1]
    hvf.info <- hvf.info[, c(1, 3)]
    axis.labels <- switch(EXPR = colnames(x = hvf.info)[2], variance.standarized = c("Average Expression",
        "Standardized Variance"), dispersion.scaled = c("Average Expression", "Dispersion"), residual_variance = c("Geometric Mean of Expression",
        "Residual Variance"))
    log <- log %||% (any(c("variance.standardized", "residual_variance") %in% colnames(x = hvf.info)))
    plot <- SingleCorPlot(data = hvf.info, col.by = var.status, pt.size = pt.size)
    plot <- plot + labs(title = NULL, x = axis.labels[1], y = axis.labels[2]) + scale_color_manual(labels = paste(c("Non-variable", "Variable"), "count:", table(var.status)), values = cols)
	if (log) {plot <- plot + scale_x_log10()}
    return(plot)
}
environment(VariableFeaturePlot.Tcells.SCT) <- environment(ElbowPlot)

cluster_umap <- function(seurat.object, npc, k.param, resolution, min.dist) {
	require(future)
	plan("multicore")
	seurat.object <- FindNeighbors(seurat.object, dims = 1:npc, k.param = k.param) %>% FindClusters(resolution = resolution, graph.name ="RNA_snn") %>% RunUMAP(min.dist = min.dist, graph="RNA_snn")
	return(seurat.object)
}

process_martin <- function(martin) {
	# Prepare Martin 2019 reference
	df.martin <- as.data.frame(martin)
	filtered.martin <- df.martin[which(df.martin$X != "28-Feb" & df.martin$X != "1-Mar"),] # need to remove the rows with "28-Feb" "1-Mar"
	rownames(filtered.martin) <- filtered.martin$X
	martin <- as.matrix(filtered.martin)
	martin <- martin[, -1]
	class(martin) <- "numeric"
	return(martin)
}

process_ssuo_Tcells <- function(Tcell) {
	# Prepare Luoma reference
	Tcell$names <- Tcell$seurat_clusters
	levels(Tcell$names) <- c("Type 1 cytokines Trm", "CD8 TRM", "IEL: CD8/gd", "Cytotoxic effector", "LP Trm", "Type 3 cytokines Trm", "Treg", "Naive/CM", "TFH","Th1 effector", "Cycling", "5", "MAIT")
	return(Tcell)
}

label_singleR <- function(seurat.object, save.diagnostics=F, save.name=NULL) {
 	require(BiocParallel)
 	pred.CPI <- SingleR(method='cluster',
		test=seurat.object[['RNA']]@data, 
		clusters= seurat.object$seurat_clusters,
		ref=as.matrix(Tcell[['RNA']]@data), 
		labels=Tcell$names,
		de.method="wilcox",
		BPPARAM=MulticoreParam())
 	seurat.object$SingleR.Luoma <- seurat.object$seurat_clusters
	levels(seurat.object$SingleR.Luoma) <- pred.CPI$labels

	pred.Martin <- SingleR(method='cluster',
		test=seurat.object[['RNA']]@data, 
		clusters= seurat.object$seurat_clusters,
		ref=martin, 
		labels=colnames(martin),
		BPPARAM=MulticoreParam())
	seurat.object$SingleR.Martin <- seurat.object$seurat_clusters
	levels(seurat.object$SingleR.Martin) <- pred.Martin$labels

	# ref <- MonacoImmuneData()
	# pred.Monaco <- SingleR(method="cluster",
	# 	test=s[['RNA']]@data,
	# 	clusters= seurat.object$seurat_clusters, 
	# 	ref=ref, 
	# 	labels=ref$label.fine, 
	# 	BPPARAM=MulticoreParam())
	# seurat.object$SingleR.Monaco <- seurat.object$seurat_clusters
	# levels(seurat.object$SingleR.Monaco) <- pred.Monaco$labels

	if(save.diagnostics) {
		saveRDS(list(pred.CPI, pred.Martin),save.name)
	}

	return(seurat.object)
}

add_labels() <- function(seurat.object, labels) {


}


save_UMAP_figures <- function(seurat.object) {
	require(gridExtra)
	filename <- paste(seurat.object,".pdf",sep="")
	pdf(filename, width = 10, height = 10)
		p.umap <- DimPlot(seurat.object, reduction = "umap", label=T)
		p1 <- DimPlot(seurat.object, reduction = "umap", label=T, group.by="SingleR.Martin")
		p2 <- DimPlot(seurat.object, reduction = "umap", label=T, group.by="SingleR.Luoma")
		p.status <- DimPlot(seurat.object, reduction = "umap", label=F, group.by="status")
		p.patient <- DimPlot(seurat.object, reduction = "umap", label=F, group.by="patient.id")

		print(p.elbow)
		print(p.umap)
		print(p.status)
		print(p.patient)
		print(grid.arrange(p1, p2, nrow = 1))
	dev.off()
}





