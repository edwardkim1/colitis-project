#############################################################
## Hypothesis Testing Functions
##
## data_summary
## wilcoxon.df
## wilcoxon.p.values
##
## edward_kim@college.harvard.edu - Dec. 20, 2019
############################################################

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

wilcoxon.df <- function(seurat.object) {
	viz.data <- prop.table(table(Idents(seurat.object),seurat.object$orig.ident), margin = 2)

	inflamed <- c(1,3,6,7,10,12,14,16,18)
	status <- rep("not inflamed", ncol(viz.data))
	status[inflamed] <- "inflamed"

	df.prop <- data.frame(cluster= rep(rownames(viz.data),ncol(viz.data)), proportion= as.vector(viz.data), status= rep(status, each=nrow(viz.data)))
	return(df.prop)
}

#### Wilcoxon Test (Paired) ####
#res <- wilcox.test(proportion ~ status, data = subset(df.prop, cluster==0), paired = TRUE)
wilcoxon.p.values <- function(clusters, df.prop) {
	p.values <- unlist(lapply(clusters,function(x) {wilcox.test(proportion ~ status, data = subset(df.prop, cluster==x), paired = TRUE)$p.value}))
	names(p.values) <- unique(df.prop$cluster)
	return(p.values)
}