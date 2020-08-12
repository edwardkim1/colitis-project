##################################
## cluster labeling: onlyTcells ##
##################################
library(SingleR)
library(BiocParallel)

s <- readRDS("saved_objects/onlyTcells.rds")

load('data/ssuo_CD3/seurat.object.RData')
Tcell <- process_ssuo_Tcells(Tcell)
martin <- read.csv(file="data/martin_categories.csv", header=TRUE,sep=",") %>% process_martin()
s <- label_singleR(s,save.diagnostics=T, save.name="saved_objects/onlyTcells_SingleR_diagnostics.rds")
require(ggpubr)
p1 <- DimPlot(s, reduction = "umap", label=T, group.by="SingleR.Martin")
p2 <- DimPlot(s, reduction = "umap", label=T, group.by="SingleR.Luoma")
p3 <- DimPlot(s, reduction = "umap", label=T, group.by="SingleR.Monaco")
ggarrange(p1,p2,p3, ncol=3, nrow=1)
ggsave("figures/onlyTcells_labeling.pdf", width= 12, height= 3, units= "in")
# Diagnostic plot