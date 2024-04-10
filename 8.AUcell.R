#reference:https://bioconductor.org/packages/release/bioc/vignettes/AUCell/inst/doc/AUCell.html

#0. Load scRNA-seq dataset and gene sets
  setwd("///////14_AUcell")
  seurat_obj =readRDS('//////seurat.rds')
  
  exprMatrix <- seurat_obj@assays$Spatial@counts
  exprMatrix <- as(exprMatrix, "dgCMatrix")
  dim(exprMatrix)
  
  library(GSEABase)
  library(AUCell)
  
  ECMM <- read.csv("ECM.csv",header = F) 
  genes <- as.vector(ECMM$V1)
  geneSetsone <- GeneSet(genes, setName="ECMscore")
  geneSetsone
  
  collagen <- read.csv("collagen.csv",header = F) 
  genes <- as.vector(collagen$V1)
  geneSetstwo <- GeneSet(genes, setName="collagenscore")
  geneSetstwo
  
  # Random
  set.seed(321)
  extraGeneSets <- c(
    GeneSet(sample(rownames(exprMatrix), 50), setName="Random (50g)"),
    GeneSet(sample(rownames(exprMatrix), 500), setName="Random (500g)"))
  
  countsPerGene <- apply(exprMatrix, 1, function(x) sum(x>0))
  # Housekeeping-like
  extraGeneSets <- c(extraGeneSets,
                     GeneSet(sample(names(countsPerGene)[which(countsPerGene>quantile(countsPerGene, probs=.95))], 100), setName="HK-like (100g)"))
  
  geneSets <- GeneSetCollection(c(geneSetsone,geneSetstwo,extraGeneSets))
  names(geneSets)
  geneSets
  
#1.Score gene signatures
  cells_AUC <- AUCell_run(exprMatrix, geneSets)
  save(cells_AUC, file="cells_AUC.RData")

#2. Determine the cells with the given gene signatures or active gene sets
  set.seed(333)
  par(mfrow=c(3,3)) 
  cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE) 

#3.Follow up examples
  library(clusterProfiler)
  library(ggplot2)
  library(Seurat)
  
  geneSet <- "ECMscore"
  aucsecm <- as.numeric(getAUC(cells_AUC)[geneSet, ])
  seurat_obj$AUCECM <- aucsecm
  
  geneSet <- "collagenscore"
  aucscollagen <- as.numeric(getAUC(cells_AUC)[geneSet, ])
  seurat_obj$AUCcollagen <- aucscollagen

    seurat_obj <- seurat_obj[,seurat_obj$slice %in% c('AM1','AM2','AM3')]
    Idents(seurat_obj) <- "toponecelltypenew"
    p1 <- SpatialDimPlot(seurat_obj, label = F, label.size = 3,
                         images = c('UM'), ncol = 1,pt.size.factor = 1.2,cols =my_color )+NoLegend()
    p2 <- SpatialDimPlot(seurat_obj, label = F, label.size = 3,
                         images = c('AM1'), ncol = 1,pt.size.factor = 1.2,cols =my_color )+NoLegend()
    p3 <- SpatialDimPlot(seurat_obj, label = F, label.size = 3,
                         images = c('AM2'), ncol = 1,pt.size.factor = 1.2,cols =my_color )+NoLegend()
    p4 <- SpatialDimPlot(seurat_obj, label = F, label.size = 3,
                         images = c('AM3'), ncol = 1,pt.size.factor = 1.2,cols =my_color )+NoLegend()
    p1|p2|p3|p4

    seurat_obj$slice <- factor(as.character(seurat_obj$slice), levels=c("AM1","AM2","AM3"))
    
    table(seurat_obj$slice)
    table(seurat_obj$toponecelltypenew)
    
    library(ggpubr)
    
    compare_means(AUCECM ~ toponecelltypenew,  data = seurat_obj@meta.data, method = "anova")
    
    seurat_obj$toponecelltypenew <- factor(seurat_obj$toponecelltypenew,levels=c('Epithelial','Stromal','Endothelial','Immune', 'SMC','PV'))
    
    # Change method to anova
    ggboxplot(seurat_obj@meta.data, x = "toponecelltypenew", y = "AUCECM",
              color = "black", palette =c('#C93F55FF', '#3C4B99FF', '#EACC62FF', '#469D76FF', '#DF9ED4FF','#924099FF'),
              add = "jitter",fill="toponecelltypenew",
              add.params=list(color = "black",size=0.1, shape = 19))+
      stat_compare_means(method = "anova")
    
    ggsave( filename="boxplot_ECM_sixcelltype_inAM.pdf",width = 3.5,height =5)
    table(seurat_obj$slice)
    table(seurat_obj$toponecelltypenew)
  




