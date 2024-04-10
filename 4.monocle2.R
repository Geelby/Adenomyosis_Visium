#reference：https://cole-trapnell-lab.github.io/monocle-release/docs/

library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
options(stringsAsFactors = F)
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
library(tidyverse)
library(patchwork)
library(monocle)
library(paletteer)
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(ggsci)

#inputfprmonocle2
  Subset <- readr::read_rds('//////episub.rds')
  setwd("///////monocle/epi/")
  celltypechoose <- "Epithelial"
  
  Idents(Subset) <- Subset$epicelltype
  table(Idents(Subset))
  
  color_region = c('#FCDACAFF','#DF837DFF','#C91105FF')
  names(color_region)=c('LDHB+ epi','MMP11+ epi','SFRP5+ epi' )
  
  my_color <- color_region

  sce <- Subset
  DefaultAssay(sce) <- "Spatial"
  levels(Idents(sce))
  
  sce
#CellDataSet
  sample_ann <-  sce@meta.data  
  head(sample_ann)
  gene_ann <- data.frame(
    gene_short_name = rownames(sce@assays$Spatial) , 
    row.names =  rownames(sce@assays$Spatial) 
  )
  head(gene_ann)
  
  pd <- new("AnnotatedDataFrame",
            data=sample_ann)
  fd <- new("AnnotatedDataFrame",
            data=gene_ann)
  
  ct=as.data.frame(sce@assays$Spatial@counts)
  ct[1:4,1:4]
  
  dim(ct)
  dim(pd)
  dim(fd)
  
  sc_cds<- newCellDataSet(
    as.matrix(ct), 
    phenoData = pd,
    featureData =fd,
    expressionFamily = negbinomial.size(),
    lowerDetectionLimit=0.5)
  sc_cds
#Estimate size factors and dispersions
  sc_cds <- estimateSizeFactors(sc_cds)
  sc_cds <- estimateDispersions(sc_cds)
#Filtering 
  sc_cds <- detectGenes(sc_cds, min_expr = 0.1)
  print(head(fData(sc_cds)))
  expressed_genes <- row.names(subset(fData(sc_cds),
                                      num_cells_expressed >= 10))
#Constructing Single Cell Trajectories
    cds <- sc_cds
  
    deg.cluster <- FindAllMarkers(cds)
    
    diff_test_res <- arrange(deg.cluster, desc(avg_log2FC))
    ordering_genes <- diff_test_res %>% subset(p_val_adj<0.05) %>% row.names()
    express_genes <- ordering_genes[1:200]
    
    cds <- setOrderingFilter(cds, express_genes)
    
    plot_ordering_genes(cds)
    
    cds <- reduceDimension(cds, max_components = 2,
                           method = 'DDRTree')
    cds <- orderCells(cds)
    
    plot_cell_trajectory(cds, color_by = "epicelltype") +scale_color_manual(values = my_color)
    
    plot_cell_trajectory(cds, color_by = "epicelltype") + facet_wrap("~epicelltype",nrow=1)+scale_color_manual(values = my_color)
    
    plot_cell_trajectory(cds, color_by = "slice") 

    plot_cell_trajectory(cds, color_by = "slice") + facet_wrap("~slice",nrow=1)
    
    plot_cell_trajectory(cds, color_by = "State") 
    
    plot_cell_trajectory(cds, color_by = "State") + facet_wrap("~State",nrow=1)
    
    plot_cell_trajectory(cds, color_by = "Pseudotime") 
    
#find pseudo-temporal gene
celltypechoose <- "Epithelial"

color_region = c('#FCDACAFF','#DF837DFF','#C91105FF')
names(color_region)=c('LDHB+ epi','MMP11+ epi','SFRP5+ epi' )

my_color <- color_region

plot_cell_trajectory(cds, color_by = "epicelltype") +scale_color_manual(values = my_color)

plot_cell_trajectory(cds, color_by = "slice") +scale_color_manual(values = slicecolor)

plot_cell_trajectory(cds, color_by = "State") 

cds <- orderCells(cds, root_state = "3")
plot_cell_trajectory(cds, color_by = "Pseudotime") 

my_cds_subset=cds

BEAM_branch1 <-   BEAM(my_cds_subset, branch_point = 1, progenitor_method = 'duplicate')

BEAM_res = BEAM_branch1

p <- plot_genes_branched_heatmap(
  my_cds_subset[row.names(subset(BEAM_res, qval < 1e-4)),],
  branch_point = 1,
  num_clusters = 4, 
  use_gene_short_name = TRUE,
  cores = 4,
  show_rownames = F,
  return_heatmap = TRUE,
  branch_colors = c("#DCDCDC", "#92C5DEFF", "#F4A582FF"),
  branch_labels = c( "UM","AM"),
  hmcols = colorRampPalette(c('#4393C3FF','white','#D6604DFF'))(67))

my_branched_heatmap <- p
head(my_branched_heatmap$annotation_row)
table(my_branched_heatmap$annotation_row$Cluster) 
my_row <- my_branched_heatmap$annotation_row
my_row <- data.frame(cluster = my_row$Cluster,
                     gene = row.names(my_row),
                     stringsAsFactors = FALSE)
write.csv(my_row,paste0(celltypechoose,'gene_cluster.csv') )

#enrichment
seurat.markers <- read_csv(paste0(celltypechoose,'gene_cluster.csv'))

  sce.markers <- seurat.markers
  library(ggpubr)
  library(clusterProfiler)
  options(stringsAsFactors = F)
  suppressMessages(library('org.Hs.eg.db'))
  keytypes(org.Hs.eg.db)
  
  go_organism <- "org.Hs.eg.db" 
  kegg_organism <- "hsa"
  
  ids = bitr(sce.markers$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb=go_organism)
  sce.markers = merge(sce.markers, ids, by.x='gene', by.y='SYMBOL')
  
  #Adding Gene Annotations:load in an annotation file
  annotations <- read.csv("///////annotation.csv")
  
  # Combine markers with gene descriptions 
  sce.markers <- sce.markers %>% 
    left_join(y = unique(annotations[, c("gene_name", "description")]),
              by = c("gene" = "gene_name"))
  
  CSV <- data.frame(sce.markers)
  write.csv(CSV,'makers_annotation.csv')  
  
#GO
    gcSample=split(sce.markers$ENTREZID, sce.markers$cluster ) 
    xx <- compareCluster(gcSample,
                         fun = "enrichGO",
                         OrgDb = "org.Hs.eg.db",
                         ont = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.05
    )

      library(org.Hs.eg.db)
      datt=as.data.frame(xx)
      ensembl=strsplit(datt$geneID,"/")
      symbol=sapply(ensembl,function(x){
        y=bitr(x, fromType="ENTREZID", toType="SYMBOL", OrgDb= go_organism )
        y=y[!duplicated(y$ENTREZID),-1]
        y=paste(y,collapse = "/")
      })
      datt$geneID=symbol
      write.csv(datt,'enrichGOBP.csv') 
