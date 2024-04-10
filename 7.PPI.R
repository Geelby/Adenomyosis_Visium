#reference：https://cytoscape.org/

library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
options(stringsAsFactors = F)
library(clustree)
library(cowplot)
library(tidyverse)
library(paletteer)
library(ggpubr)

    seurat_integrated <- readRDS('//////.rds')
    celltypechoose <- "SFRP5epi"
    setwd('../invade_env/SFRP5epi/')
    
    DefaultAssay(seurat_integrated)<- "Spatial"  
    Idents(seurat_integrated) <- seurat_integrated$slice
    table(Idents(seurat_integrated))
    levels(Idents(seurat_integrated))
    
    Idents(seurat_integrated) <- seurat_integrated$celltypenewaddsub
    table(Idents(seurat_integrated))
    levels(Idents(seurat_integrated))
  
    seurat_integrated<- seurat_integrated[,Idents(seurat_integrated) %in% c('SFRP5+ epi')]
    table(Idents(seurat_integrated))
    #DEG
      DefaultAssay(seurat_integrated)<- "Spatial"  
      
      Idents(seurat_integrated) <- seurat_integrated$celltypenewaddsub
      table(Idents(seurat_integrated))
      levels(Idents(seurat_integrated))
      
      Idents(seurat_integrated) <- seurat_integrated$slice
      table(Idents(seurat_integrated))
      
      seurat.markers <- FindAllMarkers(object = seurat_integrated, only.pos = TRUE, 
                                       min.pct = 0.25, 
                                       thresh.use = 0.25)
      
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
      
      #sce.markerssave <- sce.markers
      sce.markers <- sce.markerssave
      sce.markers=sce.markers[(sce.markers$avg_log2FC>0.263034) & (sce.markers$p_val_adj<0.05),]
   
    #GO
      gcSample=split(sce.markers$ENTREZID, sce.markers$cluster) 
      xx <- compareCluster(gcSample,
                           fun = "enrichGO",
                           OrgDb = "org.Hs.eg.db",
                           ont = "BP",
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.05)

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

    #CYTOSCAPE
      termm <- read_csv('epi_IHH.csv')
      
      seurat.markers$termm = "NA"
      for(i in 1:nrow(termm)){
        seurat.markers[which( seurat.markers$...1  == termm$gene[i]),'termm'] <- termm$term[i]}
      table(seurat.markers$termm)
      write.csv(seurat.markers,file='epi_seurat.markers_addterm.csv')
      
      seurat.markers <- seurat.markers[seurat.markers$termm != 'NA',]
      write.csv(seurat.markers,file='epi_seurat.markers_onlyterm.csv')
      
      table(seurat.markers$termm)
      termm <- termm[(termm$gene %in% seurat.markers$...1),]
      table(termm$term)
      
      string=read.table('string_interactions_short.tsv',header = T)
      string$interaction='gene'
      head(string)
      
      AA <- read_csv('genechoose.csv')
      termm <- termm[(termm$gene %in% AA$...2),]
      
      GOPPI=termm[,c(2,1)]
      head(GOPPI)
      
      colnames(GOPPI)=c('node1','node2')
      GOPPI$combined_score=NA
      GOPPI$interaction='termm'
      string=rbind(string,GOPPI)
      write.table(string,'string.csv',sep = ',',row.names = F,quote = F)
      node=data.frame(node=unique(c(string$node1,string$node2)))
      node$Log2FC=seurat.markers$avg_log2FC[match(node$node,seurat.markers$...1)]
      
      top5 <- read_csv('enrichGOBP.csv')
      node$P=top5$p.adjust[match(node$node,top5$Description)]

      node$category[node$node%in%seurat.markers$...1]='gene'
      node$category[node$node%in%top5$Description]='GO'
      write.table(node,'node.csv',quote = F,sep = ',',row.names = F)
    
    
