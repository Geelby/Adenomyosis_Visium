#reference:https://junjunlab.github.io/scRNAtoolVis-manual/jjvolcano.html

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
  library(ggpubr)

seurat_oobj <- readr::read_rds('//////.rds')

table(seurat_oobj$celltypenewaddsub)
table(seurat_oobj$intercelltypetwo)

seurat_obj1 <- seurat_oobj

##DEG
    seurat_integrated <- seurat_obj1
    celltypechoose <- "SFRP5epi"
    dir.create(celltypechoose)
    setwd(celltypechoose)
    
      DefaultAssay(seurat_integrated)<- "Spatial"  
      Idents(seurat_integrated) <- seurat_integrated$slice
      table(Idents(seurat_integrated))
      levels(Idents(seurat_integrated))

      Idents(seurat_integrated) <- seurat_integrated$celltypenewaddsub
      table(Idents(seurat_integrated))
      levels(Idents(seurat_integrated))
      
      seurat_integrated<- seurat_integrated[,Idents(seurat_integrated) %in% c('SFRP5+ epi')]
      table(Idents(seurat_integrated))
      Idents(seurat_integrated) <- seurat_integrated$intercelltypetwo
      table(Idents(seurat_integrated))
      table(seurat_integrated$intercelltypetwo)
      
      seurat.markers <- FindMarkers(object = seurat_integrated,  ident.1 = "AM", ident.2 = "UM",
                                       min.pct = 0.25, 
                                       thresh.use = 0.25)
      write.csv(seurat.markers,file='findmarker_seurat.markers.csv')

###jjvocano
  library(scRNAtoolVis)
  SFRP5epi <- read_csv('SFRP5epi/findmarker_seurat.markers.csv')
  SFRP5epi$subcelltype <- 'SFRP5+ epi'
  IGFBP3stro <- read_csv('IGFBP3stro/findmarker_seurat.markers.csv')
  IGFBP3stro$subcelltype <- 'IGFBP3+ stro'
  CNN1stro <- read_csv('CNN1stro/findmarker_seurat.markers.csv')
  CNN1stro$subcelltype <- 'CNN1+ stro'
  endo <- read_csv('endo/findmarker_seurat.markers.csv')
  endo$subcelltype <- 'Endothelial'
  immune <- read_csv('immune/findmarker_seurat.markers.csv')
  immune$subcelltype <- 'Immune'
  DESSMC <- read_csv('DESSMC/findmarker_seurat.markers.csv')
  DESSMC$subcelltype <- 'DES+ SMC'
  PV <- read_csv('PV/findmarker_seurat.markers.csv')
  PV$subcelltype <- 'PV'
  ESR1SMC <- read_csv('ESR1SMC/findmarker_seurat.markers.csv')
  ESR1SMC$subcelltype <- 'ESR1+ SMC'

  marker <-  rbind(SFRP5epi,IGFBP3stro,CNN1stro,endo,immune,DESSMC,ESR1SMC,PV)
  
  table(marker$subcelltype)
  
  marker=marker[(marker$p_val_adj<0.01),]
  table(marker$subcelltype)

  marker$subcelltype <- factor(as.character(marker$subcelltype), levels=c('SFRP5+ epi',
                                                                               'IGFBP3+ stro','CNN1+ stro',
                                                                               'Endothelial','Immune', 
                                                                               'DES+ SMC','ESR1+ SMC','PV'))
  marker$cluster <- marker$subcelltype
  
  marker$gene <- marker$...1
  
  mygene <- c('IHH','PTCH1','SOX17','HHIP',
              'IGFBP4','MMP11','ID1','TRH','C3',
              'SFRP1','CALR','CCDC80','CAS6',
              'GAS6','JUN','EMILIN1',
              'CCL5','CCL19',
              'ESR1')
  
  jjVolcano(diffData = marker,tile.col = color_celltype,size = 3.5,myMarkers = mygene,
            fontface = 'italic',log2FC.cutoff = 0.5)

#immune_proEMT
  seurat_integrated <- readRDS('/////.rds')
  setwd('figures/invade_env/')
  celltypechoose <- "immune"
  setwd(celltypechoose)
  
  DefaultAssay(seurat_integrated)<- "Spatial"  
  Idents(seurat_integrated) <- seurat_integrated$slice
  table(Idents(seurat_integrated))
  levels(Idents(seurat_integrated))
  
  Idents(seurat_integrated) <- seurat_integrated$celltypenewaddsub
  table(Idents(seurat_integrated))
  levels(Idents(seurat_integrated))

  seurat_integrated<- seurat_integrated[,Idents(seurat_integrated) %in% c('Immune')]
  table(Idents(seurat_integrated))

  g<-c('CCL19','CXCL14','NRG1','PGF')
  
  g <- unique(g)
  
  table(seurat_integrated$celltypenewaddsub)
  table(seurat_integrated$slice)
  
  Idents(seurat_integrated) <- seurat_integrated$slice
  table(Idents(seurat_integrated))

  seurat_integrated <- ScaleData(seurat_integrated)
  vln.df <- seurat_integrated@assays$Spatial@scale.data %>%
    t() %>%
    as.data.frame()%>%
    dplyr::select(g) %>% 
    rownames_to_column("CB") %>% 
    mutate(cluster = seurat_integrated$slice)%>%
    pivot_longer(cols = 2:(ncol(.)-1),
                 names_to = "gene",
                 values_to = "exp") %>% 
    mutate(gene = factor(gene,levels = g))

  my_color = slicecolor

  Idents(seurat_integrated) <- seurat_integrated$slice
  table(Idents(seurat_integrated))
  
  p1 <- ggplot(vln.df,aes(cluster,exp),color=factor(cluster))+
    geom_violin(aes(fill=cluster),scale = "width",draw_quantiles = c(0.25, 0.5, 0.75))+
    scale_fill_manual(values = my_color)+
    facet_grid(gene~.,scales = "free_y")+
    scale_y_continuous(expand = c(0,0))+
    theme_bw()+
    theme(
      panel.grid = element_blank(),
      axis.title.x.bottom = element_blank(),
      axis.ticks.x.bottom = element_blank(),
      axis.title.y.left = element_blank(),
      axis.ticks.y.left = element_blank(),
      axis.text.y.left = element_blank(),
      legend.position = "none",
      panel.spacing.y = unit(0, "cm"),
      strip.text.y = element_text(angle=0,size = 12,hjust = 0),
      strip.background.y = element_blank()
    )
  p1
  
  #pvalue
  expr <- seurat_integrated@assays$Spatial@scale.data
  gene_name <- c('PGF')
  gene_expression <- expr %>% 
    .[gene_name,] %>% 
    as.data.frame()
  colnames(gene_expression) <- paste0(gene_name)
  identical(colnames(seurat_integrated),row.names(gene_expression))
  seurat_integrated$PGF <- gene_expression[,paste0(gene_name)]  
  identical(seurat_integrated@meta.data[,paste0(gene_name)],gene_expression[,paste0(gene_name)])
  meta <- seurat_integrated@meta.data
  
  library(ggpubr)
  compare_means(PGF ~ slice,  data = seurat_integrated@meta.data, method = "t.test")
  seurat_integrated@meta.data$slice <- factor(seurat_integrated@meta.data$slice,levels=c('UM','AM1'))
  ggboxplot(seurat_integrated@meta.data, x = "slice", y = "PGF",
            color = "black", palette =c('#92C5DEFF','#FDDBC7FF'),
            add = "jitter",fill="slice",
            add.params=list(color = "black",size=0.3, shape = 19))+
    stat_compare_means(method = "t.test")
  
  ggsave( filename="boxplot_PGF_immune_env_UM_AM1.pdf",width = 3,height =5)
  table(seurat_integrated$celltypenewaddsub)
  table(seurat_integrated$slice)








