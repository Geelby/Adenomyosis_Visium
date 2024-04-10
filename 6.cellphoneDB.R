#reference：https://github.com/ventolab/CellphoneDB/tree/master/notebooks
  
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
  
#dataprepare
    sooo =readRDS('//////.rds')
    seurat_oobj <- sooo
    table(seurat_oobj@meta.data$celltypenewaddsub)
    
    sooo$annotation <- paste0(sooo$slice, '_', sooo$celltypenewaddsub)
    
    Idents(sooo) = sooo$celltypenewaddsub
    
    table(sooo$celltypenewaddsub)
    table(Idents(sooo))
    
    table(sooo$annotation)
    
    seurat_obj <- sooo
    
    ###1. Load seurat object
    
    Idents(seurat_oobj) <- 'slice'

    so <- seurat_oobj
    
    ###AM1
      dir.create("AM1")  
      dir.create("AM1/counts_mtx")  
      setwd('AM1/counts_mtx')
      
      so<- so[,so$slice %in% c('AM1')]
      
      Idents(so) = so$annotation
      
      table(so$slice)
      table(Idents(so))
      
      ###2. Write gene expression in mtx format
      # Save normalised counts - NOT scaled!
      writeMM(so@assays$Spatial@data, file = 'matrix.mtx')
      # save gene and cell names
      write(x = rownames(so@assays$Spatial@data), file = "features.tsv")
      write(x = colnames(so@assays$Spatial@data), file = "barcodes.tsv")
      
      ###3. Generate your meta
      table(so@meta.data$annotation)
      
      so@meta.data$Cell = rownames(so@meta.data)
      df = so@meta.data[, c('Cell', 'annotation')]
      write.table(df, file ='../meta.tsv', sep = '\t', quote = F, row.names = F)
      
      table(Idents(so))
      ###4. Compute DEGs (optional)
      
      p <- SpatialDimPlot(so, label = F, label.size = 3,
                          images = c("UM","AM1","AM2","AM3"), ncol = 2,crop = F)
      
      p

      ggsave(plot=p, filename="../check.png",width = 15,height = 15)
      
      DEGs <- FindAllMarkers(object = so, logfc.threshold = 0.2,test.use = "LR",verbose = T,only.pos = T,random.seed = 1,min.pct = 0.1,return.thresh = 0.05)
      
      'DKK1' %in% rownames(so@assays$Spatial@counts)
      
      fDEGs = subset(DEGs, p_val_adj < 0.05 & avg_log2FC > 0.1)
      
      # 1st column = cluster; 2nd column = gene 
      fDEGs = fDEGs[, c('cluster', 'gene', 'p_val_adj', 'p_val', 'avg_log2FC', 'pct.1', 'pct.2')] 
      write.table(fDEGs, file ='../DEGs.tsv', sep = '\t', quote = F, row.names = F)
      
      head(fDEGs)

########draw picture
  library(ktplots)
  library(Seurat)
  library(SeuratData)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  
  # one you have set that up correctly, you can then read in the files.
  AM3_means <- read.delim("///////lesion_env/AM3/output/method2/statistical_analysis_means_03_10_2024_09:24:50.txt", check.names = FALSE)
  AM3_pvals <- read.delim("///////lesion_env/AM3/output/method2/statistical_analysis_pvalues_03_10_2024_09:24:50.txt", check.names = FALSE)
  AM3_decon <- read.delim("///////lesion_env/AM3/output/method2/statistical_analysis_deconvoluted_03_10_2024_09:24:50.txt", check.names = FALSE)
  
  AM2_means <- read.delim("///////lesion_env/AM2/output/method2/statistical_analysis_means_03_10_2024_09:23:54.txt", check.names = FALSE)
  AM2_pvals <- read.delim("///////lesion_env/AM2/output/method2/statistical_analysis_pvalues_03_10_2024_09:23:54.txt", check.names = FALSE)
  AM2_decon <- read.delim("///////lesion_env/AM2/output/method2/statistical_analysis_deconvoluted_03_10_2024_09:23:54.txt", check.names = FALSE)
  
  AM1_means <- read.delim("///////lesion_env/AM1/output/method2/statistical_analysis_means_03_10_2024_09:22:05.txt", check.names = FALSE)
  AM1_pvals <- read.delim("///////lesion_env/AM1/output/method2/statistical_analysis_pvalues_03_10_2024_09:22:05.txt", check.names = FALSE)
  AM1_decon <- read.delim("///////lesion_env/AM1/output/method2/statistical_analysis_deconvoluted_03_10_2024_09:22:05.txt", check.names = FALSE)
  
  
  means <- combine_cpdb(AM1_means,AM2_means, AM3_means)
  pvals <- combine_cpdb(AM1_pvals, AM2_pvals,AM3_pvals)
  decon <- combine_cpdb(AM1_decon,AM2_decon, AM3_decon)
  
  sooo =readRDS('///////.rds')
  seurat_oobj <- sooo
  
  table(seurat_oobj@meta.data$celltypenewaddsub)
  
  Idents(seurat_oobj) = seurat_oobj$celltypenewaddsub
  
  color_celltype = c('#FCDACAFF','#DF837DFF','#C91105FF',
                     '#C6DBEFFF','#6BAED6FF','#08519CFF',
                     '#EACC62FF', '#469D76FF', 
                     '#E8D8E8FF','#583070FF',
                     '#924099FF')
  
  names(color_celltype)=c('LDHB+ epi','MMP11+ epi','SFRP5+ epi',
                          'TXN+ stro','IGFBP3+ stro','CNN1+ stro',
                          'Endothelial','Immune', 
                          'DES+ SMC','ESR1+ SMC','PV')
  
  my_color = color_celltype
  
  p2 <- SpatialDimPlot(seurat_oobj, label = F, label.size = 3,
                       images = c("UM","AM1","AM2","AM3"), ncol = 2,crop = F,cols =my_color )
  p2
  
  sooo <- seurat_oobj
  
  sooo
  
  table(sooo$slice)
  table(Idents(sooo))
  
  sooo$annotation <- paste0(sooo$slice, '_', sooo$celltypenewaddsub)
  
  Idents(sooo) = sooo$celltypenewaddsub
  
  table(sooo$celltypenewaddsub)
  table(Idents(sooo))
  
  table(sooo$annotation)
  
  seurat_obj <- sooo
  
  so<- sooo
  
  so$annotation <- factor(as.character(so$annotation), levels=c('AM1_LDHB+ epi', 'AM1_MMP11+ epi','AM1_SFRP5+ epi',"AM1_TXN+ stro","AM1_IGFBP3+ stro","AM1_CNN1+ stro",'AM1_Endothelial','AM1_Immune','AM1_DES+ SMC','AM1_ESR1+ SMC','AM1_PV',
                                                                'AM2_LDHB+ epi', 'AM2_MMP11+ epi','AM2_SFRP5+ epi',"AM2_TXN+ stro","AM2_IGFBP3+ stro","AM2_CNN1+ stro",'AM2_Endothelial','AM2_Immune','AM2_DES+ SMC','AM2_ESR1+ SMC','AM2_PV',
                                                                'AM3_LDHB+ epi', 'AM3_MMP11+ epi','AM3_SFRP5+ epi',"AM3_TXN+ stro","AM3_IGFBP3+ stro","AM3_CNN1+ stro",'AM3_Endothelial','AM3_Immune','AM3_DES+ SMC','AM3_ESR1+ SMC','AM3_PV'))
  
  AA <- c('TGFB1-TGFbeta-receptor1','TGFB1-TGFBR3','TGFB1-TGFbeta-receptor2',
    'FN1-integrin-a11b1-complex','FN1-integrin-a5b1-complex','FN1-integrin-aVb5-complex',
    'FN1-integrin-aVb1-complex','FN1-integrin-a3b1-complex',
    'PDGFB-PDGFR-complex','PDGFB-PDGFRA','PDGFB-PDGFRB',
    'VWF-Glycoprotein-Ib-complex','VWF-integrin-aVb3-complex',
    'LAMA2-ADGRG6','PLAU-PLAUR','CD24-SELP')
  
  a <-  plot_cpdb(cell_type1 = 'CNN1+ stro', cell_type2 = 'SFRP5+ epi|IGFBP3+ stro|CNN1+ stro|Endothelial|Immune|ESR1+ SMC|PV', scdata = so,
                  idents = 'celltypenewaddsub', # column name where the cell ids are located in the metadata
                  split.by = 'slice', # column name where the grouping column is. Optional.
                  means = means, pvals = pvals,highlight_size = 1,
                  genes = c('TGFB1','FN1','PDGFB','VWF','LAMA2','PLAU','CD24'),
                  return_table = TRUE,col_option = viridis::viridis(50))
  
  
  library(RColorBrewer)
  display.brewer.all() 
  
#beauty
  df <- a[a$Var1 %in% AA,]
  
  BB <- c('AM1_SFRP5+ epi-AM1_CNN1+ stro',
          'AM1_IGFBP3+ stro-AM1_CNN1+ stro',
          'AM1_CNN1+ stro-AM1_CNN1+ stro',
          'AM1_Endothelial-AM1_CNN1+ stro',
          'AM1_Immune-AM1_CNN1+ stro',
          'AM1_ESR1+ SMC-AM1_CNN1+ stro',
          'AM1_PV-AM1_CNN1+ stro',
          'AM2_SFRP5+ epi-AM2_CNN1+ stro',
          'AM2_IGFBP3+ stro-AM2_CNN1+ stro',
          'AM2_CNN1+ stro-AM2_CNN1+ stro',
          'AM2_Endothelial-AM2_CNN1+ stro',
          'AM2_Immune-AM2_CNN1+ stro',
          'AM2_ESR1+ SMC-AM2_CNN1+ stro',
          'AM2_PV-AM2_CNN1+ stro',
          'AM3_SFRP5+ epi-AM3_CNN1+ stro',
          'AM3_IGFBP3+ stro-AM3_CNN1+ stro',
          'AM3_CNN1+ stro-AM3_CNN1+ stro',
          'AM3_Endothelial-AM3_CNN1+ stro',
          'AM3_Immune-AM3_CNN1+ stro',
          'AM3_ESR1+ SMC-AM3_CNN1+ stro',
          'AM3_PV-AM3_CNN1+ stro')
  
  df <- df[df$Var2 %in% BB,]
  
  df$Var1 <- factor(df$Var1,levels=rev(AA))
  df$Var2 <- factor(df$Var2,levels=BB)
  
  if ((length(standard_scale) > 0 && standard_scale) | 
      (length(scale) > 0 && scale) | (length(scale) < 
                                      1 && length(standard_scale) < 1)) {
    if (!(length(p.adjust.method) > 0 && p.adjust.method != 
        "none") ) {
      g <- ggplot(df, aes(x = Var2, y = Var1, color = -log10(pvals), 
                          fill = scaled_means, size = scaled_means))
    }
  }

  if (!is.null(highlight_size)) {
    g <- g + geom_point(pch = 21, na.rm = TRUE, stroke = highlight_size)
  }

  g <- g + theme_bw() + theme(axis.text.x = element_text(angle = 45, 
                                                         hjust = 0, color = "#000000"), axis.text.y = element_text(color = "#000000"), 
                              axis.ticks = element_blank(), axis.title.x = element_blank(), 
                              axis.title.y = element_blank()) + scale_x_discrete(position = "top") + 
    scale_color_gradientn(colors = highlight, na.value = "white") + 
    scale_radius(range = c(0, max_size))
  if (!noir) {
    if (length(col_option) != 1)  {
      g <- g + scale_fill_gradientn(colors = c("white", 
                                               (grDevices::colorRampPalette(c('#4393C3FF','#F7F7F7FF','#F4A582FF','#D6604DFF'))(100))), 
                                    na.value = "white")
    }
  }

 a <- g+theme(axis.text.x = element_text(angle = 90))+small_guide() 
 
