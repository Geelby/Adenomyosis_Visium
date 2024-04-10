#reference：https://hbctraining.github.io/scRNA-seq_online/schedule/links-to-lessons.html

library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(clustree)
library(cowplot)
library(tidyverse)
library(paletteer)
library(Matrix.utils)
library(edgeR)
library(Matrix)
library(reshape2)
library(S4Vectors)
library(SingleCellExperiment)
library(pheatmap)
library(apeglm)
library(png)
library(DESeq2)
library(RColorBrewer)
library(data.table)

#integrate single-cell data
  pro <- readRDS('/////')
  amsc <- readRDS('/////')
  
  dim(pro)
  dim(amsc)

  pro <- SCTransform(pro)
  amsc <- SCTransform(amsc)

  integ_features <- SelectIntegrationFeatures(object.list = c(pro,amsc),
                                                  nfeatures = 3000)
  # Prepare the SCT list object for integration
  split_seurat <- PrepSCTIntegration(object.list = c(pro,amsc),
                                         anchor.features = integ_features)
  # Find best buddies - can take a while to run
  integ_anchors <- FindIntegrationAnchors(object.list = split_seurat,
                                              normalization.method = "SCT",
                                              anchor.features = integ_features)
  # Integrate across conditions
  seurat_integrated <- IntegrateData(anchorset = integ_anchors,
                                         normalization.method = "SCT")
      
  table(seurat_integrated@meta.data$sample)

  #UMAP
  seurat_integrated <- RunPCA(object = seurat_integrated)

  # Plot the elbow plot
  ElbowPlot(object = seurat_integrated, 
                ndims = 40)
      
  # Run UMAP
  seurat_integrated <- RunUMAP(seurat_integrated,dims = 1:5,reduction = "pca")
      
  Idents(seurat_integrated) <- seurat_integrated$sample
  table(Idents(seurat_integrated))
      
  # Plot UMAP
  Idents(seurat_integrated) <- seurat_integrated$sample
      
  library(ggsci)
      
  DimPlot(seurat_integrated, reduction = "umap", group.by = "sample",
              label = F,label.box = F,repel = TRUE)+ scale_color_npg()+ 
        theme_dr(xlength = 0.2, 
                 ylength = 0.2,
                 arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
        theme(panel.grid = element_blank(),
              axis.title = element_text(face = 2,hjust = 0.03))
      
#findcluster

     DimHeatmap(seurat_integrated, 
                 dims = 1:20, 
                 cells = 500, 
                 balanced = TRUE)
     
     print(x = seurat_integrated[["pca"]], 
            dims = 1:10, 
            nfeatures = 5)
      
    ElbowPlot(object = seurat_integrated, 
                ndims = 40)
      
    # Determine the K-nearest neighbor graph
      seurat_integrated <- FindNeighbors(object = seurat_integrated, 
                                         dims = 1:5)
      
    # Determine the clusters for various resolutions  
      
      seurat_integrated <- FindClusters(object = seurat_integrated,
                                        resolution = c(0.01,0.05,0.1,0.2,0.3, 0.4, 0.5,0.6,0.7,0.8,0.85,0.9))
      
    # Explore resolutions
      
      res=c(0.01,0.05,0.1,0.2,0.3, 0.4, 0.5,0.6,0.7,0.8,0.85,0.9)
      
      DiffResolution_plot <- lapply(res, function(i){
        p <- DimPlot( seurat_integrated, reduction = "umap",label=T,group.by = paste("integrated_snn_res.",i,sep=""))+
          labs(title=paste("Resolution = ",i,sep=""))
        return(p)
      })
      p_dim <- cowplot::plot_grid(plotlist = DiffResolution_plot,ncol=3)
      p_dim

      # Assign identity of clusters
      seurat_integrated$seurat_clusters<- seurat_integrated$integrated_snn_res.0.7
      
      # Plot the UMAP
      Idents(seurat_integrated) <- "seurat_clusters"
      DimPlot(seurat_integrated,
              reduction = "umap",
              label = TRUE,
              label.size = 5
             )
     
    #Exploring known cell type markers
        DimPlot(object = seurat_integrated, 
                reduction = "umap", 
                label = TRUE) + NoLegend()
        
        # Select the RNA counts slot to be the default assay
        DefaultAssay(seurat_integrated) <- "RNA"
        
        # Normalize RNA data for visualization purposes
        seurat_integrated <- NormalizeData(seurat_integrated, verbose = FALSE)
        
        genes_to_check = c('DCN','PCOLCE','DKK1',
                           'PRL','IGFBP1',
                           'EPCAM','KRT8','WFDC2',
                           'SCGB2A2','SLC18A2','PAEP',
                           'FOXJ1','PIFO',
                           'CD3D','CD3E','CD3G',
                           'NKG7','NCAM1',
                           'ACTA2','RGS5',
                           'PECAM1','VWF',
                           'CD14','CD163','MS4A6A',
                           'CD19', 'CD79A', 'MS4A1',
                           'TPSAB1','TPSB2','PTPRC')
        p_markers <- DotPlot(seurat_integrated, features = genes_to_check,
                             assay='RNA'  )  + coord_flip()
        p_markers
        
        table(seurat_integrated@meta.data$sample)
        
        seurat_integrated$seurat_clusters <- seurat_integrated$integrated_snn_res.0.7
        Idents(seurat_integrated) <- seurat_integrated$seurat_clusters
        table(Idents(seurat_integrated))

        celltype=data.frame(ClusterID=0:22,
                            celltype= NA) 


        celltype[celltype$ClusterID %in% c(0,2,9,14,21),2]='Stromal' 
        celltype[celltype$ClusterID %in% c(6,12,13,18,20,22),2]='Epithelial' 
        celltype[celltype$ClusterID %in% c(4,8,10,16),2]='Endothelial' 
        celltype[celltype$ClusterID %in% c(7,15),2]='Immune' 
        celltype[celltype$ClusterID %in% c(1,5,19),2]='PV' 
        celltype[celltype$ClusterID %in% c(3,11),2]='SMC' 
        celltype[celltype$ClusterID %in% c(17),2]='unknown' 

        celltype
        table(celltype$celltype)
        seurat_integrated@meta.data$celltype = "NA"
        for(i in 1:nrow(celltype)){
          seurat_integrated@meta.data[which(seurat_integrated@meta.data$integrated_snn_res.0.7 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
        table(seurat_integrated@meta.data$celltype)
        
        seurat_integrated<- seurat_integrated[,!(seurat_integrated$celltype %in% c('unknown'))]
        
        seurat_integrated$celltype <- factor(seurat_integrated$celltype,levels=c('Epithelial','Stromal','Endothelial','Immune', 'SMC','PV'))
      
      #dotplot

        DefaultAssay(seurat_integrated) <- "RNA"
        th2=theme(axis.text=element_text(size=8,face="bold"),axis.title=element_blank(),legend.position="right",legend.text=element_text(size=5.5),legend.title=element_text(size=6,face="bold"),axis.text.x=element_text(angle = 90, hjust = 1, vjust = .5))
        
        seurat_integrated <- NormalizeData(seurat_integrated, verbose = FALSE)
        
        genes_to_check = rev(c(
                           'EPCAM','KRT8','WFDC2',
                           'IGF1','PCOLCE','SFRP1',
                           'PECAM1','VWF',
                           'PTPRC',
                           'ACTA2','MYH11',
                           'RGS5'))
        
        seurat_integrated$celltype <- factor(seurat_integrated$celltype,levels=rev(c('Epithelial','Stromal','Endothelial','Immune', 'SMC','PV')))
        
        p <- DotPlot(seurat_integrated, features = unique(genes_to_check),
                     assay='RNA' ,group.by = 'celltype' ,cols=c("white","#666666FF"),scale = T)   +th2 
        
        p

      #umap plot
        color_region = c('#C93F55FF', '#3C4B99FF', '#EACC62FF', '#469D76FF', '#DF9ED4FF','#924099FF')
        names(color_region)=c('Epithelial','Stromal','Endothelial','Immune', 'SMC','PV')
      
          library(tidydr)
          library(ggplot2)

          DimPlot(seurat_integrated, reduction = "umap", group.by = "celltype",cols=color_region,
                  label = F,label.box = F,repel = TRUE)+ 
            theme_dr(xlength = 0.2, 
                     ylength = 0.2,
                     arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
            theme(panel.grid = element_blank(),
                  axis.title = element_text(face = 2,hjust = 0.03))
      
      #Featureplot
         genes_to_check2 <- c("EPCAM")
        seurat_integrated <- ScaleData(seurat_integrated)
        p <- FeaturePlot(seurat_integrated,reduction = "umap",features = genes_to_check2, slot = "scale.data",
                         cols = c("#DCDCDC","#F8C0C8FF","#C93F55FF"),
                         combine = FALSE,pt.size=0.001)
        legend <- get_legend(p[[1]])
        for(i in 1:length(p) ){
          p[[i]] <- p[[i]] + NoLegend() + NoAxes()
        }
        p_list <- cowplot::plot_grid(plotlist = p,ncol = 1)
        p_list
        p_legend <- p_list+legend+plot_layout(widths = c(6, 1))
        p_legend
        ggsave("EPCAM_Featureplots.pdf",width=6,height = 6)
        
        genes_to_check2 <- c("PCOLCE")
        p <- FeaturePlot(seurat_integrated,reduction = "umap",features = genes_to_check2, slot = "scale.data",
                         cols = c("#DCDCDC","#A0B7D8FF","#3C4B99FF"),
                         combine = FALSE,pt.size=0.001)
        legend <- get_legend(p[[1]])
        for(i in 1:length(p) ){
          p[[i]] <- p[[i]] + NoLegend() + NoAxes()
        }
        p_list <- cowplot::plot_grid(plotlist = p,ncol = 1)
        p_list
        p_legend <- p_list+legend+plot_layout(widths = c(6, 1))
        p_legend
        ggsave("PCOLCE_Featureplots.pdf",width=6,height = 6)
        
        genes_to_check2 <- c("PTPRC")
        p <- FeaturePlot(seurat_integrated,reduction = "umap",features = genes_to_check2, slot = "scale.data",
                         cols = c("#DCDCDC","#CCDDAAFF","#469D76FF"),
                         combine = FALSE,pt.size=0.001)
        legend <- get_legend(p[[1]])
        for(i in 1:length(p) ){
          p[[i]] <- p[[i]] + NoLegend() + NoAxes()
        }
        p_list <- cowplot::plot_grid(plotlist = p,ncol = 1)
        p_list
        p_legend <- p_list+legend+plot_layout(widths = c(6, 1))
        p_legend
        ggsave("PTPRC_Featureplots.pdf",width=6,height = 6)
        
        genes_to_check2 <- c("PECAM1")
        p <- FeaturePlot(seurat_integrated,reduction = "umap",features = genes_to_check2, slot = "scale.data",
                         cols = c("#DCDCDC",'#F9E7C2FF','#EACC62FF'),
                         combine = FALSE,pt.size=0.001)
        legend <- get_legend(p[[1]])
        for(i in 1:length(p) ){
          p[[i]] <- p[[i]] + NoLegend() + NoAxes()
        }
        p_list <- cowplot::plot_grid(plotlist = p,ncol = 1)
        p_list
        p_legend <- p_list+legend+plot_layout(widths = c(6, 1))
        p_legend
        ggsave("PECAM1_Featureplots.pdf",width=6,height = 6)
        
        genes_to_check2 <- c("ACTA2")
        p <- FeaturePlot(seurat_integrated,reduction = "umap",features = genes_to_check2, slot = "scale.data",
                         cols = c("#DCDCDC",'#F5D2E6FF','#DF9ED4FF'),
                         combine = FALSE,pt.size=0.001)
        legend <- get_legend(p[[1]])
        for(i in 1:length(p) ){
          p[[i]] <- p[[i]] + NoLegend() + NoAxes()
        }
        p_list <- cowplot::plot_grid(plotlist = p,ncol = 1)
        p_list
        p_legend <- p_list+legend+plot_layout(widths = c(6, 1))
        p_legend
        ggsave("ACTA2_Featureplots.pdf",width=6,height = 6)
        
        genes_to_check2 <- c("RGS5")
        p <- FeaturePlot(seurat_integrated,reduction = "umap",features = genes_to_check2, slot = "scale.data",
                         cols = c("#DCDCDC",'#CABEE9FF','#924099FF'),
                         combine = FALSE,pt.size=0.001)
        legend <- get_legend(p[[1]])
        for(i in 1:length(p) ){
          p[[i]] <- p[[i]] + NoLegend() + NoAxes()
        }
        p_list <- cowplot::plot_grid(plotlist = p,ncol = 1)
        p_list
        p_legend <- p_list+legend+plot_layout(widths = c(6, 1))
        p_legend
        ggsave("RGS5_Featureplots.pdf",width=6,height = 6)

      #findmarker
        DefaultAssay(seurat_integrated)<- "RNA"  
        Idents(seurat_integrated) <- seurat_integrated$celltype
        table(Idents(seurat_integrated))
        levels(Idents(seurat_integrated))
        
        seurat.markers <- FindAllMarkers(object = seurat_integrated, only.pos = TRUE, 
                                         min.pct = 0.25, 
                                         thresh.use = 0.25)
        write.csv(seurat.markers,file='seurat.markers.csv')
        
        library(dplyr) 
        seurat_integrated <- ScaleData(seurat_integrated)
        top10 <- seurat.markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)
        p <- DoHeatmap(subset(seurat_integrated, downsample = 1000),assay = "RNA",
                       features = rev(unique(top10$gene)),size=3,
                       group.colors =color_region )+
          scale_fill_gradient2(
                                limits = c(-3,3),
                                 low = ("#4393C3FF"),
                               mid = "white",
                               high = ("#D6604DFF"),
                               midpoint = 0)
        p
        
#prepare data for cell2location
        library(SeuratDisk)
        file <- seurat_integrated
        main.loom <- as.loom(x = file, filename = "/////.loom", verbose = FALSE)
        write.csv(file@meta.data,'////.csv') 
        