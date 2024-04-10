#reference：https://geneswitches.ddnetbio.com/

library(GeneSwitches)
library(SingleCellExperiment)

dir.create("geneswitch")
setwd("geneswitch")

    Subset <- readr::read_rds('//////episub.rds')
    
    Idents(Subset)="slice"
    table(Idents(Subset))
    
    Idents(Subset) <- Subset$epicelltype
    table(Idents(Subset))

    SpatialDimPlot(Subset, label = F, label.size = 3,
                   images = c('UM','AM1','AM2','AM3'), ncol = 3,pt.size.factor = 1.2,cols =my_color )+plot_layout(guides = 'collect')
    
    
    sce <- Subset
    DefaultAssay(sce) <- "Spatial"
    levels(Idents(sce))
    
    sce
    
    sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
    
    logexpdata <- sce@assays$Spatial@data
    logexpdata <- as(logexpdata, 'matrix')
    
  dim(logexpdata)
  logexpdata[0:5,0:5]
  
  cds <- my_cds_subset
  
  plot_monocle_State(cds)
  
  plot_cell_trajectory(cds, color_by = "Pseudotime") 
  plot_cell_trajectory(cds, color_by = "slice") + facet_wrap("~slice",nrow=1)
  
  sce_p1 <- convert_monocle2(monocle2_obj = cds, 
                             states = c(3,2), expdata = logexpdata)
  sce_p2 <- convert_monocle2(monocle2_obj = cds, 
                             states = c(3,1), expdata = logexpdata)
  
  #sce_p2_AM
    
    ### check the threshold for binarization
    h <- hist(assays(sce_p2)$expdata, breaks = 200, plot = FALSE)
    {plot(h, freq = FALSE, xlim = c(0,2), ylim = c(0,1), main = "Histogram of gene expression",
          xlab = "Gene expression", col = "darkgoldenrod2", border = "grey")
      abline(v=0.2, col="blue")}
    
    ###In this example, we choose 0.2 (blue line, also set as default) as the threshold.
    bn_cutoff <- 0.2
    sce_p2 <- binarize_exp(sce_p2, fix_cutoff = TRUE, binarize_cutoff = bn_cutoff, ncores = 4)
    
    ## fit logistic regression and find the switching pseudo-time point for each gene
    ## with downsampling. This step takes less than 1 mins
    sce_p2 <- find_switch_logistic_fastglm(sce_p2, downsample = TRUE, show_warning = FALSE)
    
    ## filter top best fitting switching genes among all the genes
    sg_allgenes <- filter_switchgenes(sce_p2, allgenes = TRUE, topnum =10000)
    
    ## combine switching genes and remove duplicated genes from sg_allgenes
    sg_vis <- sg_allgenes
    
    my_color = c('#000000')
    
    names(my_color)=c('All genes')
    
    plot_timeline_ggplot(sg_vis, timedata = sce_p2$Pseudotime, txtsize = 3)+
      geom_hline(yintercept = 0.5,lty=3,col="black",lwd=0.5)+
      scale_color_manual(values = my_color)
    
    #beauty
      tml=sg_vis
      timedata=sce_p2$Pseudotime
      iffulltml = TRUE
      txtsize = 3
      color_by = "feature_type"
      
      tml <- as.data.frame(tml)
      tml <- tml[order(tml$switch_at_time), ]
      tml$direction_num <- -1
      if ("up" %in% tml$direction) {
        tml[tml$direction == "up", ]$direction_num <- 1
      }
      tml$color_by <- as.factor(tml[, color_by])
      tml$feature_name <- rownames(tml)
      head(tml)
      if (iffulltml) {
        pseudotime_step <- (max(timedata) - min(timedata))/4
        pseudotime_range <- seq(min(timedata), max(timedata), 
                                by = pseudotime_step)
        pseudotime_df <- data.frame(pseudotime_range, pseudotime_format = round(pseudotime_range, 
                                                                                1))
      }
      
      for_label <- tml[tml$geneID %in% c('IHH','LINC01480','OVGP1','CA12','MMP11','FOXJ1','CLU','HSPA1A','MTRNRAL1','SOX17','SFRP1'),]
      
      tml_plot <- ggplot(tml, aes(x = switch_at_time, y = pseudoR2s * 
                                    direction_num, col = color_by, label = feature_name)) + 
        geom_point(size = txtsize/8) + xlab("Pseudo-timeline") + 
        ylab("Quality of fitting (R^2)") + theme_classic()+
        geom_hline(yintercept = 0, color = "black", 
                   size = 0.6)+
        geom_label(data = pseudotime_df, aes(x = pseudotime_range, y = 0, label = pseudotime_format), size = (txtsize - 0.5), color = "black")+
        geom_hline(yintercept = 0.5,lty=3,col="black",lwd=0.5)+
        scale_color_manual(values = my_color)+
        geom_point(size = 2, shape = 1, data = for_label) +
        ggrepel::geom_label_repel(
          aes(label = geneID),
          data = for_label,
          color="black"
        )
      
      tml_plot
    
    