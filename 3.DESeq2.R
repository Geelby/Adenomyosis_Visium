#reference：https://hbctraining.github.io/scRNA-seq_online/lessons/pseudobulk_DESeq2_scrnaseq.html

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

#epi
    setwd("figures/ST_epi_sub/forIPA/")  
    seurat_integrated <- readRDS('//////')
    seurat_integrated@meta.data$groupnew = "NA"
    
    seurat_integrated@meta.data[which(seurat_integrated@meta.data$slice == "AM1"),'groupnew'] <- "AM"
    seurat_integrated@meta.data[which(seurat_integrated@meta.data$slice == "AM2"),'groupnew'] <- "AM"
    seurat_integrated@meta.data[which(seurat_integrated@meta.data$slice == "AM3"),'groupnew'] <- "AM"
    seurat_integrated@meta.data[which(seurat_integrated@meta.data$slice == "UM"),'groupnew'] <- "UM"
    
    table(seurat_integrated@meta.data$groupnew)
    
    table(seurat_integrated$slice)
    seurat_integrated$slice <- factor(seurat_integrated$slice,levels=c('AM1','AM2','AM3','UM'))
    seurat_integrated$groupnew <- factor(seurat_integrated$groupnew,levels=c('UM','AM'))
    
    table(seurat_integrated$epicelltype)
      
      # Bring in Seurat object
      seurat <- seurat_integrated
      
      # Extract raw counts and metadata to create SingleCellExperiment object
      counts <- seurat@assays$Spatial@counts 
      
      metadata <- seurat@meta.data
      
      Idents(seurat) <- "epicelltype"
      table( Idents(seurat))
      # Set up metadata as desired for aggregation and DE analysis
      metadata$cluster_id <- factor(seurat@active.ident)
      
      # Create single cell experiment object
      sce <- SingleCellExperiment(assays = list(counts = counts), 
                                  colData = metadata)
      
      # Explore the raw counts for the dataset
      ## Check the assays present
      assays(sce)
      
      ## Check the counts matrix
      dim(counts(sce))
      counts(sce)[1:6, 1:6]
      
      # Explore the cellular metadata for the dataset
      
      dim(colData(sce))
      head(colData(sce))
      
      # Extract unique names of clusters (= levels of cluster_id factor variable)
      cluster_names <- levels(colData(sce)$cluster_id)
      cluster_names
      
      # Total number of clusters
      length(cluster_names)
      
      # Extract unique names of samples (= levels of sample_id factor variable)
      sample_names <- levels(colData(sce)$slice)
      sample_names
      
      # Total number of samples
      length(sample_names)
      
      # Subset metadata to include only the variables you want to aggregate across (here, we want to aggregate by sample and by cluster)
      groups <- colData(sce)[, c("cluster_id", "slice")]
      head(groups)
      
      # Aggregate across cluster-sample groups
      # transposing row/columns to have cell_ids as row names matching those of groups
      aggr_counts <- aggregate.Matrix(t(counts(sce)), 
                                      groupings = groups, fun = "sum") 
      
      # Explore output matrix
      class(aggr_counts)
      dim(aggr_counts)
      aggr_counts[1:6, 1:6]
      
      # Transpose aggregated matrix to have genes as rows and samples as columns
      aggr_counts <- t(aggr_counts)
      aggr_counts[1:6, 1:6]
      
      # Understanding tstrsplit()
      
      ## Exploring structure of function output (list)
      tstrsplit(colnames(aggr_counts), "_") %>% str()
      
      ## Comparing the first 10 elements of our input and output strings
      head(colnames(aggr_counts), n = 10)
      head(tstrsplit(colnames(aggr_counts), "_")[[1]], n = 10)
      
      # Using which() to look up tstrsplit() output
      oneM1_idx <- which(tstrsplit(colnames(aggr_counts), "_")[[1]] == "LDHB+ epi")
      oneM1_idx
      
      colnames(aggr_counts)[oneM1_idx]
      aggr_counts[1:10, oneM1_idx]
      
      # As a reminder, we stored our cell types in a vector called cluster_names
      cluster_names
      
      # Loop over all cell types to extract corresponding counts, and store information in a list
      
      ## Initiate empty list
      counts_ls <- list()
      
      for (i in 1:length(cluster_names)) {
        
        ## Extract indexes of columns in the global matrix that match a given cluster
        column_idx <- which(tstrsplit(colnames(aggr_counts), "_")[[1]] == cluster_names[i])
        
        ## Store corresponding sub-matrix as one element of a list
        counts_ls[[i]] <- aggr_counts[, column_idx]
        names(counts_ls)[i] <- cluster_names[i]
        
      }
      
      # Explore the different components of the list
      str(counts_ls)
      
      # Reminder: explore structure of metadata
      head(colData(sce))
      
      # Extract sample-level variables
      metadata <- colData(sce) %>% 
        as.data.frame() %>% 
        dplyr::select(slice, sampleid,groupnew)
      
      dim(metadata)
      head(metadata)
      
      # Exclude duplicated rows
      metadata <- metadata[!duplicated(metadata), ]
      
      dim(metadata)
      head(metadata)
      
      # Rename rows
      rownames(metadata) <- metadata$slice
      head(metadata)
      
      metadata$slice <- factor(metadata$slice,levels=c('AM1','AM2','AM3','UM'))
      metadata$groupnew <- factor(metadata$groupnew,levels=c('UM','AM'))
      
      levels(metadata$slice)
      str(metadata$slice)
      levels(metadata$groupnew)
      str(metadata$groupnew)
      
      # Number of cells per sample and cluster
      t <- table(colData(sce)$slice,
                 colData(sce)$cluster_id)
      t
      
      # Creating metadata list
      
      ## Initiate empty list
      metadata_ls <- list()
      
      for (i in 1:length(counts_ls)) {
        
        ## Initiate a data frame for cluster i with one row per sample (matching column names in the counts matrix)
        df <- data.frame(cluster_sample_id = colnames(counts_ls[[i]]))
        
        ## Use tstrsplit() to separate cluster (cell type) and sample IDs
        df$cluster_id <- tstrsplit(df$cluster_sample_id, "_")[[1]]
        df$sample_id  <- tstrsplit(df$cluster_sample_id, "_")[[2]]
        
        
        ## Retrieve cell count information for this cluster from global cell count table
        idx <- which(colnames(t) == unique(df$cluster_id))
        cell_counts <- t[, idx]
        
        ## Remove samples with zero cell contributing to the cluster
        cell_counts <- cell_counts[cell_counts > 0]
        
        ## Match order of cell_counts and sample_ids
        sample_order <- match(df$sample_id, names(cell_counts))
        cell_counts <- cell_counts[sample_order]
        
        ## Append cell_counts to data frame
        df$cell_count <- cell_counts
        
        rownames(df) <- df$sample_id
        
        ## Join data frame (capturing metadata specific to cluster) to generic metadata
        df <- merge(df, metadata, by.x='sample_id', by.y='slice')
        
        
        ## Update rownames of metadata to match colnames of count matrix, as needed later for DE
        rownames(df) <- df$cluster_sample_id
        
        ## Store complete metadata for cluster i in list
        metadata_ls[[i]] <- df
        names(metadata_ls)[i] <- unique(df$cluster_id)
        
      }
      
      # Explore the different components of the list
      str(metadata_ls)
      
      # Select cell type of interest
      cluster_names
      
      # Double-check that both lists have same names
      all(names(counts_ls) == names(metadata_ls))

      ######LDHB+ epi
        dir.create("LDHBepi")
        setwd("LDHBepi") 
        
        
        idx <- which(names(counts_ls) == "LDHB+ epi")
        cluster_counts <- counts_ls[[idx]]
        cluster_metadata <- metadata_ls[[idx]]
        
        # Check contents of extracted objects
        head(cluster_counts)
        head(cluster_metadata)
        
        # Check matching of matrix columns and metadata rows
        all(colnames(cluster_counts) == rownames(cluster_metadata))
        
        # Create DESeq2 object        
        dds <- DESeqDataSetFromMatrix(cluster_counts, 
                                      colData = cluster_metadata, 
                                      design = ~ groupnew)
        
        # Transform counts for data visualization
        rld <- rlog(dds, blind=TRUE)
        
        # Plot PCA
        DESeq2::plotPCA(rld, ntop = 500, intgroup = "groupnew")
        
        DESeq2::plotPCA(rld, ntop = 500, intgroup = "cell_count")
        
        DESeq2::plotPCA(rld, ntop = 500, intgroup = "sample_id")
        
        # Extract the rlog matrix from the object and compute pairwise correlation values
        rld_mat <- assay(rld)
        rld_cor <- cor(rld_mat)
        
        # Plot heatmap
        pheatmap(rld_cor, annotation = cluster_metadata[, c("groupnew","sample_id"), drop=F])
        
        # Run DESeq2 differential expression analysis
        dds <- DESeq(dds)
        
        # Plot dispersion estimates
        plotDispEsts(dds)
        
        # Check the coefficients for the comparison
        resultsNames(dds)
        
        # Generate results object
        res <- results(dds, 
                       name = "groupnew_AM_vs_UM",
                       alpha = 0.05)
        
        # Shrink the log2 fold changes to be more appropriate using the apeglm method - should cite [paper]() when using this method
        res <- lfcShrink(dds, 
                         coef = "groupnew_AM_vs_UM",
                         res=res,
                         type = "apeglm")
        
        # Turn the DESeq2 results object into a tibble for use with tidyverse functions
        res_tbl <- res %>%
          data.frame() %>%
          rownames_to_column(var = "gene") %>%
          as_tibble() %>%
          arrange(padj)
        
        # Check results output
        res_tbl 
        
        # Write all results to file
        write.csv(res_tbl,
                  paste0(unique(cluster_metadata$cluster_id), "_", 
                         levels(cluster_metadata$groupnew)[2], "_vs_", levels(cluster_metadata$groupnew)[1], "_all_genes.csv"),
                  quote = FALSE, 
                  row.names = FALSE)
        
        # Set thresholds
        padj_cutoff <- 0.05
        
        # Subset the significant results
        sig_res <- dplyr::filter(res_tbl, padj < padj_cutoff) %>%
          dplyr::arrange(padj)
        
        # Check significant genes output
        sig_res
        
        # Write significant results to file
        write.csv(sig_res,
                  paste0(unique(cluster_metadata$cluster_id), "_", 
                         levels(cluster_metadata$groupnew)[2], "_vs_", levels(cluster_metadata$groupnew)[1], "_signif_genes.csv"),
                  quote = FALSE, 
                  row.names = FALSE)
        
#SMC
        setwd("figures/ST_SMC_sub/Deseq2/")  
        
        seurat_integrated <- readRDS('//////')
        seurat_integrated@meta.data$groupnew = "NA"
        
        seurat_integrated@meta.data[which(seurat_integrated@meta.data$slice == "AM1"),'groupnew'] <- "AM_EU"
        seurat_integrated@meta.data[which(seurat_integrated@meta.data$slice == "AM2"),'groupnew'] <- "AM_EC"
        seurat_integrated@meta.data[which(seurat_integrated@meta.data$slice == "AM3"),'groupnew'] <- "AM_EC"
        seurat_integrated@meta.data[which(seurat_integrated@meta.data$slice == "UM"),'groupnew'] <- "UM_EU"
        
        table(seurat_integrated@meta.data$SMCcelltype)
        
        seurat_integrated$slice <- factor(seurat_integrated$slice,levels=c('AM1','AM2','AM3','UM'))
        seurat_integrated$groupnew <- factor(seurat_integrated$groupnew,levels=c("UM_EU","AM_EU","AM_EC"))
        
        table(seurat_integrated$SMCcelltype)
        table(seurat_integrated$SMCcelltype)
        
        # Bring in Seurat object
        seurat <- seurat_integrated
        
        # Extract raw counts and metadata to create SingleCellExperiment object
        counts <- seurat@assays$Spatial@counts 
        
        metadata <- seurat@meta.data
        
        Idents(seurat) <- "SMCcelltype"
        table( Idents(seurat))
        # Set up metadata as desired for aggregation and DE analysis
        metadata$cluster_id <- factor(seurat@active.ident)
        
        # Create single cell experiment object
        sce <- SingleCellExperiment(assays = list(counts = counts), 
                                    colData = metadata)
        
        # Explore the raw counts for the dataset
        ## Check the assays present
        assays(sce)
        
        ## Check the counts matrix
        dim(counts(sce))
        counts(sce)[1:6, 1:6]
        
        # Explore the cellular metadata for the dataset
        
        dim(colData(sce))
        head(colData(sce))
        
        # Extract unique names of clusters (= levels of cluster_id factor variable)
        cluster_names <- levels(colData(sce)$cluster_id)
        cluster_names
        
        # Total number of clusters
        length(cluster_names)
        
        # Extract unique names of samples (= levels of sample_id factor variable)
        sample_names <- levels(colData(sce)$slice)
        sample_names
        
        # Total number of samples
        length(sample_names)
        
        # Subset metadata to include only the variables you want to aggregate across (here, we want to aggregate by sample and by cluster)
        groups <- colData(sce)[, c("cluster_id", "slice")]
        head(groups)
        
        # Aggregate across cluster-sample groups
        # transposing row/columns to have cell_ids as row names matching those of groups
        aggr_counts <- aggregate.Matrix(t(counts(sce)), 
                                        groupings = groups, fun = "sum") 
        
        # Explore output matrix
        class(aggr_counts)
        dim(aggr_counts)
        aggr_counts[1:6, 1:6]
        
        # Transpose aggregated matrix to have genes as rows and samples as columns
        aggr_counts <- t(aggr_counts)
        aggr_counts[1:6, 1:6]
        
        # Understanding tstrsplit()
        
        ## Exploring structure of function output (list)
        tstrsplit(colnames(aggr_counts), "_") %>% str()
        
        ## Comparing the first 10 elements of our input and output strings
        head(colnames(aggr_counts), n = 10)
        head(tstrsplit(colnames(aggr_counts), "_")[[1]], n = 10)
        
        # Using which() to look up tstrsplit() output
        oneM1_idx <- which(tstrsplit(colnames(aggr_counts), "_")[[1]] == "ESR1+ SMC")
        oneM1_idx
        
        colnames(aggr_counts)[oneM1_idx]
        aggr_counts[1:10, oneM1_idx]
        
        # As a reminder, we stored our cell types in a vector called cluster_names
        cluster_names
        
        # Loop over all cell types to extract corresponding counts, and store information in a list
        
        ## Initiate empty list
        counts_ls <- list()
        
        for (i in 1:length(cluster_names)) {
          
          ## Extract indexes of columns in the global matrix that match a given cluster
          column_idx <- which(tstrsplit(colnames(aggr_counts), "_")[[1]] == cluster_names[i])
          
          ## Store corresponding sub-matrix as one element of a list
          counts_ls[[i]] <- aggr_counts[, column_idx]
          names(counts_ls)[i] <- cluster_names[i]
          
        }
        
        # Explore the different components of the list
        str(counts_ls)
        
        # Reminder: explore structure of metadata
        head(colData(sce))
        
        # Extract sample-level variables
        metadata <- colData(sce) %>% 
          as.data.frame() %>% 
          dplyr::select(slice, sampleid,groupnew)
        
        dim(metadata)
        head(metadata)
        
        # Exclude duplicated rows
        metadata <- metadata[!duplicated(metadata), ]
        
        dim(metadata)
        head(metadata)
        
        # Rename rows
        rownames(metadata) <- metadata$slice
        head(metadata)
        
        metadata$slice <- factor(metadata$slice,levels=c('UM','AM1','AM2','AM3'))
        metadata$groupnew <- factor(metadata$groupnew,levels=c("UM_EU","AM_EU","AM_EC"))
        
        levels(metadata$slice)
        str(metadata$slice)
        levels(metadata$groupnew)
        str(metadata$groupnew)
        
        # Number of cells per sample and cluster
        t <- table(colData(sce)$slice,
                   colData(sce)$cluster_id)
        t
        
        # Creating metadata list
        
        ## Initiate empty list
        metadata_ls <- list()
        
        for (i in 1:length(counts_ls)) {
          
          ## Initiate a data frame for cluster i with one row per sample (matching column names in the counts matrix)
          df <- data.frame(cluster_sample_id = colnames(counts_ls[[i]]))
          
          ## Use tstrsplit() to separate cluster (cell type) and sample IDs
          df$cluster_id <- tstrsplit(df$cluster_sample_id, "_")[[1]]
          df$sample_id  <- tstrsplit(df$cluster_sample_id, "_")[[2]]
          
          
          ## Retrieve cell count information for this cluster from global cell count table
          idx <- which(colnames(t) == unique(df$cluster_id))
          cell_counts <- t[, idx]
          
          ## Remove samples with zero cell contributing to the cluster
          cell_counts <- cell_counts[cell_counts > 0]
          
          ## Match order of cell_counts and sample_ids
          sample_order <- match(df$sample_id, names(cell_counts))
          cell_counts <- cell_counts[sample_order]
          
          ## Append cell_counts to data frame
          df$cell_count <- cell_counts
          
          rownames(df) <- df$sample_id
          
          ## Join data frame (capturing metadata specific to cluster) to generic metadata
          df <- merge(df, metadata, by.x='sample_id', by.y='slice')
          
          
          ## Update rownames of metadata to match colnames of count matrix, as needed later for DE
          rownames(df) <- df$cluster_sample_id
          
          ## Store complete metadata for cluster i in list
          metadata_ls[[i]] <- df
          names(metadata_ls)[i] <- unique(df$cluster_id)
          
        }
        
        # Explore the different components of the list
        str(metadata_ls)
        
        # Select cell type of interest
        cluster_names
        
        # Double-check that both lists have same names
        all(names(counts_ls) == names(metadata_ls))
        
        # Load DEGreport
        library(DEGreport)
        
        clustx <- "ESR1+ SMC"
        
        # Extract counts matrix and metadata for cluster x
        idx <- which(names(counts_ls) == clustx)
        cluster_counts <- counts_ls[[idx]]
        cluster_metadata <- metadata_ls[[idx]]
        
        # Print error message if sample names do not match
        all(colnames(cluster_counts) == rownames(cluster_metadata)) 
        
        # Run DESeq2
        dds <- DESeqDataSetFromMatrix(cluster_counts, 
                                      colData = cluster_metadata, 
                                      design = ~ groupnew)
        dds_lrt <- DESeq(dds, test = "LRT", reduced = ~ 1)
        
        # Extract results
        res_LRT <- results(dds_lrt)
        
        # Create a tibble for LRT results
        res_LRT_tb <- res_LRT %>%
          data.frame() %>%
          rownames_to_column(var = "gene") %>% 
          as_tibble()
        
        # Save all results
        write.csv(res_LRT_tb,
                  paste0( clustx, "_LRT_all_genes.csv"),
                  quote = FALSE, 
                  row.names = FALSE)
        
        # Subset to return genes with padj < 0.05
        sigLRT_genes <- res_LRT_tb %>% 
          filter(padj < 0.05)
        
        # Save significant results
        write.csv(sigLRT_genes,
                  paste0(clustx, "_LRT_signif_genes.csv"),
                  quote = FALSE, 
                  row.names = FALSE)
        
        # Transform counts for data visualization
        rld <- rlog(dds_lrt, blind = TRUE)
        
        # Extract the rlog matrix from the object and compute pairwise correlation values
        rld_mat <- assay(rld)
        rld_cor <- cor(rld_mat)
        
        # Obtain rlog values for those significant genes
        cluster_rlog <- rld_mat[sigLRT_genes$gene, ]
        cluster_meta_sig <- cluster_metadata[which(rownames(cluster_metadata) %in% colnames(cluster_rlog)), ]
        
        cluster_groups <- degPatterns(cluster_rlog, metadata = cluster_meta_sig,
                                      time = "groupnew", col = "cluster_id", plot = F)
        library(ggthemes)
        
        p <- degPlotCluster(cluster_groups$normalized, "groupnew", "cluster_id",boxes = F,points = F,smooth = F)
        p
        p+ theme_tufte()+theme(axis.text.x = element_text(angle = 45,vjust = 0.85,hjust = 0.75))+ scale_color_manual(values=c('#E8D8E8FF'))
        
        # Use the `degPatterns` function from DEGreport package to show gene clusters across sample groups
        cluster_groups <- degPatterns(cluster_rlog, metadata = cluster_meta_sig,
                                      time = "groupnew", col = NULL)
        
        # Save what is stored in the `df` component
        write.csv(cluster_groups$df,
                  paste0(clustx, "_LRT_DEgene_groups.csv"),
                  quote = FALSE, 
                  row.names = FALSE)
        
        #DEG
        seurat.markers <- read_csv("ESR1+ SMC_LRT_DEgene_groups.csv")
        sce.markers <- seurat.markers
        library(ggpubr)
        library(clusterProfiler)
        options(stringsAsFactors = F)
        suppressMessages(library('org.Hs.eg.db'))
        keytypes(org.Hs.eg.db)
        
        go_organism <- "org.Hs.eg.db" 
        kegg_organism <- "hsa"
        
        ids = bitr(sce.markers$genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=go_organism)
        sce.markers = merge(sce.markers, ids, by.x='genes', by.y='SYMBOL')
        
        #Adding Gene Annotations:load in an annotation file
        annotations <- read.csv("///////annotation.csv")
        
        # Combine markers with gene descriptions 
        sce.markers <- sce.markers %>% 
          left_join(y = unique(annotations[, c("gene_name", "description")]),
                    by = c("genes" = "gene_name"))
        
        CSV <- data.frame(sce.markers)
        write.csv(CSV,'ESR1SMC_makers_annotation.csv')  
        
        #GO
        gcSample=split(sce.markers$ENTREZID, sce.markers$cluster ) 
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
          y=paste(y,collapse = "/")})
        
        datt$geneID=symbol
        write.csv(datt,'enrichGOBP.csv') 
        
        highcol='#DF9ED4FF'
        lowcol='#F3E2ECFF'
        
        p <- dotplot(xx,showCategory = c('response to estrogen','collagen fibril organization','collagen metabolic process',
                                         'extracellular matrix organization','extracellular structure organization',
                                         'negative regulation of locomotion','negative regulation of cell migration','negative regulation of peptidase activity','response to hypoxia'))
        p + theme(axis.text.x = element_text(
          angle = 45,
          vjust = 0.5, hjust = 0.5
        ))+ scale_color_continuous(low=highcol, high=lowcol)+ scale_y_discrete(labels = function(x) str_wrap(x, width = 60) )+theme_few()
        
        
        
        
        
        
        
        
        


