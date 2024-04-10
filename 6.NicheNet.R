#reference：https://github.com/saeyslab/nichenetr/blob/master/vignettes/differential_nichenet.md

library(nichenetr)
library(RColorBrewer)
library(tidyverse)
library(Seurat) 

  seurat_oobj= readr::read_rds('data/ST_new_invade_env_UMenv.rds')
  
  sooo <- seurat_oobj
  
  sooo
  
  table(sooo$slice)
  table(Idents(sooo))
  
  sooo$annotation <- paste0(sooo$slice, '_', sooo$celltypenewaddsub)
  
  Idents(sooo) = sooo$celltypenewaddsub
  
  table(sooo$celltypenewaddsub)
  table(Idents(sooo))
  
  seurat_obj <- sooo
  
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
  
  p2 <- SpatialDimPlot(seurat_obj, label = F, label.size = 3,
                       images = c("UM","AM1","AM2","AM3"), ncol = 2,crop = F,cols =my_color )
  p2

#Differential NicheNet analysis between niches of interest
  seurat_obj <- sooo
  DimPlot(seurat_obj, group.by = "annotation", label = TRUE) 
  seurat_obj = SetIdent(seurat_obj, value = "annotation")
  table(Idents(seurat_obj))
  
  #Read in the NicheNet ligand-receptor network and ligand-target matrix
  ligand_target_matrix = readRDS("datanichenetr/ligand_target_matrix.rds")
  ligand_target_matrix[1:5,1:5] # target genes in rows, ligands in columns
  
  lr_network = readRDS("datanichenetr/lr_network.rds")
  lr_network = lr_network %>% mutate(bonafide = ! database %in% c("ppi_prediction","ppi_prediction_go"))
  lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% distinct(ligand, receptor, bonafide)
  
  head(lr_network)
  
  organism = "human"
  
  ##1. Define the niches/microenvironments of interest
  niches = list(
    "UM_niche" = list(
      "sender" = c("UM_SFRP5+ epi","UM_Endothelial","UM_DES+ SMC","UM_PV","UM_Immune","UM_CNN1+ stro","UM_IGFBP3+ stro",'UM_ESR1+ SMC'),
      "receiver" = c("UM_SFRP5+ epi")),
    "AM1_niche" = list(
      "sender" = c("AM1_SFRP5+ epi","AM1_Endothelial","AM1_IGFBP3+ stro","AM1_CNN1+ stro","AM1_ESR1+ SMC","AM1_Immune",'AM1_DES+ SMC','AM1_PV'),
      "receiver" = c("AM1_SFRP5+ epi")))

  ##2.Calculate differential expression between the niches
  assay_oi = "SCT" # other possibilities: RNA,...
 # seurat_obj = PrepSCTFindMarkers(seurat_obj, assay = "SCT", verbose = TRUE)
  
  DE_sender = calculate_niche_de(seurat_obj = seurat_obj %>% subset(features = lr_network$ligand %>% intersect(rownames(seurat_obj))), niches = niches, type = "sender", assay_oi = assay_oi) # only ligands important for sender cell types
  DE_receiver = calculate_niche_de(seurat_obj = seurat_obj %>% subset(features = lr_network$receptor %>% unique()), niches = niches, type = "receiver", assay_oi = assay_oi) # only receptors now, later on: DE analysis to find targets
  
  DE_sender = DE_sender %>% mutate(avg_log2FC = ifelse(avg_log2FC == Inf, max(avg_log2FC[is.finite(avg_log2FC)]), ifelse(avg_log2FC == -Inf, min(avg_log2FC[is.finite(avg_log2FC)]), avg_log2FC)))
  DE_receiver = DE_receiver %>% mutate(avg_log2FC = ifelse(avg_log2FC == Inf, max(avg_log2FC[is.finite(avg_log2FC)]), ifelse(avg_log2FC == -Inf, min(avg_log2FC[is.finite(avg_log2FC)]), avg_log2FC)))
  
 # Process DE results:
  expression_pct = 0.10
  DE_sender_processed = process_niche_de(DE_table = DE_sender, niches = niches, expression_pct = expression_pct, type = "sender")
  DE_receiver_processed = process_niche_de(DE_table = DE_receiver, niches = niches, expression_pct = expression_pct, type = "receiver")
  
  #Combine sender-receiver DE based on L-R pairs:
  specificity_score_LR_pairs = "min_lfc"
  DE_sender_receiver = combine_sender_receiver_de(DE_sender_processed, DE_receiver_processed, lr_network, specificity_score = specificity_score_LR_pairs)
  
  ##3. Optional: Calculate differential expression between the different spatial regions
  ##4. Calculate ligand activities and infer active ligand-target links
  
  lfc_cutoff = 0.15 # recommended for 10x as min_lfc cutoff. 
  specificity_score_targets = "min_lfc"
  
  DE_receiver_targets = calculate_niche_de_targets(seurat_obj = seurat_obj, niches = niches, lfc_cutoff = lfc_cutoff, expression_pct = expression_pct, assay_oi = assay_oi) 
  DE_receiver_processed_targets = process_receiver_target_de(DE_receiver_targets = DE_receiver_targets, niches = niches, expression_pct = expression_pct, specificity_score = specificity_score_targets)
  
  background = DE_receiver_processed_targets  %>% pull(target) %>% unique()
  geneset_UM = DE_receiver_processed_targets %>% filter(receiver == niches$UM_niche$receiver & target_score >= lfc_cutoff & target_significant == 1 & target_present == 1) %>% pull(target) %>% unique()
  geneset_AM1 = DE_receiver_processed_targets %>% filter(receiver == niches$AM1_niche$receiver & target_score >= lfc_cutoff & target_significant == 1 & target_present == 1) %>% pull(target) %>% unique()
  
  geneset_UM %>% setdiff(rownames(ligand_target_matrix))
  geneset_AM1 %>% setdiff(rownames(ligand_target_matrix))
  
  length(geneset_UM)
  length(geneset_AM1)
  
  top_n_target = 250
  
  niche_geneset_list = list(
    "UM_niche" = list(
      "receiver" = "UM_SFRP5+ epi",
      "geneset" = geneset_UM,
      "background" = background),
    "AM1_niche" = list(
      "receiver" = "AM1_SFRP5+ epi",
      "geneset" = geneset_AM1 ,
      "background" = background) )
  
  ligand_activities_targets = get_ligand_activities_targets(niche_geneset_list = niche_geneset_list, ligand_target_matrix = ligand_target_matrix, top_n_target = top_n_target)
  
  
  ##5. Calculate (scaled) expression of ligands, receptors and targets across cell types of interest (log expression values and expression fractions)
  features_oi = union(lr_network$ligand, lr_network$receptor) %>% union(ligand_activities_targets$target) %>% setdiff(NA)
  
  dotplot = suppressWarnings(Seurat::DotPlot(seurat_obj %>% subset(idents = niches %>% unlist() %>% unique()), features = features_oi, assay = assay_oi))
  exprs_tbl = dotplot$data %>% as_tibble()
  exprs_tbl = exprs_tbl %>% rename(celltype = id, gene = features.plot, expression = avg.exp, expression_scaled = avg.exp.scaled, fraction = pct.exp) %>%
    mutate(fraction = fraction/100) %>% as_tibble() %>% select(celltype, gene, expression, expression_scaled, fraction) %>% distinct() %>% arrange(gene) %>% mutate(gene = as.character(gene))
  
  exprs_tbl_ligand = exprs_tbl %>% filter(gene %in% lr_network$ligand) %>% rename(sender = celltype, ligand = gene, ligand_expression = expression, ligand_expression_scaled = expression_scaled, ligand_fraction = fraction) 
  exprs_tbl_receptor = exprs_tbl %>% filter(gene %in% lr_network$receptor) %>% rename(receiver = celltype, receptor = gene, receptor_expression = expression, receptor_expression_scaled = expression_scaled, receptor_fraction = fraction)
  exprs_tbl_target = exprs_tbl %>% filter(gene %in% ligand_activities_targets$target) %>% rename(receiver = celltype, target = gene, target_expression = expression, target_expression_scaled = expression_scaled, target_fraction = fraction)
  
  exprs_tbl_ligand = exprs_tbl_ligand %>%  mutate(scaled_ligand_expression_scaled = scale_quantile_adapted(ligand_expression_scaled)) %>% mutate(ligand_fraction_adapted = ligand_fraction) %>% mutate_cond(ligand_fraction >= expression_pct, ligand_fraction_adapted = expression_pct)  %>% mutate(scaled_ligand_fraction_adapted = scale_quantile_adapted(ligand_fraction_adapted))
  exprs_tbl_receptor = exprs_tbl_receptor %>% mutate(scaled_receptor_expression_scaled = scale_quantile_adapted(receptor_expression_scaled))  %>% mutate(receptor_fraction_adapted = receptor_fraction) %>% mutate_cond(receptor_fraction >= expression_pct, receptor_fraction_adapted = expression_pct)  %>% mutate(scaled_receptor_fraction_adapted = scale_quantile_adapted(receptor_fraction_adapted))

  ##6. Expression fraction and receptor
  exprs_sender_receiver = lr_network %>% 
    inner_join(exprs_tbl_ligand, by = c("ligand")) %>% 
    inner_join(exprs_tbl_receptor, by = c("receptor")) %>% inner_join(DE_sender_receiver %>% distinct(niche, sender, receiver))
  
  ligand_scaled_receptor_expression_fraction_df = exprs_sender_receiver %>% group_by(ligand, receiver) %>% mutate(rank_receptor_expression = dense_rank(receptor_expression), rank_receptor_fraction  = dense_rank(receptor_fraction)) %>% mutate(ligand_scaled_receptor_expression_fraction = 0.5*( (rank_receptor_fraction / max(rank_receptor_fraction)) + ((rank_receptor_expression / max(rank_receptor_expression))) ) )  %>% distinct(ligand, receptor, receiver, ligand_scaled_receptor_expression_fraction, bonafide) %>% distinct() %>% ungroup() 
  
  ##7. Prioritization of ligand-receptor and ligand-target links
  prioritizing_weights = c("scaled_ligand_score" = 5,
                           "scaled_ligand_expression_scaled" = 1,
                           "ligand_fraction" = 1,
                           "scaled_ligand_score_spatial" = 0, 
                           "scaled_receptor_score" = 0.5,
                           "scaled_receptor_expression_scaled" = 0.5,
                           "receptor_fraction" = 1, 
                           "ligand_scaled_receptor_expression_fraction" = 1,
                           "scaled_receptor_score_spatial" = 0,
                           "scaled_activity" = 0,
                           "scaled_activity_normalized" = 1,
                           "bona_fide" = 1)
  
  output = list(DE_sender_receiver = DE_sender_receiver, ligand_scaled_receptor_expression_fraction_df = ligand_scaled_receptor_expression_fraction_df, sender_spatial_DE_processed = sender_spatial_DE_processed, receiver_spatial_DE_processed = receiver_spatial_DE_processed,
                ligand_activities_targets = ligand_activities_targets, DE_receiver_processed_targets = DE_receiver_processed_targets, exprs_tbl_ligand = exprs_tbl_ligand,  exprs_tbl_receptor = exprs_tbl_receptor, exprs_tbl_target = exprs_tbl_target)
  prioritization_tables = get_prioritization_tables(output, prioritizing_weights)
  
  prioritization_tables$prioritization_tbl_ligand_receptor %>% filter(receiver == niches[[1]]$receiver) %>% head(10)
  prioritization_tables$prioritization_tbl_ligand_target %>% filter(receiver == niches[[1]]$receiver) %>% head(10)
  
  prioritization_tables$prioritization_tbl_ligand_receptor %>% filter(receiver == niches[[2]]$receiver) %>% head(10)
  prioritization_tables$prioritization_tbl_ligand_target %>% filter(receiver == niches[[2]]$receiver) %>% head(10)
  
  prioritization_tables$prioritization_tbl_ligand_receptor = prioritization_tables$prioritization_tbl_ligand_receptor %>% mutate(receiver = factor(receiver, levels = c("UM_SFRP5+ epi","AM1_SFRP5+ epi")), niche = factor(niche, levels = c("UM_niche","AM1_niche"))) 
  prioritization_tables$prioritization_tbl_ligand_target = prioritization_tables$prioritization_tbl_ligand_target %>% mutate(receiver = factor(receiver, levels = c("UM_SFRP5+ epi","AM1_SFRP5+ epi")), niche = factor(niche, levels = c("UM_niche","AM1_niche"))) 
  
  ##8. Visualization of the Differential NicheNet output
  top_ligand_niche_df = prioritization_tables$prioritization_tbl_ligand_receptor %>% select(niche, sender, receiver, ligand, receptor, prioritization_score) %>% group_by(ligand) %>% top_n(1, prioritization_score) %>% ungroup() %>% select(ligand, receptor, niche) %>% rename(top_niche = niche)
  top_ligand_receptor_niche_df = prioritization_tables$prioritization_tbl_ligand_receptor %>% select(niche, sender, receiver, ligand, receptor, prioritization_score) %>% group_by(ligand, receptor) %>% top_n(1, prioritization_score) %>% ungroup() %>% select(ligand, receptor, niche) %>% rename(top_niche = niche)
  
  ligand_prioritized_tbl_oi = prioritization_tables$prioritization_tbl_ligand_receptor %>% select(niche, sender, receiver, ligand, prioritization_score) %>% group_by(ligand, niche) %>% top_n(1, prioritization_score) %>% ungroup() %>% distinct() %>% inner_join(top_ligand_niche_df) %>% filter(niche == top_niche) %>% group_by(niche) %>% top_n(50, prioritization_score) %>% ungroup() # get the top50 ligands per niche

  write.csv(ligand_prioritized_tbl_oi,'/SFRP5epi/receiver_ligand_prioritized_tbl_oi.csv') 
  write.csv(prioritization_tables$prioritization_tbl_ligand_receptor,'/SFRP5epi/receiver_prioritization_tbl_ligand_receptor.csv') 
  write.csv(prioritization_tables$prioritization_tbl_ligand_target,'/SFRP5epi/receiver_prioritization_tbl_ligand_target.csv') 
  
##
  receiver_oi = "AM1_SFRP5+ epi" 
  
  filtered_ligands = ligand_prioritized_tbl_oi %>% filter(receiver == receiver_oi) %>% pull(ligand) %>% unique()
  
  prioritized_tbl_oi = prioritization_tables$prioritization_tbl_ligand_receptor %>% filter(ligand %in% filtered_ligands) %>% select(niche, sender, receiver, ligand,  receptor, ligand_receptor, prioritization_score) %>% distinct() %>% inner_join(top_ligand_receptor_niche_df) %>% group_by(ligand) %>% filter(receiver == receiver_oi) %>% top_n(2, prioritization_score) %>% ungroup() 
  
  lfc_plot = make_ligand_receptor_lfc_plot(receiver_oi, prioritized_tbl_oi, prioritization_tables$prioritization_tbl_ligand_receptor, plot_legend = FALSE, heights = NULL, widths = NULL)
  lfc_plot
  
  ##Ligand expression, activity and target genes
  exprs_activity_target_plot = make_ligand_activity_target_exprs_plot(receiver_oi, prioritized_tbl_oi,  prioritization_tables$prioritization_tbl_ligand_receptor,  prioritization_tables$prioritization_tbl_ligand_target, output$exprs_tbl_ligand,  output$exprs_tbl_target, lfc_cutoff, ligand_target_matrix, plot_legend = FALSE, heights = NULL, widths = NULL)
  p <- exprs_activity_target_plot$combined_plot
  p
  
  #top20
  filtered_ligands = ligand_prioritized_tbl_oi %>% filter(receiver == receiver_oi) %>% top_n(20, prioritization_score) %>% pull(ligand) %>% unique()
  
  prioritized_tbl_oi = prioritization_tables$prioritization_tbl_ligand_receptor %>% filter(ligand %in% filtered_ligands) %>% select(niche, sender, receiver, ligand,  receptor, ligand_receptor, prioritization_score) %>% distinct() %>% inner_join(top_ligand_receptor_niche_df) %>% group_by(ligand) %>% filter(receiver == receiver_oi) %>% top_n(2, prioritization_score) %>% ungroup() 
  
  exprs_activity_target_plot = make_ligand_activity_target_exprs_plot(receiver_oi, prioritized_tbl_oi,  prioritization_tables$prioritization_tbl_ligand_receptor,  prioritization_tables$prioritization_tbl_ligand_target, output$exprs_tbl_ligand,  output$exprs_tbl_target, lfc_cutoff, ligand_target_matrix, plot_legend = FALSE, heights = NULL, widths = NULL)
  exprs_activity_target_plot$combined_plot
  
  p2 <- exprs_activity_target_plot$combined_plot
  p2
  
  #Circos plot of prioritized ligand-receptor pairs
  receiver_oi = "AM1_SFRP5+ epi"
  
  filtered_ligands = ligand_prioritized_tbl_oi %>% filter(receiver == receiver_oi) %>% top_n(30, prioritization_score) %>% pull(ligand) %>% unique()
  
  prioritized_tbl_oi = prioritization_tables$prioritization_tbl_ligand_receptor %>% filter(ligand %in% filtered_ligands) %>% select(niche, sender, receiver, ligand,  receptor, ligand_receptor, prioritization_score) %>% distinct() %>% inner_join(top_ligand_receptor_niche_df) %>% group_by(ligand) %>% filter(receiver == receiver_oi) %>% top_n(2, prioritization_score) %>% ungroup() 

  colors_sender=c(
    "AM1_SFRP5+ epi"='#C91105FF',
    "AM1_CNN1+ stro"= '#08519CFF',
    "AM1_Endothelial"='#EACC62FF',
    "AM1_Immune" = '#469D76FF',
    "AM1_DES+ SMC" = '#E8D8E8FF',
    'AM1_ESR1+ SMC'='#583070FF',
    "AM1_PV" = '#924099FF'
  )
  colors_receiver=c("AM1_SFRP5+ epi"='#C91105FF')
  
  write.csv(prioritized_tbl_oi,'figures/invade_env/IHH_circos_more.csv')
  
  prioritized_tbl_oi <- read_csv('figures/invade_env/IHH_circos_more_EMT.csv')
  
  prioritized_tbl_oi$sender <- factor(prioritized_tbl_oi$sender,levels=c("AM1_SFRP5+ epi",
                                                                         "AM1_CNN1+ stro",
                                                                         "AM1_Endothelial","AM1_Immune",
                                                                         "AM1_DES+ SMC","AM1_ESR1+ SMC",'AM1_PV'))
  
  p <-  make_circos_lr(prioritized_tbl_oi, colors_sender, colors_receiver)
  
  p1 <- p$p_legend
  p1


