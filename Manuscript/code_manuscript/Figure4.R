#!/usr/bin/env Rscript

#########
# NOTES #
#########
# This code was used to produce Figure 4 of the ensemblex manuscript

########
# Main #
########
## Load libraries
packages <- c('dplyr', 'tidyr' 'pdfCluster', 'data.table','readr','lubridate', 'tidyverse', 'moments', 'mousetrap', 'usethis', 'devtools', 'desc', 'kneedle', 'ROCit', 'ggplot2', 'factoextra', 'ggpubr', 'ComplexUpset', 'Seurat', 'Matrix', 'scCustomize')
lapply(packages, library, character.only = TRUE)

output_dir <- "~/ensemblex_manuscript/Figure4"

########################
# create Seurat object #
########################
## code
    data_dir <- "~/raw_feature_bc_matrix/"
    data <- Read10X(data.dir = data_dir)
    seurat_object = CreateSeuratObject(counts = data$'Gene Expression')
    seurat_object[['CMO']] = CreateAssayObject(counts = data$'Multiplexing Capture')
    saveRDS(seurat_object, '~/ensemblex_manuscript/Figure4/cell_ranger_multi.rds')
##

#############################
# Demultiplex with HTOdemux #
#############################
## code
    output_dir <- "~/ensemblex_manuscript/Figure4"
    seurat_object<-readRDS('~/ensemblex_manuscript/Figure4/cell_ranger_multi.rds')

    confidence <- fread("~/20k_NSCLC_DTC_3p_nextgem_intron_Multiplex_multiplexing_analysis_assignment_confidence_table.csv")
    cells <- fread("~/CellRanger_multi/20k_NSCLC_DTC_3p_nextgem_intron_Multiplex_multiplexing_analysis_assignment_confidence_table.csv",select = c("Barcodes"))
    seurat_object_use <- subset(seurat_object, cells = cells$Barcodes)

    features = c("CMO301", "CMO302", "CMO303", "CMO304","CMO306", "CMO307", "CMO308")
    DefaultAssay(seurat_object_use)<-"CMO"
    keep <- rownames(seurat_object_use)[rownames(seurat_object_use) %in% features]
    counts <- GetAssayData(seurat_object_use, assay = "CMO")
    counts <- counts[(which(rownames(counts) %in% keep)),]
    seurat_object_use[["CMO"]] <- subset(seurat_object_use[["CMO"]], features = rownames(counts)) 

    seurat_object_use = subset(x = seurat_object_use, subset = nCount_CMO > 0) 

    seurat_object_norm <- NormalizeData(seurat_object_use, assay = "CMO", normalization.method = "CLR")
    seurat_object_demux <- HTODemux(seurat_object_norm, assay = "CMO", positive.quantile = 0.99)

    ## CMO Ridgeplot
    seu_temp <- seurat_object_demux
    DefaultAssay(seu_temp)
    seu_temp <- ScaleData(seu_temp)
    RidgePlot(seu_temp, assay = "CMO", features = rownames(seu_temp[["CMO"]]), group.by = "hash.ID", ncol =4, cols = c("black","#CAB2D6","#A6CEE3","#FDBF6F", "#FB9A99", "#1F78B4", "#B2DF8A","#33A02C","lightgrey" ))
    ggsave(paste(output_dir,'/Ridge_HTO.pdf',sep=""), width = 10, height = 7)

    ## Seurat analysis
    seurat_object_demux_neg <- subset(seurat_object_demux, CMO_classification.global != "Negative")
    keep_barcodes <- unlist(seurat_object_demux_neg@assays$CMO@counts@Dimnames[2])
    keep_barcodes <- data.frame(keep_barcodes)
    DefaultAssay(seurat_object_demux) <- "CMO"
    seurat_object_demux <- subset(seurat_object_demux, cells = keep_barcodes$keep_barcodes)

    ## Calculate a tSNE embedding of the CMO data
    DefaultAssay(seurat_object_demux) <- "CMO"
    seurat_object_demux <- ScaleData(seurat_object_demux, features = rownames(seurat_object_demux),
        verbose = FALSE)
    seurat_object_demux <- RunPCA(seurat_object_demux, features = rownames(seurat_object_demux), approx = FALSE)
    seurat_object_demux <- RunTSNE(seurat_object_demux, dims = 1:7, perplexity = 100)

    ## Plot CMO tSNE
    DimPlot(seurat_object_demux) +
        theme_void() +
        theme(legend.position = "none") +
        scale_colour_manual(values = c("Doublet" = "black", "CMO303" = "#CAB2D6", "CMO307" = "#A6CEE3", "CMO302" =  "#FDBF6F", "CMO306" =  "#FB9A99", "CMO308" = "#1F78B4", "CMO301" =  "#B2DF8A", "CMO304" =  "#33A02C"))
    ggsave(paste(output_dir,"/TSNE_CMO.pdf", sep=""),width = 5, height = 5, dpi = 300)

    ## Add ensemblex metadata
    ensemblex <- read.delim("~/EID/PWE_GBD_EID_droplet_assignment.csv", sep=",")
    ensemblex <- ensemblex %>% dplyr::select(barcode, assignment_R_consensus)
    rownames(ensemblex) <- unlist(seurat_object_demux@assays$CMO@counts@Dimnames[2])
    seurat_object_demux <- AddMetaData(seurat_object_demux, ensemblex)

    ## ensemblex tSNE all cells
    DimPlot(seurat_object_demux, group.by = "assignment_R_consensus") +
        theme_void() +
        theme(legend.position = "none",
        plot.title = element_blank()) +
        scale_colour_manual(values = c("Doublet" = "black", "CMO303" = "#CAB2D6", "CMO307" = "#A6CEE3", "CMO302" =  "#FDBF6F", "CMO306" =  "#FB9A99", "CMO308" = "#1F78B4", "CMO301" =  "#B2DF8A", "CMO304" =  "#33A02C"))
    ggsave(paste(output_dir,"/TSNE_ensemblex_all_cells.pdf", sep=""),width = 5, height = 5, dpi = 300)

    ### Evaluate doublet
    seurat_object_demux_doublet <- seurat_object_demux

    ## Add ensemblex metadata
    ensemblex <- read.delim("~/EID/PWE_GBD_EID_droplet_assignment.csv", sep=",")
    ensemblex <- ensemblex %>% dplyr::select(barcode, assignment_R_consensus, truth)
    rownames(ensemblex) <- unlist(seurat_object_demux_doublet@assays$CMO@counts@Dimnames[2])

    ## Evaluate
    ensemblex$eval <- "singlet"
    ensemblex$eval[ensemblex$assignment_R_consensus == "Doublet" & ensemblex$truth == "Doublet"] <- "TP"
    ensemblex$eval[ensemblex$assignment_R_consensus != "Doublet" & ensemblex$truth == "Doublet"] <- "FN"
    ensemblex$eval[ensemblex$assignment_R_consensus == "Doublet" & ensemblex$truth != "Doublet"] <- "FP"

    seurat_object_demux_doublet <- AddMetaData(seurat_object_demux_doublet, ensemblex)
    seurat_object_demux_doublet <- subset(seurat_object_demux_doublet, eval != "singlet")

    ## Plot
    DimPlot(seurat_object_demux_doublet, group.by = "eval") +
        theme_void() +
        theme(legend.position = "none",
        plot.title = element_blank()) +
        scale_colour_manual(values = c("#fdda54","#f64941",  "#a6e0e9", "grey"))
    ggsave(paste(output_dir,"/TSNE_ensemblex_doublet.pdf", sep=""),width = 3, height = 3, dpi = 300)


    ### Evaluate singlet
    seurat_object_demux_singlet <- seurat_object_demux

    ## Add ensemblex metadata
    ensemblex <- read.delim("~/EID/PWE_GBD_EID_droplet_assignment.csv", sep=",")
    ensemblex <- ensemblex %>% dplyr::select(barcode, assignment_R_consensus, truth)
    rownames(ensemblex) <- unlist(seurat_object_demux_singlet@assays$CMO@counts@Dimnames[2])

    ## Evaluate
    ensemblex$eval <- "doublet"
    ensemblex$eval[ensemblex$assignment_R_consensus == ensemblex$truth & ensemblex$truth != "Doublet"] <- "TP"
    ensemblex$eval[ensemblex$assignment_R_consensus != ensemblex$truth & ensemblex$truth != "Doublet"] <- "FP"

    seurat_object_demux_singlet <- AddMetaData(seurat_object_demux_singlet, ensemblex)
    seurat_object_demux_singlet <- subset(seurat_object_demux_singlet, eval != "doublet")

    ## Plot
    DimPlot(seurat_object_demux_singlet, group.by = "eval") +
        theme_void() +
        theme(legend.position = "none",
        plot.title = element_blank()) +
        scale_colour_manual(values = c("#f64941", "#a6e0e9"))
    ggsave(paste(output_dir,"/TSNE_ensemblex_singlet.pdf", sep=""),width = 3, height = 3, dpi = 300)
##

################################
# Singlet classification rates #
################################
## code
    ## Load in ensemblex outputs
    ensemblex <- read.delim("~/EID/PWE_GBD_EID_droplet_assignment.csv", sep=",")

    ### Compute how many true singlets we have
    ensemblex_singlets <- subset(ensemblex, truth != "Doublet")
    n_singlets = nrow(ensemblex_singlets)

    ## ensemblex
    ensemblex$eval_ensemblex <- "Doublet"
    ensemblex$eval_ensemblex[ensemblex$assignment_R_consensus != "Doublet" & ensemblex$truth != "Doublet" & ensemblex$truth == ensemblex$assignment_R_consensus] <- "TP"
    ensemblex$eval_ensemblex[ensemblex$assignment_R_consensus != "Doublet" & ensemblex$truth != "Doublet" & ensemblex$truth != ensemblex$assignment_R_consensus] <- "FP"
    df <- data.frame(table(ensemblex$eval_ensemblex))
    df <- subset(df, Var1 != "Doublet")
    df$Freq <- df$Freq/n_singlets
    df$Trool <- "ensemblex"
    df_ensemblex <- df

    ## Demuxalot
    ensemblex$eval_ensemblex <- "Doublet"
    ensemblex$eval_ensemblex[ensemblex$demuxalot != "Doublet" & ensemblex$truth != "Doublet" & ensemblex$truth == ensemblex$demuxalot] <- "TP"
    ensemblex$eval_ensemblex[ensemblex$demuxalot != "Doublet" & ensemblex$truth != "Doublet" & ensemblex$truth != ensemblex$demuxalot] <- "FP"
    df <- data.frame(table(ensemblex$eval_ensemblex))
    df <- subset(df, Var1 != "Doublet")
    df$Freq <- df$Freq/n_singlets
    df$Trool <- "demuxalot"
    df_demuxalot <- df

    ## Demuxlet
    ensemblex$eval_ensemblex <- "Doublet"
    ensemblex$eval_ensemblex[ensemblex$freemuxlet_sample != "Doublet" & ensemblex$truth != "Doublet" & ensemblex$truth == ensemblex$freemuxlet_sample] <- "TP"
    ensemblex$eval_ensemblex[ensemblex$freemuxlet_sample != "Doublet" & ensemblex$truth != "Doublet" & ensemblex$truth != ensemblex$freemuxlet_sample] <- "FP"
    df <- data.frame(table(ensemblex$eval_ensemblex))
    df <- subset(df, Var1 != "Doublet")
    df$Freq <- df$Freq/n_singlets
    df$Trool <- "freemuxlet"
    df_freemuxlet <- df

    ## Souporcell
    ensemblex$eval_ensemblex <- "Doublet"
    ensemblex$eval_ensemblex[ensemblex$souporcell_sample != "Doublet" & ensemblex$truth != "Doublet" & ensemblex$truth == ensemblex$souporcell_sample] <- "TP"
    ensemblex$eval_ensemblex[ensemblex$souporcell_sample != "Doublet" & ensemblex$truth != "Doublet" & ensemblex$truth != ensemblex$souporcell_sample] <- "FP"
    df <- data.frame(table(ensemblex$eval_ensemblex))
    df <- subset(df, Var1 != "Doublet")
    df$Freq <- df$Freq/n_singlets
    df$Trool <- "souporcell"
    df_souporcell <- df

    ## Vireo
    ensemblex$eval_ensemblex <- "Doublet"
    ensemblex$eval_ensemblex[ensemblex$vireo_sample != "Doublet" & ensemblex$truth != "Doublet" & ensemblex$truth == ensemblex$vireo_sample] <- "TP"
    ensemblex$eval_ensemblex[ensemblex$vireo_sample != "Doublet" & ensemblex$truth != "Doublet" & ensemblex$truth != ensemblex$vireo_sample] <- "FP"
    df <- data.frame(table(ensemblex$eval_ensemblex))
    df <- subset(df, Var1 != "Doublet")
    df$Freq <- df$Freq/n_singlets
    df$Trool <- "vireo"
    df_vireo <- df

    ## Merge
    df_merge <- rbind(df_ensemblex, df_demuxalot, df_freemuxlet, df_souporcell, df_vireo)
    df_merge$Trool <- factor(df_merge$Trool, levels = c("ensemblex",  "demuxalot",  "freemuxlet", "souporcell", "vireo"))
    df_merge$Var1 <- factor(df_merge$Var1, levels = c("TP", "FP"))

    ## Plot
    doublet_TP <- ggplot(df_merge, aes(x = Freq, y = Trool, fill = Trool, label = round(Freq, digits = 3))) +
        theme_bw() +
        theme(legend.position = "none",
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 12, face = "bold", colour = "black"),
        panel.grid = element_blank(),
        strip.background=element_rect(colour="black",
                                    fill="white")) +
        scale_x_continuous(expand = c(0,0.01), limits = c(0,1.15)) +
        scale_y_discrete(labels = c("ensemblex",  "Demuxalot",  "Demuxlet", "Souporcell", "Vireo")) +
        geom_bar(stat = "identity") +
        geom_text(hjust=-0.25) +
        xlab("Classification rate") + 
        facet_wrap(.~Var1, scales = "free_y", ncol = 1)+
        scale_fill_manual(values = c("black", "#d95f02", "#e6ab02", "#7570b3", "#66a61e"))
    ggsave(paste(output_dir,"/Singlet_rate.pdf", sep=""),width =6, height = 2.75, dpi = 300)
##

################################
# Doublet classification rates #
################################
## code
    ## Load in ensemblex outputs
    ensemblex <- read.delim("~/EID/PWE_GBD_EID_droplet_assignment.csv", sep=",")

    ### Compute how many true doublets we have
    ensemblex_doublet <- subset(ensemblex, truth == "Doublet")

    n_doublets = nrow(ensemblex_doublet)
    n_droplets = nrow(ensemblex)

    ## ensemblex
    ensemblex$eval_ensemblex <- "singlet"
    ensemblex$eval_ensemblex[ensemblex$assignment_R_consensus == "Doublet" & ensemblex$truth == "Doublet"] <- "TP"
    ensemblex$eval_ensemblex[ensemblex$assignment_R_consensus != "Doublet" & ensemblex$truth == "Doublet"] <- "FN"
    ensemblex$eval_ensemblex[ensemblex$assignment_R_consensus == "Doublet" & ensemblex$truth != "Doublet"] <- "FP"
    df <- data.frame(table(ensemblex$eval_ensemblex))
    df <- subset(df, Var1 != "singlet")
    df$Freq[df$Var1 == "TP"] <- df$Freq[df$Var1 == "TP"]/n_doublets
    df$Freq[df$Var1 == "FP"] <- df$Freq[df$Var1 == "FP"]/n_droplets
    df$Trool <- "ensemblex"
    df_ensemblex <- df

    ## Demuxalot
    ensemblex$eval_ensemblex <- "singlet"
    ensemblex$eval_ensemblex[ensemblex$demuxalot == "Doublet" & ensemblex$truth == "Doublet"] <- "TP"
    ensemblex$eval_ensemblex[ensemblex$demuxalot != "Doublet" & ensemblex$truth == "Doublet"] <- "FN"
    ensemblex$eval_ensemblex[ensemblex$demuxalot == "Doublet" & ensemblex$truth != "Doublet"] <- "FP"
    df <- data.frame(table(ensemblex$eval_ensemblex))
    df <- subset(df, Var1 != "singlet")
    df$Freq[df$Var1 == "TP"] <- df$Freq[df$Var1 == "TP"]/n_doublets
    df$Freq[df$Var1 == "FP"] <- df$Freq[df$Var1 == "FP"]/n_droplets
    df$Trool <- "demuxalot"
    df_demuxalot <- df

    ## Demuxlet
    ensemblex$eval_ensemblex <- "singlet"
    ensemblex$eval_ensemblex[ensemblex$freemuxlet == "Doublet" & ensemblex$truth == "Doublet"] <- "TP"
    ensemblex$eval_ensemblex[ensemblex$freemuxlet != "Doublet" & ensemblex$truth == "Doublet"] <- "FN"
    ensemblex$eval_ensemblex[ensemblex$freemuxlet == "Doublet" & ensemblex$truth != "Doublet"] <- "FP"
    df <- data.frame(table(ensemblex$eval_ensemblex))
    df <- subset(df, Var1 != "singlet")
    df$Freq[df$Var1 == "TP"] <- df$Freq[df$Var1 == "TP"]/n_doublets
    df$Freq[df$Var1 == "FP"] <- df$Freq[df$Var1 == "FP"]/n_droplets
    df$Trool <- "freemuxlet"
    df_freemuxlet <- df

    ## Souporcell
    ensemblex$eval_ensemblex <- "singlet"
    ensemblex$eval_ensemblex[ensemblex$souporcell == "Doublet" & ensemblex$truth == "Doublet"] <- "TP"
    ensemblex$eval_ensemblex[ensemblex$souporcell != "Doublet" & ensemblex$truth == "Doublet"] <- "FN"
    ensemblex$eval_ensemblex[ensemblex$souporcell == "Doublet" & ensemblex$truth != "Doublet"] <- "FP"
    df <- data.frame(table(ensemblex$eval_ensemblex))
    df <- subset(df, Var1 != "singlet")
    df$Freq[df$Var1 == "TP"] <- df$Freq[df$Var1 == "TP"]/n_doublets
    df$Freq[df$Var1 == "FP"] <- df$Freq[df$Var1 == "FP"]/n_droplets
    df$Trool <- "souporcell"
    df_souporcell <- df

    ## Vireo
    ensemblex$eval_ensemblex <- "singlet"
    ensemblex$eval_ensemblex[ensemblex$vireo == "Doublet" & ensemblex$truth == "Doublet"] <- "TP"
    ensemblex$eval_ensemblex[ensemblex$vireo != "Doublet" & ensemblex$truth == "Doublet"] <- "FN"
    ensemblex$eval_ensemblex[ensemblex$vireo == "Doublet" & ensemblex$truth != "Doublet"] <- "FP"
    df <- data.frame(table(ensemblex$eval_ensemblex))
    df <- subset(df, Var1 != "singlet")
    df$Freq[df$Var1 == "TP"] <- df$Freq[df$Var1 == "TP"]/n_doublets
    df$Freq[df$Var1 == "FP"] <- df$Freq[df$Var1 == "FP"]/n_droplets
    df$Trool <- "vireo"
    df_vireo <- df

    ## Merge
    df_merge <- rbind(df_ensemblex, df_demuxalot, df_freemuxlet, df_souporcell, df_vireo)
    df_merge$Trool <- factor(df_merge$Trool, levels = c("ensemblex",  "demuxalot",  "freemuxlet", "souporcell", "vireo"))
    df_merge <- subset(df_merge, Var1 != "FN")
    df_merge$Var1 <- factor(df_merge$Var1, levels = c("TP", "FP"))

    ## Plot
    doublet_TP <- ggplot(df_merge, aes(x = Freq, y = Trool, fill = Trool, label = round(Freq, digits = 3))) +
        theme_bw() +
        theme(legend.position = "none",
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 12, face = "bold", colour = "black"),
        panel.grid = element_blank(),
        strip.background=element_rect(colour="black",
                                        fill="white")) +
        scale_x_continuous(expand = c(0,0.01), limits = c(0,0.8)) +
        scale_y_discrete(labels = c("ensemblex",  "Demuxalot",  "Demuxlet", "Souporcell", "Vireo")) +
        geom_bar(stat = "identity") +
        geom_text(hjust=-0.25) +
        xlab("Classification rate") + 
        facet_wrap(.~Var1, scales = "free_y", ncol = 1)+
        scale_fill_manual(values = c("black", "#d95f02", "#e6ab02", "#7570b3", "#66a61e"))
    ggsave(paste(output_dir,"/Doublet_rate.pdf", sep=""),width =6, height = 2.75, dpi = 300)
##

##############################
# Proportion of usable cells #
##############################
## code
    eval_df <- read.delim("~/EID/PWE_GBD_EID_droplet_assignment.csv", sep=",")
    test_prop <- eval_df
            
    ## Demuxalot
    test_prop$demuxalot_test <- test_prop$demuxalot
    test_prop$demuxalot_test[test_prop$demuxalot_max < 0.9] <- "unassigned"
    prop_demuxalot  <- subset(test_prop, test_prop$demuxalot_test != "unassigned" & test_prop$demuxalot_test !="Doublet") 
    prop_usable_demuxalot <- nrow(prop_demuxalot)/nrow(test_prop)
    prop_usable_demuxalot

    ## Demuxlet
    prop_demuxlet  <- subset(test_prop, test_prop$freemuxlet_sample != "unassigned" & test_prop$freemuxlet_sample !="Doublet") 
    prop_usable_demuxlet <- nrow(prop_demuxlet)/nrow(test_prop)
    prop_usable_demuxlet

    ## Vireo
    prop_vireo  <- subset(test_prop, test_prop$vireo_sample != "unassigned" & test_prop$vireo_sample !="Doublet") 
    prop_usable_vireo <- nrow(prop_vireo)/nrow(test_prop)
    prop_usable_vireo

    ## Souporcell
    prop_souporcell  <- subset(test_prop, test_prop$souporcell_sample != "unassigned" & test_prop$souporcell_sample !="Doublet") 
    prop_usable_souporcell <- nrow(prop_souporcell)/nrow(test_prop)
    prop_usable_souporcell

    ## ensemblex
    test_prop$assignment_R_consensus[test_prop$singlet_confidence < 1.0] <- "unassigned"
    prop_ensemblex  <- subset(test_prop, test_prop$assignment_R_consensus != "unassigned" & test_prop$assignment_R_consensus !="Doublet") 
    prop_usable_ensemblex <- nrow(prop_ensemblex)/nrow(test_prop)
    prop_usable_ensemblex

    ## Dataframe
    tool <- c("ensemblex",  "Demuxalot", "Demuxlet", "Souporcell", "Vireo")
    prop_usable <- c(prop_usable_ensemblex, prop_usable_demuxalot, prop_usable_demuxlet, prop_usable_souporcell, prop_usable_vireo)
    prop_usable_df <- data.frame(tool, prop_usable)
    prop_usable_df$tool <- factor(prop_usable_df$tool, levels = c("ensemblex",  "Demuxalot", "Demuxlet", "Souporcell", "Vireo"))
##

###########################
# Usable cells error rate #
###########################
## code
    eval_df <- read.delim("~/EID/PWE_GBD_EID_droplet_assignment.csv", sep=",")
    test_prop <- eval_df
                
    ## Demuxalot
    test_prop$demuxalot_test <- test_prop$demuxalot
    test_prop$demuxalot_test[test_prop$demuxalot_max < 0.9] <- "unassigned"
    prop_demuxalot  <- subset(test_prop, test_prop$demuxalot_test != "unassigned" & test_prop$demuxalot_test !="Doublet") 
    prop_demuxalot_correct <- subset(prop_demuxalot, demuxalot_test == truth)
    prop_usable_error_demuxalot <- 1- nrow(prop_demuxalot_correct)/nrow(prop_demuxalot)
    prop_usable_error_demuxalot

    ## Demuxlet
    prop_demuxlet  <- subset(test_prop, test_prop$freemuxlet_sample != "unassigned" & test_prop$freemuxlet_sample !="Doublet") 
    prop_demuxlet_correct <- subset(prop_demuxlet, freemuxlet_sample == truth)
    prop_usable_error_demuxlet <- 1-nrow(prop_demuxlet_correct)/nrow(prop_demuxlet)
    prop_usable_error_demuxlet

    ## Vireo
    prop_vireo  <- subset(test_prop, test_prop$vireo_sample != "unassigned" & test_prop$vireo_sample !="Doublet") 
    prop_vireo_correct <- subset(prop_vireo, vireo_sample == truth)
    prop_usable_error_vireo <- 1-nrow(prop_vireo_correct)/nrow(prop_vireo)
    prop_usable_error_vireo

    ## Souporcell
    prop_souporcell  <- subset(test_prop, test_prop$souporcell_sample != "unassigned" & test_prop$souporcell_sample !="Doublet") 
    prop_souporcell_correct <- subset(prop_souporcell, souporcell_sample == truth)
    prop_usable_error_souporcell <- 1- nrow(prop_souporcell_correct)/nrow(prop_souporcell)
    prop_usable_error_souporcell

    ## ensemblex
    test_prop$assignment_R_consensus[test_prop$singlet_confidence < 1.0] <- "unassigned"
    prop_ensemblex  <- subset(test_prop, test_prop$assignment_R_consensus != "unassigned" & test_prop$assignment_R_consensus !="Doublet") 
    prop_ensemblex_correct <- subset(prop_ensemblex, assignment_R_consensus == truth)
    prop_usable_error_ensemblex <- 1- nrow(prop_ensemblex_correct)/nrow(prop_ensemblex)
    prop_usable_error_ensemblex
        
    ## Plot
    tool <- c("ensemblex",  "Demuxalot", "Demuxlet", "Souporcell", "Vireo")
    prop_usable_error <- c(prop_usable_error_ensemblex, prop_usable_error_demuxalot, prop_usable_error_demuxlet, prop_usable_error_souporcell, prop_usable_error_vireo)
    prop_usable_error_df <- data.frame(tool, prop_usable_error)
    prop_usable_error_df$tool <- factor(prop_usable_error_df$tool, levels = c("ensemblex",  "Demuxalot", "Demuxlet", "Souporcell", "Vireo"))
##

#######################
# Adjusted Rand Index #
#######################
## code
    eval_df <- read.delim("~/EID/PWE_GBD_EID_droplet_assignment.csv", sep=",")
    test_prop <- eval_df

    ### Apply threshold for ensemblex and Demuxalot and make sure best guess is doublet even if "unassigned"
    ## Demuxalot
    test_prop$demuxalot_test <- test_prop$demuxalot
    test_prop$demuxalot_test[test_prop$demuxalot_max < 0.9 & test_prop$demuxalot != "Doublet"] <- "unassigned"

    ## ensemblex
    test_prop$assignment_R_consensus[test_prop$singlet_confidence < 1.0 & test_prop$assignment_R_consensus != "Doublet"] <- "unassigned"

    ## Souporcell
    test_prop$souporcell_sample[test_prop$souporcell == "Doublet"] <- "Doublet"

    ## Demuxlet
    test_prop$freemuxlet_sample[test_prop$freemuxlet == "Doublet"] <- "Doublet"

    ## Vireo
    test_prop$vireo_sample[test_prop$vireo == "Doublet"] <- "Doublet"

    ## calculate ARI
    ARI_vireo_truth <- adj.rand.index(test_prop$vireo_sample, test_prop$truth) 
    ARI_ensemblex_truth <- adj.rand.index(test_prop$assignment_R_consensus, test_prop$truth) 
    ARI_demuxalot_truth <- adj.rand.index(test_prop$demuxalot_test, test_prop$truth) 
    ARI_demuxlet_truth <- adj.rand.index(test_prop$freemuxlet_sample, test_prop$truth) 
    ARI_souporcell_truth <- adj.rand.index(test_prop$souporcell_sample, test_prop$truth) 

    ## Dataframe
    tool <- c("ensemblex",  "Demuxalot", "Demuxlet", "Souporcell", "Vireo")
    ARI_thresh <- c(ARI_ensemblex_truth, ARI_demuxalot_truth, ARI_demuxlet_truth, ARI_souporcell_truth, ARI_vireo_truth)
    ARI_threshold_df <- data.frame(tool, ARI_thresh)
    ARI_threshold_df$tool <- factor(ARI_threshold_df$tool, levels = c("ensemblex", "Demuxalot", "Demuxlet", "Souporcell", "Vireo"))
##

######################
# Balanced accuracy #
######################
## code
    eval_df <- read.delim("~/EID/PWE_GBD_EID_droplet_assignment.csv", sep=",")
    eval_df_ba <- eval_df

    ### Apply threshold for ensemblex and Demuxalot and make sure best guess is doublet even if "unassigned"
    ## Demuxalot
    eval_df_ba$demuxalot_test <- eval_df_ba$demuxalot
    eval_df_ba$demuxalot_test[eval_df_ba$demuxalot_max < 0.9 & eval_df_ba$demuxalot != "Doublet"] <- "unassigned"

    ## ensemblex
    eval_df_ba$assignment_R_consensus[eval_df_ba$singlet_confidence < 1.0 & eval_df_ba$assignment_R_consensus != "Doublet"] <- "unassigned"

    ## Souporcell
    eval_df_ba$souporcell_sample[eval_df_ba$souporcell == "Doublet"] <- "Doublet"

    ## Demuxlet
    eval_df_ba$freemuxlet_sample[eval_df_ba$freemuxlet == "Doublet"] <- "Doublet"

    ## Vireo
    eval_df_ba$vireo_sample[eval_df_ba$vireo == "Doublet"] <- "Doublet"

    ### Compute balanced accuracy
    ## Vireo
    eval_df_ba$vireo_eval[str_detect(eval_df_ba$truth, "D") & eval_df_ba$truth == eval_df_ba$vireo_sample] <- "TN"
    eval_df_ba$vireo_eval[str_detect(eval_df_ba$truth, "C") & eval_df_ba$truth == eval_df_ba$vireo_sample] <- "TP"
    eval_df_ba$vireo_eval[str_detect(eval_df_ba$truth, "C") & eval_df_ba$truth != eval_df_ba$vireo_sample] <- "FP"
    eval_df_ba$vireo_eval[str_detect(eval_df_ba$truth, "D") & eval_df_ba$truth != eval_df_ba$vireo_sample] <- "FN"
    df <- eval_df_ba
    df_summary <- data.frame(table(df$vireo_eval))
    Var1 <- c("FN", "TN", "TP", "FP")
    Freq<- c(0,0,0,0)
    filler_frame <- data.frame(Var1, Freq)
    df_summary <- rbind(df_summary, filler_frame)
    df_summary <- df_summary %>% dplyr::group_by(Var1) %>% 
        dplyr::summarise(Freq = sum(Freq)) %>% as.data.frame()
    BA <- ((df_summary$Freq[df_summary$Var1 == "TP"]/(df_summary$Freq[df_summary$Var1 == "TP"] + df_summary$Freq[df_summary$Var1 == "FN"])) + 
            (df_summary$Freq[df_summary$Var1 == "TN"]/(df_summary$Freq[df_summary$Var1 == "TN"] + df_summary$Freq[df_summary$Var1 == "FP"])))/2
    vireo_BA <- BA

    ## Demuxalot
    eval_df_ba$demuxalot_eval[str_detect(eval_df_ba$truth, "D") & eval_df_ba$truth == eval_df_ba$demuxalot_test] <- "TN"
    eval_df_ba$demuxalot_eval[str_detect(eval_df_ba$truth, "C") & eval_df_ba$truth == eval_df_ba$demuxalot_test] <- "TP"
    eval_df_ba$demuxalot_eval[str_detect(eval_df_ba$truth, "C") & eval_df_ba$truth != eval_df_ba$demuxalot_test] <- "FP"
    eval_df_ba$demuxalot_eval[str_detect(eval_df_ba$truth, "D") & eval_df_ba$truth != eval_df_ba$demuxalot_test] <- "FN"
    df <- eval_df_ba
    df_summary <- data.frame(table(df$demuxalot_eval))
    Var1 <- c("FN", "TN", "TP", "FP")
    Freq<- c(0,0,0,0)
    filler_frame <- data.frame(Var1, Freq)
    df_summary <- rbind(df_summary, filler_frame)
    df_summary <- df_summary %>% dplyr::group_by(Var1) %>% 
        dplyr::summarise(Freq = sum(Freq)) %>% as.data.frame()
    BA <- ((df_summary$Freq[df_summary$Var1 == "TP"]/(df_summary$Freq[df_summary$Var1 == "TP"] + df_summary$Freq[df_summary$Var1 == "FN"])) + 
            (df_summary$Freq[df_summary$Var1 == "TN"]/(df_summary$Freq[df_summary$Var1 == "TN"] + df_summary$Freq[df_summary$Var1 == "FP"])))/2
    demuxalot_BA <- BA

    ## ensemblex
    eval_df_ba$ensemblex_eval[str_detect(eval_df_ba$truth, "D") & eval_df_ba$truth == eval_df_ba$assignment_R_consensus] <- "TN"
    eval_df_ba$ensemblex_eval[str_detect(eval_df_ba$truth, "C") & eval_df_ba$truth == eval_df_ba$assignment_R_consensus] <- "TP"
    eval_df_ba$ensemblex_eval[str_detect(eval_df_ba$truth, "C") & eval_df_ba$truth != eval_df_ba$assignment_R_consensus] <- "FP"
    eval_df_ba$ensemblex_eval[str_detect(eval_df_ba$truth, "D") & eval_df_ba$truth != eval_df_ba$assignment_R_consensus] <- "FN"
    df <- eval_df_ba
    df_summary <- data.frame(table(df$ensemblex_eval))
    Var1 <- c("FN", "TN", "TP", "FP")
    Freq<- c(0,0,0,0)
    filler_frame <- data.frame(Var1, Freq)
    df_summary <- rbind(df_summary, filler_frame)
    df_summary <- df_summary %>% dplyr::group_by(Var1) %>% 
        dplyr::summarise(Freq = sum(Freq)) %>% as.data.frame()
    BA <- ((df_summary$Freq[df_summary$Var1 == "TP"]/(df_summary$Freq[df_summary$Var1 == "TP"] + df_summary$Freq[df_summary$Var1 == "FN"])) + 
            (df_summary$Freq[df_summary$Var1 == "TN"]/(df_summary$Freq[df_summary$Var1 == "TN"] + df_summary$Freq[df_summary$Var1 == "FP"])))/2
    ensemblex_BA <- BA

    ## Demuxlet
    eval_df_ba$demuxlet_eval[str_detect(eval_df_ba$truth, "D") & eval_df_ba$truth == eval_df_ba$freemuxlet_sample] <- "TN"
    eval_df_ba$demuxlet_eval[str_detect(eval_df_ba$truth, "C") & eval_df_ba$truth == eval_df_ba$freemuxlet_sample] <- "TP"
    eval_df_ba$demuxlet_eval[str_detect(eval_df_ba$truth, "C") & eval_df_ba$truth != eval_df_ba$freemuxlet_sample] <- "FP"
    eval_df_ba$demuxlet_eval[str_detect(eval_df_ba$truth, "D") & eval_df_ba$truth != eval_df_ba$freemuxlet_sample] <- "FN"
    df <- eval_df_ba
    df_summary <- data.frame(table(df$demuxlet_eval))
    Var1 <- c("FN", "TN", "TP", "FP")
    Freq<- c(0,0,0,0)
    filler_frame <- data.frame(Var1, Freq)
    df_summary <- rbind(df_summary, filler_frame)
    df_summary <- df_summary %>% dplyr::group_by(Var1) %>% 
        dplyr::summarise(Freq = sum(Freq)) %>% as.data.frame()
    BA <- ((df_summary$Freq[df_summary$Var1 == "TP"]/(df_summary$Freq[df_summary$Var1 == "TP"] + df_summary$Freq[df_summary$Var1 == "FN"])) + 
            (df_summary$Freq[df_summary$Var1 == "TN"]/(df_summary$Freq[df_summary$Var1 == "TN"] + df_summary$Freq[df_summary$Var1 == "FP"])))/2
    demuxlet_BA <- BA

    ## Souporcell
    eval_df_ba$souporcell_eval[str_detect(eval_df_ba$truth, "D") & eval_df_ba$truth == eval_df_ba$souporcell_sample] <- "TN"
    eval_df_ba$souporcell_eval[str_detect(eval_df_ba$truth, "C") & eval_df_ba$truth == eval_df_ba$souporcell_sample] <- "TP"
    eval_df_ba$souporcell_eval[str_detect(eval_df_ba$truth, "C") & eval_df_ba$truth != eval_df_ba$souporcell_sample] <- "FP"
    eval_df_ba$souporcell_eval[str_detect(eval_df_ba$truth, "D") & eval_df_ba$truth != eval_df_ba$souporcell_sample] <- "FN"
    df <- eval_df_ba
    df_summary <- data.frame(table(df$souporcell_eval))
    Var1 <- c("FN", "TN", "TP", "FP")
    Freq<- c(0,0,0,0)
    filler_frame <- data.frame(Var1, Freq)
    df_summary <- rbind(df_summary, filler_frame)
    df_summary <- df_summary %>% dplyr::group_by(Var1) %>% 
        dplyr::summarise(Freq = sum(Freq)) %>% as.data.frame()
    BA <- ((df_summary$Freq[df_summary$Var1 == "TP"]/(df_summary$Freq[df_summary$Var1 == "TP"] + df_summary$Freq[df_summary$Var1 == "FN"])) + 
            (df_summary$Freq[df_summary$Var1 == "TN"]/(df_summary$Freq[df_summary$Var1 == "TN"] + df_summary$Freq[df_summary$Var1 == "FP"])))/2
    souporcell_BA <- BA

    ## Dataframe
    tool <- c("ensemblex",  "Demuxalot", "Demuxlet", "Souporcell", "Vireo")
    BA_thresh <- c(ensemblex_BA, demuxalot_BA, demuxlet_BA,souporcell_BA, vireo_BA)
    BA_threshold_df <- data.frame(tool, BA_thresh)
    BA_threshold_df$tool <- factor(BA_threshold_df$tool, levels = c("ensemblex",  "Demuxalot", "Demuxlet", "Souporcell", "Vireo"))
##

#####################################
# Matthew's Correlation Coefficient #
#####################################
## code
    eval_df <- read.delim("~/EID/PWE_GBD_EID_droplet_assignment.csv", sep=",")
    eval_df_ba <- eval_df

    ### Apply threshold for ensemblex and Demuxalot and make sure best guess is doublet even if "unassigned"
    ## Demuxalot
    eval_df_ba$demuxalot_test <- eval_df_ba$demuxalot
    eval_df_ba$demuxalot_test[eval_df_ba$demuxalot_max < 0.9 & eval_df_ba$demuxalot != "Doublet"] <- "unassigned"

    ## ensemblex
    eval_df_ba$assignment_R_consensus[eval_df_ba$singlet_confidence < 1.0 & eval_df_ba$assignment_R_consensus != "Doublet"] <- "unassigned"

    ## Souporcell
    eval_df_ba$souporcell_sample[eval_df_ba$souporcell == "Doublet"] <- "Doublet"

    ## Demuxlet
    eval_df_ba$freemuxlet_sample[eval_df_ba$freemuxlet == "Doublet"] <- "Doublet"

    ## Vireo
    eval_df_ba$vireo_sample[eval_df_ba$vireo == "Doublet"] <- "Doublet"

    ### Compute Matthew's Correlation Coefficient
    ## Vireo
    eval_df_ba$vireo_eval[str_detect(eval_df_ba$truth, "D") & eval_df_ba$truth == eval_df_ba$vireo_sample] <- "TN"
    eval_df_ba$vireo_eval[str_detect(eval_df_ba$truth, "C") & eval_df_ba$truth == eval_df_ba$vireo_sample] <- "TP"
    eval_df_ba$vireo_eval[str_detect(eval_df_ba$truth, "C") & eval_df_ba$truth != eval_df_ba$vireo_sample] <- "FP"
    eval_df_ba$vireo_eval[str_detect(eval_df_ba$truth, "D") & eval_df_ba$truth != eval_df_ba$vireo_sample] <- "FN"
    df <- eval_df_ba
    df_summary <- data.frame(table(df$vireo_eval))
    Var1 <- c("FN", "TN", "TP", "FP")
    Freq<- c(0,0,0,0)
    filler_frame <- data.frame(Var1, Freq)
    df_summary <- rbind(df_summary, filler_frame)
    df_summary <- df_summary %>% dplyr::group_by(Var1) %>% 
        dplyr::summarise(Freq = sum(Freq)) %>% as.data.frame()
    MCC <- ((as.numeric((df_summary$Freq[df_summary$Var1 == "TP"])*as.numeric(df_summary$Freq[df_summary$Var1 == "TN"]))-(as.numeric(df_summary$Freq[df_summary$Var1 == "FP"])*as.numeric(df_summary$Freq[df_summary$Var1 == "FN"])))/
                sqrt((as.numeric(df_summary$Freq[df_summary$Var1 == "TP"])+as.numeric(df_summary$Freq[df_summary$Var1 == "FP"]))*(as.numeric(df_summary$Freq[df_summary$Var1 == "TP"])+as.numeric(df_summary$Freq[df_summary$Var1 == "FN"]))*(as.numeric(df_summary$Freq[df_summary$Var1 == "TN"])+as.numeric(df_summary$Freq[df_summary$Var1 == "FP"]))*(as.numeric(df_summary$Freq[df_summary$Var1 == "TN"])+as.numeric(df_summary$Freq[df_summary$Var1 == "FN"]))))
    vireo_MCC <- MCC

    ## Demuxalot
    eval_df_ba$demuxalot_eval[str_detect(eval_df_ba$truth, "D") & eval_df_ba$truth == eval_df_ba$demuxalot_test] <- "TN"
    eval_df_ba$demuxalot_eval[str_detect(eval_df_ba$truth, "C") & eval_df_ba$truth == eval_df_ba$demuxalot_test] <- "TP"
    eval_df_ba$demuxalot_eval[str_detect(eval_df_ba$truth, "C") & eval_df_ba$truth != eval_df_ba$demuxalot_test] <- "FP"
    eval_df_ba$demuxalot_eval[str_detect(eval_df_ba$truth, "D") & eval_df_ba$truth != eval_df_ba$demuxalot_test] <- "FN"
    df <- eval_df_ba
    df_summary <- data.frame(table(df$demuxalot_eval))
    Var1 <- c("FN", "TN", "TP", "FP")
    Freq<- c(0,0,0,0)
    filler_frame <- data.frame(Var1, Freq)
    df_summary <- rbind(df_summary, filler_frame)
    df_summary <- df_summary %>% dplyr::group_by(Var1) %>% 
        dplyr::summarise(Freq = sum(Freq)) %>% as.data.frame()
    MCC <- ((as.numeric((df_summary$Freq[df_summary$Var1 == "TP"])*as.numeric(df_summary$Freq[df_summary$Var1 == "TN"]))-(as.numeric(df_summary$Freq[df_summary$Var1 == "FP"])*as.numeric(df_summary$Freq[df_summary$Var1 == "FN"])))/
                sqrt((as.numeric(df_summary$Freq[df_summary$Var1 == "TP"])+as.numeric(df_summary$Freq[df_summary$Var1 == "FP"]))*(as.numeric(df_summary$Freq[df_summary$Var1 == "TP"])+as.numeric(df_summary$Freq[df_summary$Var1 == "FN"]))*(as.numeric(df_summary$Freq[df_summary$Var1 == "TN"])+as.numeric(df_summary$Freq[df_summary$Var1 == "FP"]))*(as.numeric(df_summary$Freq[df_summary$Var1 == "TN"])+as.numeric(df_summary$Freq[df_summary$Var1 == "FN"]))))
    demuxalot_MCC <- MCC

    ## ensemblex
    eval_df_ba$ensemblex_eval[str_detect(eval_df_ba$truth, "D") & eval_df_ba$truth == eval_df_ba$assignment_R_consensus] <- "TN"
    eval_df_ba$ensemblex_eval[str_detect(eval_df_ba$truth, "C") & eval_df_ba$truth == eval_df_ba$assignment_R_consensus] <- "TP"
    eval_df_ba$ensemblex_eval[str_detect(eval_df_ba$truth, "C") & eval_df_ba$truth != eval_df_ba$assignment_R_consensus] <- "FP"
    eval_df_ba$ensemblex_eval[str_detect(eval_df_ba$truth, "D") & eval_df_ba$truth != eval_df_ba$assignment_R_consensus] <- "FN"
    df <- eval_df_ba
    df_summary <- data.frame(table(df$ensemblex_eval))
    Var1 <- c("FN", "TN", "TP", "FP")
    Freq<- c(0,0,0,0)
    filler_frame <- data.frame(Var1, Freq)
    df_summary <- rbind(df_summary, filler_frame)
    df_summary <- df_summary %>% dplyr::group_by(Var1) %>% 
        dplyr::summarise(Freq = sum(Freq)) %>% as.data.frame()
    MCC <- ((as.numeric((df_summary$Freq[df_summary$Var1 == "TP"])*as.numeric(df_summary$Freq[df_summary$Var1 == "TN"]))-(as.numeric(df_summary$Freq[df_summary$Var1 == "FP"])*as.numeric(df_summary$Freq[df_summary$Var1 == "FN"])))/
                sqrt((as.numeric(df_summary$Freq[df_summary$Var1 == "TP"])+as.numeric(df_summary$Freq[df_summary$Var1 == "FP"]))*(as.numeric(df_summary$Freq[df_summary$Var1 == "TP"])+as.numeric(df_summary$Freq[df_summary$Var1 == "FN"]))*(as.numeric(df_summary$Freq[df_summary$Var1 == "TN"])+as.numeric(df_summary$Freq[df_summary$Var1 == "FP"]))*(as.numeric(df_summary$Freq[df_summary$Var1 == "TN"])+as.numeric(df_summary$Freq[df_summary$Var1 == "FN"]))))
    ensemblex_MCC <- MCC

    ## Demuxlet
    eval_df_ba$demuxlet_eval[str_detect(eval_df_ba$truth, "D") & eval_df_ba$truth == eval_df_ba$freemuxlet_sample] <- "TN"
    eval_df_ba$demuxlet_eval[str_detect(eval_df_ba$truth, "C") & eval_df_ba$truth == eval_df_ba$freemuxlet_sample] <- "TP"
    eval_df_ba$demuxlet_eval[str_detect(eval_df_ba$truth, "C") & eval_df_ba$truth != eval_df_ba$freemuxlet_sample] <- "FP"
    eval_df_ba$demuxlet_eval[str_detect(eval_df_ba$truth, "D") & eval_df_ba$truth != eval_df_ba$freemuxlet_sample] <- "FN"
    df <- eval_df_ba
    df_summary <- data.frame(table(df$demuxlet_eval))
    Var1 <- c("FN", "TN", "TP", "FP")
    Freq<- c(0,0,0,0)
    filler_frame <- data.frame(Var1, Freq)
    df_summary <- rbind(df_summary, filler_frame)
    df_summary <- df_summary %>% dplyr::group_by(Var1) %>% 
        dplyr::summarise(Freq = sum(Freq)) %>% as.data.frame()
    MCC <- ((as.numeric((df_summary$Freq[df_summary$Var1 == "TP"])*as.numeric(df_summary$Freq[df_summary$Var1 == "TN"]))-(as.numeric(df_summary$Freq[df_summary$Var1 == "FP"])*as.numeric(df_summary$Freq[df_summary$Var1 == "FN"])))/
                sqrt((as.numeric(df_summary$Freq[df_summary$Var1 == "TP"])+as.numeric(df_summary$Freq[df_summary$Var1 == "FP"]))*(as.numeric(df_summary$Freq[df_summary$Var1 == "TP"])+as.numeric(df_summary$Freq[df_summary$Var1 == "FN"]))*(as.numeric(df_summary$Freq[df_summary$Var1 == "TN"])+as.numeric(df_summary$Freq[df_summary$Var1 == "FP"]))*(as.numeric(df_summary$Freq[df_summary$Var1 == "TN"])+as.numeric(df_summary$Freq[df_summary$Var1 == "FN"]))))
    demuxlet_MCC <- MCC

    ## souporcell
    eval_df_ba$souporcell_eval[str_detect(eval_df_ba$truth, "D") & eval_df_ba$truth == eval_df_ba$souporcell_sample] <- "TN"
    eval_df_ba$souporcell_eval[str_detect(eval_df_ba$truth, "C") & eval_df_ba$truth == eval_df_ba$souporcell_sample] <- "TP"
    eval_df_ba$souporcell_eval[str_detect(eval_df_ba$truth, "C") & eval_df_ba$truth != eval_df_ba$souporcell_sample] <- "FP"
    eval_df_ba$souporcell_eval[str_detect(eval_df_ba$truth, "D") & eval_df_ba$truth != eval_df_ba$souporcell_sample] <- "FN"
    df <- eval_df_ba
    df_summary <- data.frame(table(df$souporcell_eval))
    Var1 <- c("FN", "TN", "TP", "FP")
    Freq<- c(0,0,0,0)
    filler_frame <- data.frame(Var1, Freq)
    df_summary <- rbind(df_summary, filler_frame)
    df_summary <- df_summary %>% dplyr::group_by(Var1) %>% 
        dplyr::summarise(Freq = sum(Freq)) %>% as.data.frame()
    MCC <- ((as.numeric((df_summary$Freq[df_summary$Var1 == "TP"])*as.numeric(df_summary$Freq[df_summary$Var1 == "TN"]))-(as.numeric(df_summary$Freq[df_summary$Var1 == "FP"])*as.numeric(df_summary$Freq[df_summary$Var1 == "FN"])))/
                sqrt((as.numeric(df_summary$Freq[df_summary$Var1 == "TP"])+as.numeric(df_summary$Freq[df_summary$Var1 == "FP"]))*(as.numeric(df_summary$Freq[df_summary$Var1 == "TP"])+as.numeric(df_summary$Freq[df_summary$Var1 == "FN"]))*(as.numeric(df_summary$Freq[df_summary$Var1 == "TN"])+as.numeric(df_summary$Freq[df_summary$Var1 == "FP"]))*(as.numeric(df_summary$Freq[df_summary$Var1 == "TN"])+as.numeric(df_summary$Freq[df_summary$Var1 == "FN"]))))
    souporcell_MCC <- MCC
        
    ## Dataframe
    tool <- c("ensemblex",  "Demuxalot", "Demuxlet", "Souporcell", "Vireo")
    MCC_thresh <- c(ensemblex_MCC, demuxalot_MCC, demuxlet_MCC, souporcell_MCC, vireo_MCC)
    MCC_threshold_df <- data.frame(tool, MCC_thresh)
    MCC_threshold_df$tool <- factor(MCC_threshold_df$tool, levels = c("ensemblex",  "Demuxalot", "Demuxlet", "Souporcell", "Vireo"))
##

#############################
# AUC for singlet detection #
#############################
## code
    eval_df_ba <- read.delim('~/confidence_PWE_GBD_EID_droplet_assignment.csv', header = T, sep = ",") 

    ### Add eval columns
    ## Vireo 
    eval_df_ba$vireo_eval <- "bad"
    eval_df_ba$vireo_eval[eval_df_ba$vireo == eval_df_ba$truth ] <- "good"

    ## Souporcell 
    eval_df_ba$souporcell_eval <- "bad"
    eval_df_ba$souporcell_eval[eval_df_ba$souporcell == eval_df_ba$truth ] <- "good"

    ## Demuxlet 
    eval_df_ba$demuxlet_eval <- "bad"
    eval_df_ba$demuxlet_eval[eval_df_ba$freemuxlet == eval_df_ba$truth ] <- "good"

    ## Demuxalot
    eval_df_ba$demuxalot_eval <- "bad"
    eval_df_ba$demuxalot_eval[eval_df_ba$demuxalot == eval_df_ba$truth ] <- "good"

    ## ensemblex 
    eval_df_ba$concensus_eval <- "bad"
    eval_df_ba$concensus_eval[eval_df_ba$assignment_R_consensus == eval_df_ba$truth ] <- "good"

    ## Subset main dataframe to only include true singlets
    df <- subset(eval_df_ba, truth != "Doublet")
    df <- subset(df, truth != "Negative")

    ## Vireo 
    roc_empirical <- rocit(score = df$vireo_max, class = df$vireo_eval, ##change_here
                    negref = "bad") 
    summary(roc_empirical)
    vireo_AUC <- roc_empirical$AUC  

    ## Demuxlet 
    roc_empirical <- rocit(score = df$freemuxlet_max_posterior, class = df$demuxlet_eval, ##change_here
                    negref = "bad") 
    summary(roc_empirical)
    freemuxlet_AUC <- roc_empirical$AUC  

    ## Demuxalot
    roc_empirical <- rocit(score = df$demuxalot_max, class = df$demuxalot_eval, ##change_here
                    negref = "bad") 
    summary(roc_empirical)
    AUC_demuxalot <- roc_empirical$AUC

    ## Souporcell
    roc_empirical <- rocit(score =  df$souporcell_log_prob_singleton, class = df$souporcell_eval, ##change_here
                    negref = "good") 
    summary(roc_empirical)
    AUC_souporcell <-roc_empirical$AUC

    ## ensemblex
    roc_empirical <- rocit(score = df$singlet_confidence, class = df$concensus_eval, ##change_here
                    negref = "bad") 
    summary(roc_empirical)
    AUC_ensemblex <-roc_empirical$AUC
    AUC_ensemblex
##