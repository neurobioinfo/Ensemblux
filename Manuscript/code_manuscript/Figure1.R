#!/usr/bin/env Rscript

#########
# NOTES #
#########
# This code was used to produce Figure 1 of the Ensemblux manuscript

########
# Main #
########
## Load libraries
packages <- c('dplyr', 'tidyr' 'pdfCluster', 'data.table','readr','lubridate', 'tidyverse', 'moments', 'mousetrap', 'usethis', 'devtools', 'desc', 'kneedle', 'ROCit', 'ggplot2', 'factoextra', 'ggpubr', 'ComplexUpset', 'Seurat', 'Matrix', 'scCustomize')
lapply(packages, library, character.only = TRUE)

## Figure 1 function
# eval_df = output of Ensemblux
# run = simulated pool ID
# sample_size = number of samples pooled
# concentration = number of cells pooled
figure_1 <- function(eval_df, run, sample_size, concentration){
    #####################
    ## Proportion correct all cells 
    #####################
    total_cell <- nrow(eval_df)
    
    ## Vireo 
    eval_df$vireo_eval <- "bad"
    eval_df$vireo_eval[eval_df$vireo_sample == eval_df$truth ] <- "good"
    
    ## Souporcell 
    eval_df$souporcell_eval <- "bad"
    eval_df$souporcell_eval[eval_df$souporcell_sample == eval_df$truth ] <- "good"

    ## Demuxlet 
    eval_df$freemuxlet_eval <- "bad"
    eval_df$freemuxlet_eval[eval_df$freemuxlet_sample == eval_df$truth ] <- "good"
    
    ## Demuxalot
    eval_df$demuxalot_eval <- "bad"
    eval_df$demuxalot_eval[eval_df$demuxalot_sample == eval_df$truth ] <- "good"
    
    bind_lim <-dplyr::select(eval_df, c("vireo_eval", "souporcell_eval", "freemuxlet_eval","demuxalot_eval"))

    data_long <- gather(bind_lim, tool, eval, vireo_eval:demuxalot_eval, factor_key=TRUE)

    test <- data_long %>% 
        dplyr::group_by(tool, eval ) %>% 
        dplyr::filter() %>% 
        dplyr::count() %>%
        dplyr::rename("Variant_N"= "n")

    test_good <- subset(test, eval == "good")
    test_good$prop <- test_good$Variant_N/total_cell
    test_good <- dplyr::select(test_good, "tool", "prop")
    test_good$run <- as.character(run)
    test_good$sample_size <- as.character(sample_size) 
    test_good$concentration <- as.character(concentration) 
    test_good$droplet_type <- "all"
    test_good
    prop_tool_all_cells <- test_good

    #####################
    ## Proportion correct singlets
    #####################
    eval_df_singlet <- subset(eval_df, truth != "doublet")
    total_cell <- nrow(eval_df_singlet)
    
    ## Vireo 
    eval_df_singlet$vireo_eval <- "bad"
    eval_df_singlet$vireo_eval[eval_df_singlet$vireo_sample == eval_df_singlet$truth ] <- "good"
    
    ## Souporcell 
    eval_df_singlet$souporcell_eval <- "bad"
    eval_df_singlet$souporcell_eval[eval_df_singlet$souporcell_sample == eval_df_singlet$truth ] <- "good"
    
    ## Demuxlet 
    eval_df_singlet$freemuxlet_eval <- "bad"
    eval_df_singlet$freemuxlet_eval[eval_df_singlet$freemuxlet_sample == eval_df_singlet$truth ] <- "good"
    
    ## Demuxalot
    eval_df_singlet$demuxalot_eval <- "bad"
    eval_df_singlet$demuxalot_eval[eval_df_singlet$demuxalot_sample == eval_df_singlet$truth ] <- "good"
    
    bind_lim <-dplyr::select(eval_df_singlet, c("vireo_eval", "souporcell_eval", "freemuxlet_eval","demuxalot_eval"))
    
    data_long <- gather(bind_lim, tool, eval, vireo_eval:demuxalot_eval, factor_key=TRUE)
    
    test <- data_long %>% 
        dplyr::group_by(tool, eval ) %>% 
        dplyr::filter() %>% 
        dplyr::count() %>%
        dplyr::rename("Variant_N"= "n")
    
    test_good <- subset(test, eval == "good")
    test_good$prop <- test_good$Variant_N/total_cell
    test_good <- dplyr::select(test_good, "tool", "prop")
    test_good$run <- as.character(run)
    test_good$sample_size <- as.character(sample_size) 
    test_good$concentration <- as.character(concentration) 
    test_good$droplet_type <- "singlet"
    test_good
    prop_tool_singlets <- test_good

    ##################### 
    ## Proportion correct doublets
    #####################
    eval_df_doublet <- subset(eval_df, truth == "doublet")
    total_cell <- nrow(eval_df_doublet)
    
    ## Vireo 
    eval_df_doublet$vireo_eval <- "bad"
    eval_df_doublet$vireo_eval[eval_df_doublet$vireo_second.x == eval_df_doublet$truth ] <- "good"
    
    ## Souporcell 
    eval_df_doublet$souporcell_eval <- "bad"
    eval_df_doublet$souporcell_eval[eval_df_doublet$souporcell_second.x == eval_df_doublet$truth ] <- "good"

    ## Demuxlet 
    eval_df_doublet$freemuxlet_eval <- "bad"
    eval_df_doublet$freemuxlet_eval[eval_df_doublet$freemuxlet_second.x == eval_df_doublet$truth ] <- "good"
    
    ## Demuxalot
    eval_df_doublet$demuxalot_eval <- "bad"
    eval_df_doublet$demuxalot_eval[eval_df_doublet$demuxalot_sample.x == eval_df_doublet$truth ] <- "good"

    bind_lim <-dplyr::select(eval_df_doublet, c("vireo_eval", "souporcell_eval", "freemuxlet_eval","demuxalot_eval"))
    
    data_long <- gather(bind_lim, tool, eval, vireo_eval:demuxalot_eval, factor_key=TRUE)
        
    test <- data_long %>% 
        dplyr::group_by(tool, eval ) %>% 
        dplyr::filter() %>% 
        dplyr::count() %>%
        dplyr::rename("Variant_N"= "n")
    
    test_good <- subset(test, eval == "good")
    test_good$prop <- test_good$Variant_N/total_cell
    test_good <- dplyr::select(test_good, "tool", "prop")
    test_good$run <- as.character(run)
    test_good$sample_size <- as.character(sample_size) 
    test_good$concentration <- as.character(concentration) 
    test_good$droplet_type <- "doublet"
    test_good
    prop_tool_doublets <- test_good

    ##################### 
    ## Proportion of all cells correctly classified by at least one tool
    #####################  
    run_frame <- eval_df
    
    number_droplets <- nrow(run_frame)
    number_doublets <- nrow(run_frame[run_frame$truth == "doublet",])
    number_singlets <- nrow(run_frame[run_frame$truth != "doublet",])
    
    ## Vireo 
    run_frame$vireo_eval <- "bad"
    run_frame$vireo_eval[run_frame$vireo_second.x == run_frame$truth ] <- "good"
    
    ## Souporcell 
    run_frame$souporcell_eval <- "bad"
    run_frame$souporcell_eval[run_frame$souporcell_second.x == run_frame$truth ] <- "good"
    
    ## Demuxlet 
    run_frame$freemuxlet_eval <- "bad"
    run_frame$freemuxlet_eval[run_frame$freemuxlet_second.x == run_frame$truth ] <- "good"
    
    ## Demuxalot
    run_frame$demuxalot_eval <- "bad"
    run_frame$demuxalot_eval[run_frame$demuxalot_sample.x == run_frame$truth ] <- "good"
        
    ## correct classification: all droplets
    vireo_list <- run_frame$barcode[run_frame$vireo_eval == "good"]
    length(vireo_list)
    souporcell_list <- run_frame$barcode[run_frame$souporcell_eval == "good"] 
    length(souporcell_list)
    freemuxlet_list <- run_frame$barcode[run_frame$freemuxlet_eval == "good"]
    length(freemuxlet_list)
    demuxalot_list <- run_frame$barcode[run_frame$demuxalot_eval == "good"]
    length(demuxalot_list)
    
    #create list
    total_list_upset <- list(  vireo_barcodes =vireo_list, 
                                souporcell_barcodes =souporcell_list,
                                freemuxlet_barcodes = freemuxlet_list,
                                demuxalot_barcodes = demuxalot_list)
    
    names(total_list_upset)= c(
        "Vireo",
        "Souporcell",
        "Freemuxlet",
        "Demuxalot",
        "Concensus")
    
    Unique_genes<- unique(unlist(total_list_upset))
    df <- data.frame(genes = Unique_genes,
                    Vireo.bin = 0,
                    Souporcell.bin = 0,
                    Freemuxlet.bin = 0,
                    Demuxalot.bin = 0)
    
    ## Vireo
    for (i in seq_along(df$genes)) {
        if (df$genes[i] %in% vireo_list) {
        df$Vireo.bin[i] <- 1
        } else {
        df$Vireo.bin[i] <- 0
        }
    }
    
    ## Souporcell
    for (i in seq_along(df$genes)) {
        if (df$genes[i] %in% souporcell_list) {
        df$Souporcell.bin[i] <- 1
        } else {
        df$Souporcell.bin[i] <- 0
        }
    }
    
    ## Demuxlet
    for (i in seq_along(df$genes)) {
        if (df$genes[i] %in% freemuxlet_list) {
        df$Freemuxlet.bin[i] <- 1
        } else {
        df$Freemuxlet.bin[i] <- 0
        }
    }
    
    ## Demuxalot
    for (i in seq_along(df$genes)) {
        if (df$genes[i] %in% demuxalot_list) {
        df$Demuxalot.bin[i] <- 1
        } else {
        df$Demuxalot.bin[i] <- 0
        }
    }
    
    colnames(df) <- c("genes", 
                        "Vireo", 
                        "Souporcell", 
                        "Freemuxlet",
                        "Demuxalot")
    
    prop_all <- nrow(df)/number_droplets
    

    ##################### 
    ## Proportion of singlets correctly classified by at least one tool
    #####################  
    run_frame <- eval_df
    
    number_droplets <- nrow(run_frame)
    number_doublets <- nrow(run_frame[run_frame$truth == "doublet",])
    number_singlets <- nrow(run_frame[run_frame$truth != "doublet",])
    run_frame <- subset(run_frame, truth != "doublet")

    ## Vireo 
    run_frame$vireo_eval <- "bad"
    run_frame$vireo_eval[run_frame$vireo_second.x == run_frame$truth ] <- "good"
    
    ## Souporcell 
    run_frame$souporcell_eval <- "bad"
    run_frame$souporcell_eval[run_frame$souporcell_second.x == run_frame$truth ] <- "good"
    
    ## Demuxlet 
    run_frame$freemuxlet_eval <- "bad"
    run_frame$freemuxlet_eval[run_frame$freemuxlet_second.x == run_frame$truth ] <- "good"
    
    ## Demuxalot
    run_frame$demuxalot_eval <- "bad"
    run_frame$demuxalot_eval[run_frame$demuxalot_sample.x == run_frame$truth ] <- "good"
    
    ## Correct classification: singlets
    vireo_list <- run_frame$barcode[run_frame$vireo_eval == "good"]
    length(vireo_list)
    souporcell_list <- run_frame$barcode[run_frame$souporcell_eval == "good"] 
    length(souporcell_list)
    freemuxlet_list <- run_frame$barcode[run_frame$freemuxlet_eval == "good"]
    length(freemuxlet_list)
    demuxalot_list <- run_frame$barcode[run_frame$demuxalot_eval == "good"]
    length(demuxalot_list)
    
    #create list
    total_list_upset <- list(  vireo_barcodes =vireo_list, 
                                souporcell_barcodes =souporcell_list,
                                freemuxlet_barcodes = freemuxlet_list,
                                demuxalot_barcodes = demuxalot_list)
    
    names(total_list_upset)= c(
        "Vireo",
        "Souporcell",
        "Freemuxlet",
        "Demuxalot")
    
    Unique_genes<- unique(unlist(total_list_upset))
    df <- data.frame(genes = Unique_genes,
                    Vireo.bin = 0,
                    Souporcell.bin = 0,
                    Freemuxlet.bin = 0,
                    Demuxalot.bin = 0)
    
    ## Vireo
    for (i in seq_along(df$genes)) {
        if (df$genes[i] %in% vireo_list) {
        df$Vireo.bin[i] <- 1
        } else {
        df$Vireo.bin[i] <- 0
        }
    }
    
    ## Souporcell
    for (i in seq_along(df$genes)) {
        if (df$genes[i] %in% souporcell_list) {
        df$Souporcell.bin[i] <- 1
        } else {
        df$Souporcell.bin[i] <- 0
        }
    }
    
    ## Demuxlet
    for (i in seq_along(df$genes)) {
        if (df$genes[i] %in% freemuxlet_list) {
        df$Freemuxlet.bin[i] <- 1
        } else {
        df$Freemuxlet.bin[i] <- 0
        }
    }
    
    ## Demuxalot
    for (i in seq_along(df$genes)) {
        if (df$genes[i] %in% demuxalot_list) {
        df$Demuxalot.bin[i] <- 1
        } else {
        df$Demuxalot.bin[i] <- 0
        }
    }
    
    colnames(df) <- c("genes", 
                        "Vireo", 
                        "Souporcell", 
                        "Freemuxlet",
                        "Demuxalot")
    
    prop_singlets <- nrow(df)/number_singlets
    

    ##################### 
    ## Proportion of doublets correctly classified by at least one tool
    ##################### 
    run_frame <- eval_df
    number_droplets <- nrow(run_frame)
    number_doublets <- nrow(run_frame[run_frame$truth == "doublet",])
    number_singlets <- nrow(run_frame[run_frame$truth != "doublet",])
    run_frame <- subset(run_frame, truth == "doublet")
    
    ## Vireo 
    run_frame$vireo_eval <- "bad"
    run_frame$vireo_eval[run_frame$vireo_second.x == run_frame$truth ] <- "good"
    
    ## Souporcell 
    run_frame$souporcell_eval <- "bad"
    run_frame$souporcell_eval[run_frame$souporcell_second.x == run_frame$truth ] <- "good"
    
    ## Demuxlet 
    run_frame$freemuxlet_eval <- "bad"
    run_frame$freemuxlet_eval[run_frame$freemuxlet_second.x == run_frame$truth ] <- "good"
    
    ## Demuxalot
    run_frame$demuxalot_eval <- "bad"
    run_frame$demuxalot_eval[run_frame$demuxalot_sample.x == run_frame$truth ] <- "good"
    
    ## correct classification: all droplets
    vireo_list <- run_frame$barcode[run_frame$vireo_eval == "good"]
    length(vireo_list)
    souporcell_list <- run_frame$barcode[run_frame$souporcell_eval == "good"] 
    length(souporcell_list)
    freemuxlet_list <- run_frame$barcode[run_frame$freemuxlet_eval == "good"]
    length(freemuxlet_list)
    demuxalot_list <- run_frame$barcode[run_frame$demuxalot_eval == "good"]
    length(demuxalot_list)
    
    #create list
    total_list_upset <- list(  vireo_barcodes =vireo_list, 
                                souporcell_barcodes =souporcell_list,
                                freemuxlet_barcodes = freemuxlet_list,
                                demuxalot_barcodes = demuxalot_list)
    
    names(total_list_upset)= c(
        "Vireo",
        "Souporcell",
        "Freemuxlet",
        "Demuxalot")
    
    
    #list_test <- as.data.frame(total_list_upset)
    Unique_genes<- unique(unlist(total_list_upset))
    df <- data.frame(genes = Unique_genes,
                    Vireo.bin = 0,
                    Souporcell.bin = 0,
                    Freemuxlet.bin = 0,
                    Demuxalot.bin = 0)
    
    
    #Vireo
    for (i in seq_along(df$genes)) {
        if (df$genes[i] %in% vireo_list) {
        df$Vireo.bin[i] <- 1
        } else {
        df$Vireo.bin[i] <- 0
        }
    }
    
    #Souporcell
    for (i in seq_along(df$genes)) {
        if (df$genes[i] %in% souporcell_list) {
        df$Souporcell.bin[i] <- 1
        } else {
        df$Souporcell.bin[i] <- 0
        }
    }
    
    #freemuxlet
    for (i in seq_along(df$genes)) {
        if (df$genes[i] %in% freemuxlet_list) {
        df$Freemuxlet.bin[i] <- 1
        } else {
        df$Freemuxlet.bin[i] <- 0
        }
    }
    
    #demuxalot
    for (i in seq_along(df$genes)) {
        if (df$genes[i] %in% demuxalot_list) {
        df$Demuxalot.bin[i] <- 1
        } else {
        df$Demuxalot.bin[i] <- 0
        }
    }
    
    colnames(df) <- c("genes", 
                        "Vireo", 
                        "Souporcell",
                        "Freemuxlet",
                        "Demuxalot")
    
    prop_doublets <- nrow(df)/number_doublets

    ##################### 
    #Proportion of droplets uniquely correctl classsified by each tool
    ##################### 
    number_cells <- nrow(eval_df)

    ## Vireo 
    eval_df$vireo_eval <- "bad"
    eval_df$vireo_eval[eval_df$vireo_second.x == eval_df$truth ] <- "good"

    ## Souporcell 
    eval_df$souporcell_eval <- "bad"
    eval_df$souporcell_eval[eval_df$souporcell_second.x == eval_df$truth ] <- "good"

    ## Demuxlet 
    eval_df$freemuxlet_eval <- "bad"
    eval_df$freemuxlet_eval[eval_df$freemuxlet_second.x == eval_df$truth ] <- "good"

    ## Demuxalot
    eval_df$demuxalot_eval <- "bad"
    eval_df$demuxalot_eval[eval_df$demuxalot_sample.x == eval_df$truth ] <- "good"

    vireo_list <- eval_df$barcode[eval_df$vireo_eval == "good"]
    length(vireo_list)
    souporcell_list <- eval_df$barcode[eval_df$souporcell_eval == "good"] 
    length(souporcell_list)
    freemuxlet_list <- eval_df$barcode[eval_df$freemuxlet_eval == "good"]
    length(freemuxlet_list)
    demuxalot_list <- eval_df$barcode[eval_df$demuxalot_eval == "good"]
    length(demuxalot_list)

    #create list
    total_list_upset <- list(  vireo_barcodes =vireo_list, 
                            souporcell_barcodes =souporcell_list,
                            freemuxlet_barcodes = freemuxlet_list,
                            demuxalot_barcodes = demuxalot_list)
    length(total_list_upset)

    names(total_list_upset)= c(
    "Vireo",
    "Souporcell",
    "Freemuxlet",
    "Demuxalot")

    Unique_genes<- unique(unlist(total_list_upset))
    length(unique(Unique_genes))


    df <- data.frame(genes = Unique_genes,
                    Vireo.bin = 0,
                    Souporcell.bin = 0,
                    Freemuxlet.bin = 0,
                    Demuxalot.bin = 0)


    ## Vireo
    df$Vireo.bin[df$genes %in% vireo_list] <- 1
    unique(df$Vireo.bin)
    ## Souporcell
    df$Souporcell.bin[df$genes %in% souporcell_list] <- 1
    unique(df$Souporcell.bin)
    ## Demuxlet
    df$Freemuxlet.bin[df$genes %in% freemuxlet_list] <- 1
    unique(df$Freemuxlet.bin)
    ## Demuxalot
    df$Demuxalot.bin[df$genes %in% demuxalot_list] <- 1
    unique(df$Demuxalot.bin)

    colnames(df) <- c("genes", 
                    "Vireo", 
                    "Souporcell",
                    "Freemuxlet",
                    "Demuxalot")

    ## Demuxalot
    demuxalot_unique <- df$genes[df$Demuxalot == 1 & df$Vireo == 0 & df$Souporcell == 0  & df$Freemuxlet == 0]
    demuxalot_unique <- length(demuxalot_unique)
    ## Demuxlet
    freemuxlet_unique <- df$genes[df$Demuxalot == 0 & df$Vireo == 0 & df$Souporcell == 0  & df$Freemuxlet == 1]
    freemuxlet_unique <- length(freemuxlet_unique)
    ## Vireo
    vireo_unique <- df$genes[df$Demuxalot == 0 & df$Vireo == 1 & df$Souporcell == 0  & df$Freemuxlet == 0]
    vireo_unique <- length(vireo_unique)
    ## Souporcell
    souporcell_unique <- df$genes[df$Demuxalot == 0 & df$Vireo == 0 & df$Souporcell == 1  & df$Freemuxlet == 0]
    souporcell_unique <- length(souporcell_unique)

    ## calculate the prop of all cells
    demuxalot_unique <- demuxalot_unique/number_cells
    freemuxlet_unique <- freemuxlet_unique/number_cells
    vireo_unique <- vireo_unique/number_cells
    souporcell_unique <- souporcell_unique/number_cells

    ## prepare a dataframe for downstream merge
    tool <- c('vireo_eval', 'souporcell_eval', 'freemuxlet_eval', 'demuxalot_eval')
    prop_all_unique_identified <- c(vireo_unique, souporcell_unique, freemuxlet_unique, demuxalot_unique)
    prop_all_unique_identified_df <- data.frame(tool, prop_all_unique_identified)


    ##################### 
    #Proportion of cells correctly classified but labelled as unassigned due insufficient assignment probabilities
    ##################### 
    eval_df_sing <- subset(eval_df, truth != "doublet")
    number_cells <- nrow(eval_df_sing)

    ### No threshold
    ## Vireo 
    eval_df_sing$vireo_eval_no_thresh <- "bad"
    eval_df_sing$vireo_eval_no_thresh[eval_df_sing$vireo_second.x == eval_df_sing$truth ] <- "good"

    ## Souporcell 
    eval_df_sing$souporcell_eval_no_thresh <- "bad"
    eval_df_sing$souporcell_eval_no_thresh[eval_df_sing$souporcell_second.x == eval_df_sing$truth ] <- "good"

    ## Demuxlet 
    eval_df_sing$freemuxlet_eval_no_thresh <- "bad"
    eval_df_sing$freemuxlet_eval_no_thresh[eval_df_sing$freemuxlet_second.x == eval_df_sing$truth ] <- "good"

    ## Demuxalot
    eval_df_sing$demuxalot_eval_no_thresh <- "bad"
    eval_df_sing$demuxalot_eval_no_thresh[eval_df_sing$demuxalot_sample.x == eval_df_sing$truth ] <- "good"


    ### Threshold
    ## Vireo 
    eval_df_sing$vireo_eval_thresh <- "bad"
    eval_df_sing$vireo_eval_thresh[eval_df_sing$vireo_sample == eval_df_sing$truth ] <- "good"
    
    ## Souporcell 
    eval_df_sing$souporcell_eval_thresh <- "bad"
    eval_df_sing$souporcell_eval_thresh[eval_df_sing$souporcell_sample == eval_df_sing$truth ] <- "good"
    
    ## Demuxlet 
    eval_df_sing$freemuxlet_eval_thresh <- "bad"
    eval_df_sing$freemuxlet_eval_thresh[eval_df_sing$freemuxlet_sample == eval_df_sing$truth ] <- "good"
    
    ## Demuxalot
    eval_df_sing$demuxalot_eval_thresh <- "bad"
    eval_df_sing$demuxalot_eval_thresh[eval_df_sing$demuxalot_sample == eval_df_sing$truth ] <- "good"

    bind_lim <-dplyr::select(eval_df_sing, c("vireo_eval_no_thresh", "souporcell_eval_no_thresh", "freemuxlet_eval_no_thresh","demuxalot_eval_no_thresh"))
    data_long <- gather(bind_lim, tool, eval, vireo_eval_no_thresh:demuxalot_eval_no_thresh, factor_key=TRUE)

    test_no_thresh <- data_long %>% 
        dplyr::group_by(tool, eval ) %>% 
        dplyr::filter() %>% 
        dplyr::count() %>%
        dplyr::rename("no_thresh"= "n")

        test_no_thresh <- subset(test_no_thresh, eval == "good")
        test_no_thresh$tool_name <- c("vireo", "souporcell", "freemuxlet", "demuxalot")

    bind_lim <-dplyr::select(eval_df_sing, c("vireo_eval_thresh", "souporcell_eval_thresh", "freemuxlet_eval_thresh","demuxalot_eval_thresh"))
    data_long <- gather(bind_lim, tool, eval, vireo_eval_thresh:demuxalot_eval_thresh, factor_key=TRUE)
        
    test_thresh <- data_long %>% 
        dplyr::group_by(tool, eval ) %>% 
        dplyr::filter() %>% 
        dplyr::count() %>%
        dplyr::rename("thresh"= "n")

        test_thresh <- subset(test_thresh, eval == "good")
        test_thresh$tool_name <- c("vireo", "souporcell", "freemuxlet", "demuxalot")

    merge_improve <- merge(test_no_thresh, test_thresh, by = "tool_name")
    merge_improve$improve <- merge_improve$no_thresh - merge_improve$thresh
    merge_improve$prop_improve <- merge_improve$improve/number_cells


    ##################### 
    #merge into summary frame for downstream plotting with all datframes.
    ##################### 

    # prop correct all cells
    prop_tool_all_cells <- prop_tool_all_cells %>% dplyr::select(tool, prop)
    prop_tool_all_cells$tool_name <- c("vireo", "souporcell", "freemuxlet", "demuxalot")
    prop_tool_all_cells <- prop_tool_all_cells[,c(3,4)]
    colnames(prop_tool_all_cells) <- c("prop_all_threshold", "tool_name")

    # prop correct singlets
    prop_tool_singlets <- prop_tool_singlets %>% dplyr::select(tool, prop)
    prop_tool_singlets$tool_name <- c("vireo", "souporcell","freemuxlet", "demuxalot")
    prop_tool_singlets <- prop_tool_singlets[,c(3,4)]
    colnames(prop_tool_singlets) <- c("prop_singlet_threshold", "tool_name")

    # prop correct doublets
    prop_tool_doublets <- prop_tool_doublets %>% dplyr::select(tool, prop)
    prop_tool_doublets$tool_name <- c("vireo", "souporcell", "freemuxlet", "demuxalot")
    prop_tool_doublets <- prop_tool_doublets[,c(3,4)]
    colnames(prop_tool_doublets) <- c("prop_doublet_no_threshold", "tool_name")

    # prop correct all at least one
    prop_all_threshold <- prop_all
    tool_name <- "at_least_one"
    at_least_one_df_all <- data.frame(prop_all_threshold, tool_name)
    prop_tool_all_cells <- rbind(prop_tool_all_cells, at_least_one_df_all)

    # prop correct singlet at least one
    prop_singlet_threshold <- prop_singlets
    tool_name <- "at_least_one"
    at_least_one_df_singlet <- data.frame(prop_singlet_threshold, tool_name)
    prop_tool_singlets <- rbind(prop_tool_singlets, at_least_one_df_singlet)

    # prop correct doublet at least one
    prop_doublet_no_threshold <- prop_doublets
    tool_name <- "at_least_one"
    at_least_one_df_doublets <- data.frame(prop_doublet_no_threshold, tool_name)
    prop_tool_doublets <- rbind(prop_tool_doublets, at_least_one_df_doublets)

    # unique identified by each tool
    prop_all_unique_identified_df$tool_name <- c("vireo", "souporcell", "freemuxlet", "demuxalot")
    prop_all_unique_identified_df <- prop_all_unique_identified_df[,c(2,3)]

    # no threshold improvement
    merge_improve <- merge_improve %>% dplyr::select(tool_name, prop_improve)

    ## merge all dataframe
    eval_merge <- merge(prop_tool_all_cells, prop_tool_singlets, by = "tool_name")
    eval_merge <- merge(eval_merge, prop_tool_doublets, by = "tool_name")
    eval_merge <- merge(eval_merge, prop_all_unique_identified_df, by = "tool_name", all = T)
    eval_merge <- merge(eval_merge, merge_improve, by = "tool_name", all = T)

    ## add run information 
    eval_merge$run <- as.character(run)
    eval_merge$sample_size <- sample_size 
    eval_merge$concentration <- concentration 
    eval_merge
}

