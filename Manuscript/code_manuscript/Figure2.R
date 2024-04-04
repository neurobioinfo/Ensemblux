#!/usr/bin/env Rscript

#########
# NOTES #
#########
# This code was used to produce Figure 2 of the Ensemblux manuscript

########
# Main #
########
## Load libraries
packages <- c('dplyr', 'tidyr' 'pdfCluster', 'data.table','readr','lubridate', 'tidyverse', 'moments', 'mousetrap', 'usethis', 'devtools', 'desc', 'kneedle', 'ROCit', 'ggplot2', 'factoextra', 'ggpubr', 'ComplexUpset', 'Seurat', 'Matrix', 'scCustomize' )
lapply(packages, library, character.only = TRUE)

## Figure 2 function
# PWE = Ensemblux probabilistic-weighted ensemble output
# GBE = Ensemblux graph-based doubet detection output
# EID = Ensemblux ebsemble-independent doublet detection output
# run_code = simulated pool ID
# size = number of samples pooled
figure_2 <- function(PWE, GBD, EID, run_code, size){
    ### evaluate PWE
        ## singlet
        PWE_singlet <- subset(PWE, truth != "doublet")
        total_cell <- nrow(PWE_singlet)
        PWE_singlet <- PWE_singlet[!duplicated(PWE_singlet$barcode), ]
        length(unique(PWE_singlet$barcode))
        PWE_singlet$ensemblux_eval <- "bad"
        PWE_singlet$ensemblux_eval[PWE_singlet$assignment_R_consensus == PWE_singlet$truth ] <- "good"
        test <- data.frame(table(PWE_singlet$ensemblux_eval))
        prop <- test$Freq[test$Var1 == "good"]
        prop_PWE_singlet <- prop/total_cell

        ## doublet
        PWE_doublet <- subset(PWE, truth == "doublet")
        total_cell <- nrow(PWE_doublet)
        PWE_doublet <- PWE_doublet[!duplicated(PWE_doublet$barcode), ]
        length(unique(PWE_doublet$barcode))
        PWE_doublet$ensemblux_eval <- "bad"
        PWE_doublet$ensemblux_eval[PWE_doublet$assignment_R_consensus == PWE_doublet$truth ] <- "good"
        test <- data.frame(table(PWE_doublet$ensemblux_eval))
        prop <- test$Freq[test$Var1 == "good"]
        prop_PWE_doublet <- prop/total_cell

    ### evaluate GBD
        ## singlet
        GBD_singlet <- subset(GBD, truth != "doublet")
        total_cell <- nrow(GBD_singlet)
        GBD_singlet <- GBD_singlet[!duplicated(GBD_singlet$barcode), ]
        length(unique(GBD_singlet$barcode))
        GBD_singlet$ensemblux_eval <- "bad"
        GBD_singlet$ensemblux_eval[GBD_singlet$assignment_R_consensus == GBD_singlet$truth ] <- "good"
        test <- data.frame(table(GBD_singlet$ensemblux_eval))
        prop <- test$Freq[test$Var1 == "good"]
        prop_GBD_singlet <- prop/total_cell

        ## doublet
        GBD_doublet <- subset(GBD, truth == "doublet")
        total_cell <- nrow(GBD_doublet)
        GBD_doublet <- GBD_doublet[!duplicated(GBD_doublet$barcode), ]
        length(unique(GBD_doublet$barcode))
        GBD_doublet$ensemblux_eval <- "bad"
        GBD_doublet$ensemblux_eval[GBD_doublet$assignment_R_consensus == GBD_doublet$truth ] <- "good"
        test <- data.frame(table(GBD_doublet$ensemblux_eval))
        prop <- test$Freq[test$Var1 == "good"]
        prop_GBD_doublet <- prop/total_cell

    ### evaluate EID
        ## singlet
        EID_singlet <- subset(EID, truth != "doublet")
        total_cell <- nrow(EID_singlet)
        EID_singlet <- EID_singlet[!duplicated(EID_singlet$barcode), ]
        length(unique(EID_singlet$barcode))
        EID_singlet$ensemblux_eval <- "bad"
        EID_singlet$ensemblux_eval[EID_singlet$assignment_R_consensus == EID_singlet$truth ] <- "good"
        test <- data.frame(table(EID_singlet$ensemblux_eval))
        prop <- test$Freq[test$Var1 == "good"]
        prop_EID_singlet <- prop/total_cell

        ## doublet
        EID_doublet <- subset(EID, truth == "doublet")
        total_cell <- nrow(EID_doublet)
        EID_doublet <- EID_doublet[!duplicated(EID_doublet$barcode), ]
        length(unique(EID_doublet$barcode))
        EID_doublet$ensemblux_eval <- "bad"
        EID_doublet$ensemblux_eval[EID_doublet$assignment_R_consensus == EID_doublet$truth ] <- "good"
        test <- data.frame(table(EID_doublet$ensemblux_eval))
        prop <- test$Freq[test$Var1 == "good"]
        prop_EID_doublet <- prop/total_cell

    ### make a summary dataframe
    tool_name <- c("PWE", "GBD", "EID")
    prop_singlet_threshold <- c(prop_PWE_singlet, prop_GBD_singlet, prop_EID_singlet)
    prop_doublet_no_threshold <- c(prop_PWE_doublet, prop_GBD_doublet, prop_EID_doublet)
    run <- c(run_code, run_code, run_code)
    sample_size <- c(size, size, size)
    sum_df <- data.frame(tool_name, prop_singlet_threshold, prop_doublet_no_threshold, run, sample_size)
    sum_df
}