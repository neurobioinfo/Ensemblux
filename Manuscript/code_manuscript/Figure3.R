#!/usr/bin/env Rscript

#########
# NOTES #
#########
# This code was used to produce Figure 3 of the Ensemblux manuscript

########
# Main #
########
## Load libraries
packages <- c('dplyr', 'tidyr' 'pdfCluster', 'data.table','readr','lubridate', 'tidyverse', 'moments', 'mousetrap', 'usethis', 'devtools', 'desc', 'kneedle', 'ROCit', 'ggplot2', 'factoextra', 'ggpubr', 'ComplexUpset', 'Seurat', 'Matrix', 'scCustomize' )
lapply(packages, library, character.only = TRUE)

## Figure 2 function
# eval_df = Ensemblux final output
figure_3 <- function(eval_df){

        #######################################
        # prop singlet correct with threshold # 
        #######################################
        test_prop <- eval_df
        
        ## Demuxalot
        test_prop$demuxalot_test <- test_prop$demuxalot
        test_prop$demuxalot_test[test_prop$demuxalot_max < 0.9] <- "unassigned"
        prop_demuxalot  <- subset(test_prop, test_prop$truth == test_prop$demuxalot_test & test_prop$truth !="doublet") 
        test_prop_singlet <- subset(test_prop, truth != "doublet")  
        prop_singlet_demuxalot <- nrow(prop_demuxalot)/nrow(test_prop_singlet)
        
        ## Demuxlet
        prop_demuxlet  <- subset(test_prop, test_prop$truth == test_prop$freemuxlet_sample & test_prop$truth !="doublet") 
        test_prop_singlet <- subset(test_prop, truth != "doublet")
        prop_singlet_demuxlet <- nrow(prop_demuxlet)/nrow(test_prop_singlet)
        
        ## Vireo
        prop_vireo  <- subset(test_prop, test_prop$truth == test_prop$vireo_sample & test_prop$truth !="doublet") 
        test_prop_singlet <- subset(test_prop, truth != "doublet")
        prop_singlet_vireo <- nrow(prop_vireo)/nrow(test_prop_singlet)
        
        ## Souporcell
        prop_souporcell  <- subset(test_prop, test_prop$truth == test_prop$souporcell_sample & test_prop$truth !="doublet") 
        test_prop_singlet <- subset(test_prop, truth != "doublet")
        prop_singlet_souporcell <- nrow(prop_souporcell)/nrow(test_prop_singlet)
        
        ## Ensemblux
        test_prop$assignment_R_consensus[test_prop$singlet_confidence < 1.0] <- "unassigned"
        prop_ensemblux  <- subset(test_prop, test_prop$truth == test_prop$assignment_R_consensus & test_prop$truth !="doublet") 
        test_prop_singlet <- subset(test_prop, truth != "doublet")
        prop_singlet_ensemblux <- nrow(prop_ensemblux)/nrow(test_prop_singlet)
            
        ## Dataframe
        tool <- c("Ensemblux", "Demuxalot", "Demuxlet", "Souporcell", "Vireo")
        prop_singlet <- c(prop_singlet_ensemblux, prop_singlet_demuxalot, prop_singlet_demuxlet, prop_singlet_souporcell, prop_singlet_vireo)
        prop_singlet_df <- data.frame(tool, prop_singlet)
        prop_singlet_df$tool <- factor(prop_singlet_df$tool, levels = c("Ensemblux", "Demuxalot", "Demuxlet", "Souporcell", "Vireo"))
        
        #######################################
        # prop doublet correct no threshold # 
        #######################################
        test_prop <- eval_df 
        
        ## Demuxalot
        prop_demuxalot  <- subset(test_prop, test_prop$truth == test_prop$demuxalot & test_prop$truth =="doublet") 
        test_prop_doublet <- subset(test_prop, truth == "doublet")
        prop_doublet_demuxalot <- nrow(prop_demuxalot)/nrow(test_prop_doublet)
        
        ## demuxlet
        prop_demuxlet  <- subset(test_prop, test_prop$truth == test_prop$freemuxlet & test_prop$truth =="doublet") 
        test_prop_doublet <- subset(test_prop, truth == "doublet")
        prop_doublet_demuxlet <- nrow(prop_demuxlet)/nrow(test_prop_doublet)
        
        ## vireo
        prop_vireo  <- subset(test_prop, test_prop$truth == test_prop$vireo & test_prop$truth =="doublet") 
        test_prop_doublet <- subset(test_prop, truth == "doublet")
        prop_doublet_vireo <- nrow(prop_vireo)/nrow(test_prop_doublet)
        
        ## souporcell
        prop_souporcell <- subset(test_prop, test_prop$truth == test_prop$souporcell & test_prop$truth =="doublet") 
        test_prop_doublet <- subset(test_prop, truth == "doublet")
        prop_doublet_souporcell <- nrow(prop_souporcell)/nrow(test_prop_doublet)
        
        ## Ensemblux
        prop_ensemblux  <- subset(test_prop, test_prop$truth == test_prop$assignment_R_consensus & test_prop$truth =="doublet") 
        test_prop_doublet <- subset(test_prop, truth == "doublet") 
        prop_doublet_ensemblux <- nrow(prop_ensemblux)/nrow(test_prop_doublet)
        
        ## Dataframe
        tool <- c("Ensemblux", "Demuxalot", "Demuxlet", "Souporcell", "Vireo")
        prop_doublet <- c(prop_doublet_ensemblux, prop_doublet_demuxalot, prop_doublet_demuxlet, prop_doublet_souporcell, prop_doublet_vireo)
        prop_doublet_df <- data.frame(tool, prop_doublet)
        prop_doublet_df$tool <- factor(prop_doublet_df$tool, levels = c("Ensemblux", "Demuxalot", "Demuxlet", "Souporcell", "Vireo"))
        
    
        ####################################################################
        # prop all cells correct (singlet threshold; doublet no threshold) # 
        ####################################################################
        test_prop <- eval_df
        
        ## Demuxalot
        test_prop$demuxalot_test <- test_prop$demuxalot
        test_prop$demuxalot_test[test_prop$demuxalot_max < 0.9 & test_prop$demuxalot != "doublet"] <- "unassigned"
        prop_demuxalot  <- subset(test_prop, test_prop$truth == test_prop$demuxalot_test) 
        prop_all_demuxalot <- nrow(prop_demuxalot)/nrow(test_prop)
        
        ## Demuxlet
        test_prop$freemuxlet_sample[test_prop$freemuxlet == "doublet"] <- "doublet"
        prop_demuxlet  <- subset(test_prop, test_prop$truth == test_prop$freemuxlet_sample) 
        prop_all_demuxlet <- nrow(prop_demuxlet)/nrow(test_prop)
        
        ## Vireo
        test_prop$vireo_sample[test_prop$vireo == "doublet"] <- "doublet"
        prop_vireo  <- subset(test_prop, test_prop$truth == test_prop$vireo_sample) 
        prop_all_vireo <- nrow(prop_vireo)/nrow(test_prop)
        
        ## Souporcell
        test_prop$souporcell_sample[test_prop$souporcell == "doublet"] <- "doublet"
        prop_souporcell  <- subset(test_prop, test_prop$truth == test_prop$souporcell_sample) 
        prop_all_souporcell <- nrow(prop_souporcell)/nrow(test_prop)
        
        ## Ensemblux
        test_prop$assignment_R_consensus[test_prop$singlet_confidence < 1.0 & test_prop$assignment_R_consensus != "doublet"] <- "unassigned"
        prop_ensemblux  <- subset(test_prop, test_prop$truth == test_prop$assignment_R_consensus) 
        prop_all_ensemblux <- nrow(prop_ensemblux)/nrow(test_prop)
        
        ## Dataframe
        tool <- c("Ensemblux", "Demuxalot", "Demuxlet", "Souporcell", "Vireo")
        prop_all <- c(prop_all_ensemblux, prop_all_demuxalot, prop_all_demuxlet, prop_all_souporcell, prop_all_vireo)
        prop_all_df <- data.frame(tool, prop_all)
        prop_all_df$tool <- factor(prop_all_df$tool, levels = c("Ensemblux", "Demuxalot", "Demuxlet", "Souporcell", "Vireo"))
        
        
        ######################################################################
        # prop usable cells (singlet above assignment probability threshold) # 
        ######################################################################
        test_prop <- eval_df
        
        ## demuxalot
        test_prop$demuxalot_test <- test_prop$demuxalot
        test_prop$demuxalot_test[test_prop$demuxalot_max < 0.9] <- "unassigned"
        prop_demuxalot  <- subset(test_prop, test_prop$demuxalot_test != "unassigned" & test_prop$demuxalot_test !="doublet") 
        prop_usable_demuxalot <- nrow(prop_demuxalot)/nrow(test_prop)
        
        ## demuxlet
        prop_demuxlet  <- subset(test_prop, test_prop$freemuxlet_sample != "unassigned" & test_prop$freemuxlet_sample !="doublet") 
        prop_usable_demuxlet <- nrow(prop_demuxlet)/nrow(test_prop)
        
        ## vireo
        prop_vireo  <- subset(test_prop, test_prop$vireo_sample != "unassigned" & test_prop$vireo_sample !="doublet") 
        prop_usable_vireo <- nrow(prop_vireo)/nrow(test_prop)
        
        ## souporcell
        prop_souporcell  <- subset(test_prop, test_prop$souporcell_sample != "unassigned" & test_prop$souporcell_sample !="doublet") 
        prop_usable_souporcell <- nrow(prop_souporcell)/nrow(test_prop)
        
        ## Ensemblux
        test_prop$assignment_R_consensus[test_prop$singlet_confidence < 1.0] <- "unassigned"
        prop_ensemblux  <- subset(test_prop, test_prop$assignment_R_consensus != "unassigned" & test_prop$assignment_R_consensus !="doublet") 
        prop_usable_ensemblux <- nrow(prop_ensemblux)/nrow(test_prop)
        
        ## Dataframe
        tool <- c("Ensemblux", "Demuxalot", "Demuxlet", "Souporcell", "Vireo")
        prop_usable <- c(prop_usable_ensemblux, prop_usable_demuxalot, prop_usable_demuxlet, prop_usable_souporcell, prop_usable_vireo)
        prop_usable_df <- data.frame(tool, prop_usable)
        prop_usable_df$tool <- factor(prop_usable_df$tool, levels = c("Ensemblux", "Demuxalot", "Demuxlet", "Souporcell", "Vireo"))


        ##########################
        # Usable cell error rate # 
        ##########################
        test_prop <- eval_df

        ## Demuxalot
        test_prop$demuxalot_test <- test_prop$demuxalot
        test_prop$demuxalot_test[test_prop$demuxalot_max < 0.9] <- "unassigned"
        prop_demuxalot  <- subset(test_prop, test_prop$demuxalot_test != "unassigned" & test_prop$demuxalot_test !="doublet") 
        prop_demuxalot_correct <- subset(prop_demuxalot, demuxalot_test == truth)
        prop_usable_error_demuxalot <- 1- nrow(prop_demuxalot_correct)/nrow(prop_demuxalot)
                
        ## demuxlet
        prop_demuxlet  <- subset(test_prop, test_prop$freemuxlet_sample != "unassigned" & test_prop$freemuxlet_sample !="doublet") 
        prop_demuxlet_correct <- subset(prop_demuxlet, freemuxlet_sample == truth)
        prop_usable_error_demuxlet <- 1-nrow(prop_demuxlet_correct)/nrow(prop_demuxlet)
        
        ## vireo
        prop_vireo  <- subset(test_prop, test_prop$vireo_sample != "unassigned" & test_prop$vireo_sample !="doublet") 
        prop_vireo_correct <- subset(prop_vireo, vireo_sample == truth)
        prop_usable_error_vireo <- 1-nrow(prop_vireo_correct)/nrow(prop_vireo)
        
        ## souporcell
        prop_souporcell  <- subset(test_prop, test_prop$souporcell_sample != "unassigned" & test_prop$souporcell_sample !="doublet")         
        prop_souporcell_correct <- subset(prop_souporcell, souporcell_sample == truth)
        prop_usable_error_souporcell <- 1- nrow(prop_souporcell_correct)/nrow(prop_souporcell)
        
        ## Ensemblux
        test_prop$assignment_R_consensus[test_prop$singlet_confidence < 1.0] <- "unassigned"
        prop_ensemblux  <- subset(test_prop, test_prop$assignment_R_consensus != "unassigned" & test_prop$assignment_R_consensus !="doublet") 
        prop_ensemblux_correct <- subset(prop_ensemblux, assignment_R_consensus == truth)
        prop_usable_error_ensemblux <- 1- nrow(prop_ensemblux_correct)/nrow(prop_ensemblux)
            
        ## Dataframe
        tool <- c("Ensemblux", "Demuxalot", "Demuxlet", "Souporcell", "Vireo")
        prop_usable_error <- c(prop_usable_error_ensemblux, prop_usable_error_demuxalot, prop_usable_error_demuxlet, prop_usable_error_souporcell, prop_usable_error_vireo)
        prop_usable_error_df <- data.frame(tool, prop_usable_error)
        prop_usable_error_df$tool <- factor(prop_usable_error_df$tool, levels = c("Ensemblux", "Demuxalot", "Demuxlet", "Souporcell", "Vireo"))
        
        #######################
        # Adjusted Rand Index # 
        #######################
        test_prop <- eval_df
        
        ### Apply threshold for Ensemblux and Demuxalot and make sure best guess is doublet even if "unassigned"
        ## Demuxalot
        test_prop$demuxalot_test <- test_prop$demuxalot
        test_prop$demuxalot_test[test_prop$demuxalot_max < 0.9 & test_prop$demuxalot != "doublet"] <- "unassigned"
        
        ## Ensemblux
        test_prop$assignment_R_consensus[test_prop$singlet_confidence < 1.0 & test_prop$assignment_R_consensus != "doublet"] <- "unassigned"
        
        ## Souporcell
        test_prop$souporcell_sample[test_prop$souporcell == "doublet"] <- "doublet"
        
        ## Demuxlet
        test_prop$freemuxlet_sample[test_prop$freemuxlet == "doublet"] <- "doublet"
        
        ## Vireo
        test_prop$vireo_sample[test_prop$vireo == "doublet"] <- "doublet"
        
        ## Calculate ARI
        ARI_vireo_truth <- adj.rand.index(test_prop$vireo_sample, test_prop$truth) 
        ARI_ensemblux_truth <- adj.rand.index(test_prop$assignment_R_consensus, test_prop$truth) 
        ARI_demuxalot_truth <- adj.rand.index(test_prop$demuxalot_test, test_prop$truth) 
        ARI_demuxlet_truth <- adj.rand.index(test_prop$freemuxlet_sample, test_prop$truth) 
        ARI_souporcell_truth <- adj.rand.index(test_prop$souporcell_sample, test_prop$truth) 

        ## Dataframe
        tool <- c("Ensemblux", "Demuxalot", "Demuxlet", "Souporcell", "Vireo")
        ARI_thresh <- c(ARI_ensemblux_truth, ARI_demuxalot_truth, ARI_demuxlet_truth, ARI_souporcell_truth, ARI_vireo_truth)
        ARI_threshold_df <- data.frame(tool, ARI_thresh)
        ARI_threshold_df$tool <- factor(ARI_threshold_df$tool, levels = c("Ensemblux", "Demuxalot", "Demuxlet", "Souporcell", "Vireo"))
        

        
        #####################
        # Balanced accuracy # 
        #####################
        test_prop <- eval_df
        
        ### Apply threshold for Ensemblux and Demuxalot and make sure best guess is doublet even if "unassigned"
        ## Demuxalot
        eval_df_ba$demuxalot_test <- eval_df_ba$demuxalot
        eval_df_ba$demuxalot_test[eval_df_ba$demuxalot_max < 0.9 & eval_df_ba$demuxalot != "doublet"] <- "unassigned"
        
        ## Ensemblux
        eval_df_ba$assignment_R_consensus[eval_df_ba$singlet_confidence < 1.0 & eval_df_ba$assignment_R_consensus != "doublet"] <- "unassigned"
        
        ## Souporcell
        eval_df_ba$souporcell_sample[eval_df_ba$souporcell == "doublet"] <- "doublet"
        
        ## Demuxlet
        eval_df_ba$freemuxlet_sample[eval_df_ba$freemuxlet == "doublet"] <- "doublet"
        
        ## Vireo
        eval_df_ba$vireo_sample[eval_df_ba$vireo == "doublet"] <- "doublet"
        
        ### Compute balanced accuracy
        ## Vireo
        eval_df_ba$vireo_eval[str_detect(eval_df_ba$truth, "d") & eval_df_ba$truth == eval_df_ba$vireo_sample] <- "TN"
        eval_df_ba$vireo_eval[str_detect(eval_df_ba$truth, "s") & eval_df_ba$truth == eval_df_ba$vireo_sample] <- "TP"
        eval_df_ba$vireo_eval[str_detect(eval_df_ba$truth, "s") & eval_df_ba$truth != eval_df_ba$vireo_sample] <- "FP"
        eval_df_ba$vireo_eval[str_detect(eval_df_ba$truth, "d") & eval_df_ba$truth != eval_df_ba$vireo_sample] <- "FN"
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
        eval_df_ba$demuxalot_eval[str_detect(eval_df_ba$truth, "d") & eval_df_ba$truth == eval_df_ba$demuxalot_test] <- "TN"
        eval_df_ba$demuxalot_eval[str_detect(eval_df_ba$truth, "s") & eval_df_ba$truth == eval_df_ba$demuxalot_test] <- "TP"
        eval_df_ba$demuxalot_eval[str_detect(eval_df_ba$truth, "s") & eval_df_ba$truth != eval_df_ba$demuxalot_test] <- "FP"
        eval_df_ba$demuxalot_eval[str_detect(eval_df_ba$truth, "d") & eval_df_ba$truth != eval_df_ba$demuxalot_test] <- "FN"
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
        
        ## Ensemblux
        eval_df_ba$ensemblux_eval[str_detect(eval_df_ba$truth, "d") & eval_df_ba$truth == eval_df_ba$assignment_R_consensus] <- "TN"
        eval_df_ba$ensemblux_eval[str_detect(eval_df_ba$truth, "s") & eval_df_ba$truth == eval_df_ba$assignment_R_consensus] <- "TP"
        eval_df_ba$ensemblux_eval[str_detect(eval_df_ba$truth, "s") & eval_df_ba$truth != eval_df_ba$assignment_R_consensus] <- "FP"
        eval_df_ba$ensemblux_eval[str_detect(eval_df_ba$truth, "d") & eval_df_ba$truth != eval_df_ba$assignment_R_consensus] <- "FN"
        df <- eval_df_ba
        df_summary <- data.frame(table(df$ensemblux_eval))
        Var1 <- c("FN", "TN", "TP", "FP")
        Freq<- c(0,0,0,0)
        filler_frame <- data.frame(Var1, Freq)
        df_summary <- rbind(df_summary, filler_frame)
        df_summary <- df_summary %>% dplyr::group_by(Var1) %>% 
            dplyr::summarise(Freq = sum(Freq)) %>% as.data.frame()
        BA <- ((df_summary$Freq[df_summary$Var1 == "TP"]/(df_summary$Freq[df_summary$Var1 == "TP"] + df_summary$Freq[df_summary$Var1 == "FN"])) + 
                (df_summary$Freq[df_summary$Var1 == "TN"]/(df_summary$Freq[df_summary$Var1 == "TN"] + df_summary$Freq[df_summary$Var1 == "FP"])))/2
        ensemblux_BA <- BA
        
        ## Demuxlet
        eval_df_ba$demuxlet_eval[str_detect(eval_df_ba$truth, "d") & eval_df_ba$truth == eval_df_ba$freemuxlet_sample] <- "TN"
        eval_df_ba$demuxlet_eval[str_detect(eval_df_ba$truth, "s") & eval_df_ba$truth == eval_df_ba$freemuxlet_sample] <- "TP"
        eval_df_ba$demuxlet_eval[str_detect(eval_df_ba$truth, "s") & eval_df_ba$truth != eval_df_ba$freemuxlet_sample] <- "FP"
        eval_df_ba$demuxlet_eval[str_detect(eval_df_ba$truth, "d") & eval_df_ba$truth != eval_df_ba$freemuxlet_sample] <- "FN"
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
        eval_df_ba$souporcell_eval[str_detect(eval_df_ba$truth, "d") & eval_df_ba$truth == eval_df_ba$souporcell_sample] <- "TN"
        eval_df_ba$souporcell_eval[str_detect(eval_df_ba$truth, "s") & eval_df_ba$truth == eval_df_ba$souporcell_sample] <- "TP"
        eval_df_ba$souporcell_eval[str_detect(eval_df_ba$truth, "s") & eval_df_ba$truth != eval_df_ba$souporcell_sample] <- "FP"
        eval_df_ba$souporcell_eval[str_detect(eval_df_ba$truth, "d") & eval_df_ba$truth != eval_df_ba$souporcell_sample] <- "FN"
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
        tool <- c("Ensemblux", "Demuxalot", "Demuxlet", "Souporcell", "Vireo")
        BA_thresh <- c(ensemblux_BA, demuxalot_BA, demuxlet_BA,souporcell_BA, vireo_BA) 
        BA_threshold_df <- data.frame(tool, BA_thresh)
        BA_threshold_df$tool <- factor(BA_threshold_df$tool, levels = c("Ensemblux", "Demuxalot", "Demuxlet", "Souporcell", "Vireo"))
         
        #####################################
        # Matthew's Correlation Coefficient # 
        #####################################
        eval_df_ba <- eval_df
        
        ### Apply threshold for Ensemblux and Demuxalot and make sure best guess is doublet even if "unassigned"
        ## Demuxalot
        eval_df_ba$demuxalot_test <- eval_df_ba$demuxalot
        eval_df_ba$demuxalot_test[eval_df_ba$demuxalot_max < 0.9 & eval_df_ba$demuxalot != "doublet"] <- "unassigned"
        
        ## Ensemblux
        eval_df_ba$assignment_R_consensus[eval_df_ba$singlet_confidence < 1.0 & eval_df_ba$assignment_R_consensus != "doublet"] <- "unassigned"
        
        ## Souporcell
        eval_df_ba$souporcell_sample[eval_df_ba$souporcell == "doublet"] <- "doublet"
        
        ## Demuxlet
        eval_df_ba$freemuxlet_sample[eval_df_ba$freemuxlet == "doublet"] <- "doublet"
        
        ## Vireo
        eval_df_ba$vireo_sample[eval_df_ba$vireo == "doublet"] <- "doublet"
        
        ## Vireo
        eval_df_ba$vireo_eval[str_detect(eval_df_ba$truth, "d") & eval_df_ba$truth == eval_df_ba$vireo_sample] <- "TN"
        eval_df_ba$vireo_eval[str_detect(eval_df_ba$truth, "s") & eval_df_ba$truth == eval_df_ba$vireo_sample] <- "TP"
        eval_df_ba$vireo_eval[str_detect(eval_df_ba$truth, "s") & eval_df_ba$truth != eval_df_ba$vireo_sample] <- "FP"
        eval_df_ba$vireo_eval[str_detect(eval_df_ba$truth, "d") & eval_df_ba$truth != eval_df_ba$vireo_sample] <- "FN"
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
        eval_df_ba$demuxalot_eval[str_detect(eval_df_ba$truth, "d") & eval_df_ba$truth == eval_df_ba$demuxalot_test] <- "TN"
        eval_df_ba$demuxalot_eval[str_detect(eval_df_ba$truth, "s") & eval_df_ba$truth == eval_df_ba$demuxalot_test] <- "TP"
        eval_df_ba$demuxalot_eval[str_detect(eval_df_ba$truth, "s") & eval_df_ba$truth != eval_df_ba$demuxalot_test] <- "FP"
        eval_df_ba$demuxalot_eval[str_detect(eval_df_ba$truth, "d") & eval_df_ba$truth != eval_df_ba$demuxalot_test] <- "FN"
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
        
        ## Ensemblux
        eval_df_ba$ensemblux_eval[str_detect(eval_df_ba$truth, "d") & eval_df_ba$truth == eval_df_ba$assignment_R_consensus] <- "TN"
        eval_df_ba$ensemblux_eval[str_detect(eval_df_ba$truth, "s") & eval_df_ba$truth == eval_df_ba$assignment_R_consensus] <- "TP"
        eval_df_ba$ensemblux_eval[str_detect(eval_df_ba$truth, "s") & eval_df_ba$truth != eval_df_ba$assignment_R_consensus] <- "FP"
        eval_df_ba$ensemblux_eval[str_detect(eval_df_ba$truth, "d") & eval_df_ba$truth != eval_df_ba$assignment_R_consensus] <- "FN"
        df <- eval_df_ba
        df_summary <- data.frame(table(df$ensemblux_eval))
        Var1 <- c("FN", "TN", "TP", "FP")
        Freq<- c(0,0,0,0)
        filler_frame <- data.frame(Var1, Freq)
        df_summary <- rbind(df_summary, filler_frame)
        df_summary <- df_summary %>% dplyr::group_by(Var1) %>% 
            dplyr::summarise(Freq = sum(Freq)) %>% as.data.frame()
        MCC <- ((as.numeric((df_summary$Freq[df_summary$Var1 == "TP"])*as.numeric(df_summary$Freq[df_summary$Var1 == "TN"]))-(as.numeric(df_summary$Freq[df_summary$Var1 == "FP"])*as.numeric(df_summary$Freq[df_summary$Var1 == "FN"])))/
                    sqrt((as.numeric(df_summary$Freq[df_summary$Var1 == "TP"])+as.numeric(df_summary$Freq[df_summary$Var1 == "FP"]))*(as.numeric(df_summary$Freq[df_summary$Var1 == "TP"])+as.numeric(df_summary$Freq[df_summary$Var1 == "FN"]))*(as.numeric(df_summary$Freq[df_summary$Var1 == "TN"])+as.numeric(df_summary$Freq[df_summary$Var1 == "FP"]))*(as.numeric(df_summary$Freq[df_summary$Var1 == "TN"])+as.numeric(df_summary$Freq[df_summary$Var1 == "FN"]))))
        ensemblux_MCC <- MCC
        
        ## Demuxlet
        eval_df_ba$demuxlet_eval[str_detect(eval_df_ba$truth, "d") & eval_df_ba$truth == eval_df_ba$freemuxlet_sample] <- "TN"
        eval_df_ba$demuxlet_eval[str_detect(eval_df_ba$truth, "s") & eval_df_ba$truth == eval_df_ba$freemuxlet_sample] <- "TP"
        eval_df_ba$demuxlet_eval[str_detect(eval_df_ba$truth, "s") & eval_df_ba$truth != eval_df_ba$freemuxlet_sample] <- "FP"
        eval_df_ba$demuxlet_eval[str_detect(eval_df_ba$truth, "d") & eval_df_ba$truth != eval_df_ba$freemuxlet_sample] <- "FN"
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
        
        ## Souporcell
        eval_df_ba$souporcell_eval[str_detect(eval_df_ba$truth, "d") & eval_df_ba$truth == eval_df_ba$souporcell_sample] <- "TN"
        eval_df_ba$souporcell_eval[str_detect(eval_df_ba$truth, "s") & eval_df_ba$truth == eval_df_ba$souporcell_sample] <- "TP"
        eval_df_ba$souporcell_eval[str_detect(eval_df_ba$truth, "s") & eval_df_ba$truth != eval_df_ba$souporcell_sample] <- "FP"
        eval_df_ba$souporcell_eval[str_detect(eval_df_ba$truth, "d") & eval_df_ba$truth != eval_df_ba$souporcell_sample] <- "FN"
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
        tool <- c("Ensemblux", "Demuxalot", "Demuxlet", "Souporcell", "Vireo")
        MCC_thresh <- c(ensemblux_MCC, demuxalot_MCC, demuxlet_MCC, souporcell_MCC, vireo_MCC)
        MCC_threshold_df <- data.frame(tool, MCC_thresh)
        MCC_threshold_df$tool <- factor(MCC_threshold_df$tool, levels = c("Ensemblux", "Demuxalot", "Demuxlet", "Souporcell", "Vireo"))

        
        #########################
        # AUC singlet detection # 
        #########################
        ## load in the dataset
        eval_df_ba <- eval_df
        
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
        
        ## Ensemblux 
        eval_df_ba$concensus_eval <- "bad"
        eval_df_ba$concensus_eval[eval_df_ba$assignment_R_consensus == eval_df_ba$truth ] <- "good"
            
        ## subset main dataframe to only include true singlets
        df <- subset(eval_df_ba, truth != "doublet")
    
        ### calculate AUC
        ## Demuxalot
        roc_empirical <- rocit(score = df$demuxalot_max, class = df$demuxalot_eval, ##change_here
                                negref = "bad") 
        summary(roc_empirical)
        AUC_demuxalot <- roc_empirical$AUC
        
        ## Demuxlet
        roc_empirical <- rocit(score = df$freemuxlet_max_posterior, class = df$demuxlet_eval, ##change_here
                                negref = "bad") 
        summary(roc_empirical)
        AUC_demuxlet <- roc_empirical$AUC
        AUC_demuxlet
        
        ## Vireo
        roc_empirical <- rocit(score = df$vireo_max, class = df$vireo_eval, ##change_here
                                negref = "bad") 
        summary(roc_empirical)
        AUC_vireo <-roc_empirical$AUC
        AUC_vireo
        
        ## Souporcell
        roc_empirical <- rocit(score =  1-(10^(df$souporcell_log_prob_singleton)), class = df$souporcell_eval, ##change_here
                                negref = "bad") 
        summary(roc_empirical)
        AUC_souporcell <-roc_empirical$AUC
        AUC_souporcell
        
        ## Ensemblux
        roc_empirical <- rocit(score = df$singlet_confidence, class = df$concensus_eval, ##change_here
                                negref = "bad") 
        summary(roc_empirical)
        AUC_ensemblux <-roc_empirical$AUC
        AUC_ensemblux
            
        ## plot
        tool <- c("Ensemblux", "Demuxalot", "Demuxlet", "Souporcell", "Vireo")
        AUC_singlet <- c(AUC_ensemblux, AUC_demuxalot, AUC_demuxlet, AUC_souporcell, AUC_vireo)
        AUC_singlet_df <- data.frame(tool, AUC_singlet)
        AUC_singlet_df$tool <- factor(AUC_singlet_df$tool, levels = c("Ensemblux", "Demuxalot", "Demuxlet", "Souporcell", "Vireo"))
}


