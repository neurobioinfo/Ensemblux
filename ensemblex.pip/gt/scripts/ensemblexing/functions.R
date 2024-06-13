###########################################################################################################################
# CONSTITUENT OUTPUT PREP
########################################################################################################################### 
demuxalot_prep <- function(demuxalot){
    ## Compute number of columns in the dataframe
    xx <- ncol(demuxalot)

    ## Change doublet column names for downstream operations
    list_col <- c((par_sample_size + 2):xx)
    for (i in list_col){
    colnames(demuxalot)[i] <- paste0("doublet_",i)
    }
    
    ## Place Demuxalot output in long format
    data_long <- gather(demuxalot, assignment, stat, 2:xx, factor_key=F)
    
    ## Identify the top sample prediction for each barcode
    result <- data_long %>% 
    dplyr::group_by(X) %>%
    dplyr::filter(stat == max(stat))
    nrow(result)
    
    ## Identify uniadentifiable cells; these are barcodes that have two or more samples assigned with the same probability
    duplicates <- subset(result,duplicated(X))
    duplicated_id <- duplicates$X
    result$assignment[result$X %in% duplicated_id] <- "unassigned"

    ## Remove duplicated barcodes to only keep one sample assignment; these barcodes are labelled as unassigned. 
    result <- subset(result,!duplicated(X))
    
    ## Set up demuxalot assignment
    result$demuxalot_sample <- result$assignment
    result$demuxalot_sample[str_detect(result$demuxalot_sample, "doublet")] <- "doublet"

    ## Identify Demuxalot doublet probability
    # Keep only doublet assignments
    data_long_2 <- data_long[str_detect(data_long$assignment, "doublet"),]
    # Keep only the top prediction for each barcode
    result_doublet <- data_long_2 %>% 
    dplyr::group_by(X) %>%
    dplyr::filter(stat == max(stat))
    # Keep only one of each barcode --> some barcodes will have multiple doublet combinations with the same doublet probability
    result_doublet <- result_doublet[!(duplicated(result_doublet$X)),]
    result_doublet <- result_doublet[,c(1,3)]
    colnames(result_doublet) <- c("X", "demuxalot_doublet_probability")
    
    ## Merge doublet probability back to main dataframe
    nrow(result) == nrow(result_doublet)
    demuxalot_2 <- merge(result, result_doublet, by = "X")
    nrow(demuxalot_2) == nrow(demuxalot)
    demuxalot <- demuxalot_2
    
    ## Print number of retained droplets
    message(paste0("Retained ", nrow(demuxalot), " droplets during Demuxalot preparation."))
    demuxalot
}

vireo_prep <- function(vireo){
    ## Vireo's sample assignment
    vireo$vireo_sample <- vireo$donor_id

    ## Prepare Vireo's best guess assignment
    vireo$vireo_best <- vireo$donor_id
    vireo$vireo_best[vireo$vireo_best == "unassigned" & vireo$prob_max > vireo$prob_doublet] <- vireo$best_singlet[vireo$vireo_best == "unassigned" & vireo$prob_max > vireo$prob_doublet]
    vireo$vireo_best[vireo$vireo_best == "unassigned" & vireo$prob_max < vireo$prob_doublet] <- "doublet"

    ## Print number of retained droplets
    message(paste0("Retained ", nrow(vireo), " droplets during Vireo preparation."))
    vireo
}

freemuxlet_prep <- function(freemuxlet){
    ## Prepare Demuxlet best guess
    freemuxlet$guess_1 <- sub(".*,(.*),.*", "\\1", freemuxlet$BEST.GUESS)
    freemuxlet$guess_2 <-sub("(.*),.*,.*", "\\1", freemuxlet$BEST.GUESS)
    freemuxlet$freemuxlet_best[freemuxlet$guess_1 != freemuxlet$guess_2] <- "doublet"  
    freemuxlet$freemuxlet_best[freemuxlet$guess_1 == freemuxlet$guess_2] <- freemuxlet$guess_1[freemuxlet$guess_1 == freemuxlet$guess_2] 
    
    ## Prepare Demuxlet sample assignment
    freemuxlet$freemuxlet_sample <- freemuxlet$freemuxlet_best
    freemuxlet$freemuxlet_sample[freemuxlet$DROPLET.TYPE == "AMB"] <- "unassigned"
    freemuxlet <- freemuxlet %>% dplyr::select(-c(guess_1, guess_2))
    
    ## Print number of retained droplets
    message(paste0("Retained ", nrow(freemuxlet), " droplets during Demuxlet preparation."))
    freemuxlet
}

souporcell_prep <- function(souporcell, demuxalot, vireo, freemuxlet){
    
    ## Number of singlet cluster identified by Souporcell
    soup_singlet_cluster <- unique(souporcell$assignment)
    soup_singlet_cluster <- soup_singlet_cluster[!grepl("/", soup_singlet_cluster)]
    
    ## Demuxalot
    demuxalot_for_soup <- demuxalot[,c(1,4)]
    nrow(demuxalot_for_soup) == nrow(souporcell)
    ## Merge Souporcell and Demuxalot
    merge_soup_demuxalot <- merge(demuxalot_for_soup, souporcell, by.x = "X", by.y = "barcode")
    nrow(merge_soup_demuxalot) == nrow(demuxalot_for_soup)
    merge_soup_demuxalot <- merge_soup_demuxalot[,c(2:4)]
    colnames(merge_soup_demuxalot) <- c("Sample_ID", "status", "assignment")

    ## Vireo
    vireo_for_soup <- vireo[,c(1,10)]
    nrow(vireo_for_soup) == nrow(souporcell)
    ## Merge Souporcell and Vireo
    merge_soup_vireo <- merge(vireo_for_soup, souporcell, by.x = "cell", by.y = "barcode")
    nrow(merge_soup_vireo) == nrow(vireo_for_soup)
    merge_soup_vireo <- merge_soup_vireo[,c(2:4)]
    colnames(merge_soup_vireo) <- c("Sample_ID", "status", "assignment")

    ## Demuxlet
    demuxlet_for_soup <- freemuxlet[,c(2,21)]
    nrow(demuxlet_for_soup) == nrow(souporcell)
    ## Merge Souporcell and Demuxlet
    merge_soup_demuxlet <- merge(demuxlet_for_soup, souporcell, by.x = "BARCODE", by.y = "barcode")
    nrow(merge_soup_demuxlet) == nrow(demuxlet_for_soup)
    merge_soup_demuxlet <- merge_soup_demuxlet[,c(2:4)]
    colnames(merge_soup_demuxlet) <- c("Sample_ID", "status", "assignment")

    ### Calculate cluster-Sample ID probabilities based on assignments of remaining constituent tools
    bind_df <- rbind(merge_soup_demuxalot, merge_soup_vireo, merge_soup_demuxlet)
    bind_df <- subset(bind_df, status != "doublet")
    bind_df <- subset(bind_df, status != "unassigned")
    bind_df <- subset(bind_df, Sample_ID != "doublet")
    souporcell_assignment = "none"
    Sample_ID = "none"
    Probability = "none"
    data_frame_final <- data.frame(souporcell_assignment, Sample_ID, Probability)

    for(i in unique(bind_df$assignment)){
        df_lim <- subset(bind_df, assignment == i)
        n_droplets <- nrow(df_lim)

        souporcell_assignment = "none"
        Sample_ID = "none"
        Probability = "none"
        data_frame <- data.frame(souporcell_assignment, Sample_ID, Probability)

        for (j in unique(bind_df$Sample_ID)){
        df_lim2 <- subset(df_lim, Sample_ID == j )
        prob <- nrow(df_lim2)/n_droplets
        
        souporcell_assignment = i
        Sample_ID = j
        Probability = prob
        data_frame_temp <- data.frame(souporcell_assignment, Sample_ID, Probability)
        data_frame <- rbind(data_frame, data_frame_temp)
        }
        data_frame_final <- rbind(data_frame_final, data_frame)
    }

    data_frame_final <- subset(data_frame_final, Probability != "none")
    data_frame_final$Probability <- as.numeric(data_frame_final$Probability)

    ## take maximum cluster-sample probability
    result <- data_frame_final %>% 
    dplyr::group_by(souporcell_assignment) %>%
    dplyr::filter(Probability == max(Probability))

    ## set duplicated clusters to unassigned
    n_occur <- data.frame(table(result$souporcell_assignment))
    n_occur <- subset(n_occur, Freq > 1)
    duplicated <- n_occur$Var1
    result$Sample_ID[result$souporcell_assignment %in% duplicated] <- "unassigned"

    #temp_df <- data.frame(result)
    result <- result[!duplicated(result$souporcell_assignment),]


    if(length(result$souporcell_assignment) == length(unique(bind_df$assignment)) ){
        result <- data.frame(result)
        result <- result[,c(1:2)]
        colnames(result) <- c("assignment","souporcell_sample")
        souporcell <- merge(souporcell, result, by = "assignment", all = T)
    }


    souporcell$souporcell_sample[souporcell$status == 'doublet'] <- 'doublet'
    souporcell$souporcell_sample[str_detect(souporcell$assignment, "/")] <- "doublet" 

    souporcell$souporcell_sample[is.na(souporcell$souporcell_sample)] <- "unassigned"
    unique(souporcell$souporcell_sample)
    
    souporcell$souporcell_best <- souporcell$souporcell_sample
    souporcell$souporcell_sample[souporcell$status == 'unassigned'] <- 'unassigned'

    if(length(unique(result$souporcell_sample)) == par_sample_size){
        message("Successfully matched all Souporcell clusters to Sample ID.")
    } else {
        message(paste0("WARNING: ensemblex failed to match all Souporcell clusters to all Sample ID. ensemblex will proceed; however, it is recommended to look into the data manually."))
    }

    ## Print number of retained droplets
    message(paste0("Retained ", nrow(souporcell), " droplets during Demuxlet preparation."))
    souporcell
}

merge_concensus <- function(vireo,souporcell, freemuxlet, demuxalot){
  ## Select important columns from Vireo
  vireo <- dplyr::select(vireo, c("cell", "prob_max", "prob_doublet", "n_vars", "best_doublet", "doublet_logLikRatio", "vireo_sample", "vireo_best" ))
  colnames(vireo) <- c("barcode",  "vireo_singlet_probability", "vireo_doublet_probability", "vireo_n_vars", "vireo_best_doublet", "vireo_doublet_logLikRatio", "vireo_assignment", "vireo_best_assignment")
  
  ## Select important columns from Souporcell
  souporcell <- dplyr::select(souporcell, c("barcode", "log_prob_singleton", "log_prob_doublet",  "souporcell_sample", "souporcell_best"))
  colnames(souporcell) <- c("barcode", "souporcell_log_probability_singlet", "souporcell_log_probability_doublet",  "souporcell_assignment", "souporcell_best_assignment")
  
  ## Select important columns from Demuxlet
  freemuxlet <- dplyr::select(freemuxlet, c("BARCODE", "NUM.SNPS", "NUM.READS", "SNG.POSTERIOR", "freemuxlet_sample", "freemuxlet_best", "DIFF.LLK.SNG.DBL"))
  colnames(freemuxlet) <- c("barcode", "demuxlet_n_snps",  "demuxlet_n_reads", "demuxlet_max_probability",  "demuxlet_assignment", "demuxlet_best_assignment", "demuxlet_DIFF_LLK_SNG_DBL")
  
  ## Select important columns from Demuxalot
  demuxalot$demuxalot_second <- demuxalot$demuxalot_sample
  demuxalot$demuxalot_sample[demuxalot$stat < 0.9] <- "unassigned"
  demuxalot <- dplyr::select(demuxalot, c("X", "stat", "demuxalot_sample", "demuxalot_second", "demuxalot_doublet_probability")) 
  colnames(demuxalot) <- c("barcode", "demuxalot_max_probability", "demuxalot_assignment", "demuxalot_best_assignment", "demuxalot_doublet_probability" ) 
  
  ## Merge dataframes
  merge_df <- merge(vireo, souporcell, by = c("barcode"), all = T)
  merge_df <- merge(merge_df, freemuxlet, by = c("barcode"), all = T)
  merge_df <- merge(merge_df, demuxalot, by = c("barcode"), all = T)
    
  ## Generate a general consensus column
  merge_lim <- dplyr::select(merge_df, c("barcode", "vireo_assignment", "souporcell_assignment", "demuxlet_assignment", "demuxalot_assignment"))
  merge_lim$general_consensus <- apply(merge_lim,1,function(x) names(which.max(table(x))))
  merge_lim$general_consensus <- sub(".*(-).*", "\\1", merge_lim$general_consensus)
  merge_lim$general_consensus[merge_lim$general_consensus == "-" ] <- "unassigned"
  merge_lim[is.na(merge_lim)] <- "unassigned"
  merge_bind <-  merge_lim
  
  ## Merge back to initial dataframe
  merge_df <- merge(merge_bind, merge_df, by = c("barcode"), all = T)

  ## Print number of droplets in merged dataframe
  message(paste0("Retained ", nrow(merge_df), " after merging output files from each constituent demultiplexing tool."))

  ## Clean up output dataframe
  merge_df <- dplyr::select(merge_df, -c("vireo_assignment.y", "souporcell_assignment.y","demuxlet_assignment.y", "demuxalot_assignment.y"))
  colnames(merge_df) <- c("barcode", "vireo_assignment", "souporcell_assignment", "demuxlet_assignment", "demuxalot_assignment", "general_consensus","vireo_singlet_probability", "vireo_doublet_probability",         
                            "vireo_n_vars", "vireo_best_doublet", "vireo_doublet_logLikRatio", "vireo_best_assignment", "souporcell_log_probability_singlet", "souporcell_log_probability_doublet", "souporcell_best_assignment", 
                            "demuxlet_n_snps", "demuxlet_n_reads", "demuxlet_max_probability", "demuxlet_best_assignment", "demuxlet_DIFF_LLK_SNG_DBL", "demuxalot_max_probability", "demuxalot_best_assignment",         
                            "demuxalot_doublet_probability")
  merge_df <- dplyr::select(merge_df, c("barcode", "vireo_assignment", "souporcell_assignment", "demuxlet_assignment", "demuxalot_assignment", "general_consensus", "vireo_best_assignment", "souporcell_best_assignment",
                                        "demuxlet_best_assignment", "demuxalot_best_assignment", "vireo_singlet_probability", "vireo_doublet_probability", "vireo_n_vars", "vireo_best_doublet", "vireo_doublet_logLikRatio",  
                                        "souporcell_log_probability_singlet", "souporcell_log_probability_doublet", "demuxlet_n_snps", "demuxlet_n_reads", "demuxlet_max_probability", "demuxlet_DIFF_LLK_SNG_DBL", 
                                        "demuxalot_max_probability", "demuxalot_doublet_probability"))
  merge_df
}


###########################################################################################################################
# PROBABILISTIC-WEIGHTED ENSEMBLE 
###########################################################################################################################
## FUNCTIONS
BA_weight_consensus <- function(merge_df,par_sample_size,par_output_dir){
  
  ## Set seed 
  set.seed(1234)
  
  ## Create an output directory of probabilistic-weighted ensmeble outputs
  dir.create(paste(par_output_dir,"/step1",sep=''))
  
  ## Rename the dataset
  eval_df <- merge_df

  #### Adjusted Rand Index between sample assignments -- here we are using the best guess from each tool 
  ## Vireo
  ARI_vireo_vireo <- adj.rand.index(eval_df$vireo_best_assignment, eval_df$vireo_best_assignment) 
  ARI_vireo_demuxlet <- adj.rand.index(eval_df$vireo_best_assignment, eval_df$demuxlet_best_assignment) 
  ARI_vireo_demuxalot <- adj.rand.index(eval_df$vireo_best_assignment, eval_df$demuxalot_best_assignment) 
  ARI_vireo_souporcell <- adj.rand.index(eval_df$vireo_best_assignment, eval_df$souporcell_best_assignment) 
  ## Demuxlet
  ARI_demuxlet_demuxlet <- adj.rand.index(eval_df$demuxlet_best_assignment, eval_df$demuxlet_best_assignment) 
  ARI_demuxlet_vireo <- adj.rand.index(eval_df$demuxlet_best_assignment, eval_df$vireo_best_assignment) 
  ARI_demuxlet_demuxalot <- adj.rand.index(eval_df$demuxlet_best_assignment, eval_df$demuxalot_best_assignment) 
  ARI_demuxlet_souporcell <- adj.rand.index(eval_df$demuxlet_best_assignment, eval_df$souporcell_best_assignment) 
  ## Demuxalot
  ARI_demuxalot_demuxalot <- adj.rand.index(eval_df$demuxalot_best_assignment, eval_df$demuxalot_best_assignment) 
  ARI_demuxalot_souporcell <- adj.rand.index(eval_df$demuxalot_best_assignment, eval_df$souporcell_best_assignment) 
  ARI_demuxalot_demuxlet <- adj.rand.index(eval_df$demuxalot_best_assignment, eval_df$demuxlet_best_assignment) 
  ARI_demuxalot_vireo <- adj.rand.index(eval_df$demuxalot_best_assignment, eval_df$vireo_best_assignment) 
  ## Souporcell
  ARI_souporcell_souporcell <- adj.rand.index(eval_df$souporcell_best_assignment, eval_df$souporcell_best_assignment) 
  ARI_souporcell_vireo <- adj.rand.index(eval_df$souporcell_best_assignment, eval_df$vireo_best_assignment) 
  ARI_souporcell_demuxlet <- adj.rand.index(eval_df$souporcell_best_assignment, eval_df$demuxlet_best_assignment) 
  ARI_souporcell_demuxalot <- adj.rand.index(eval_df$souporcell_best_assignment, eval_df$demuxalot_best_assignment) 
  
  ## Produce data frame
  tool_1 <- c("Vireo","Vireo","Vireo","Vireo",
              "Demuxlet", "Demuxlet", "Demuxlet", "Demuxlet",
              "Demuxalot","Demuxalot","Demuxalot","Demuxalot",
              "Souporcell","Souporcell","Souporcell", "Souporcell")
  tool_2 <- c("Vireo", "Demuxlet", "Demuxalot", "Souporcell",
              "Demuxlet", "Vireo","Demuxalot","Souporcell",
              "Demuxalot", "Souporcell","Demuxlet","Vireo",
              "Souporcell","Vireo","Demuxlet","Demuxalot")
  ARI <- c(ARI_vireo_vireo, ARI_vireo_demuxlet, ARI_vireo_demuxalot, ARI_vireo_souporcell,
           ARI_demuxlet_demuxlet, ARI_demuxlet_vireo, ARI_demuxlet_demuxalot, ARI_demuxlet_souporcell,
           ARI_demuxalot_demuxalot, ARI_demuxalot_souporcell, ARI_demuxalot_demuxlet, ARI_demuxalot_vireo,
           ARI_souporcell_souporcell, ARI_souporcell_vireo, ARI_souporcell_demuxlet, ARI_souporcell_demuxalot)
  ARI_df <- data.frame(tool_1,tool_2,ARI )
  
  ## Plot ARI heatmap
  ggplot(ARI_df, aes(x = tool_1, y = tool_2, fill = ARI, label = round(ARI, digits = 3) )) +geom_tile() + theme_bw() +
    scale_fill_gradient(low="white", high="darkblue") +
    xlab("Demultiplexing tool") +
    ylab("Demultiplexing tool") +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0)) + geom_text()
  ggsave(paste(par_output_dir,"/step1","/ARI_demultiplexing_tools.pdf", sep=""))
  
  #### Compute proxy balanced accuracies for each tool using consensus cells for the remaining tools
  ### Vireo 
  ## Compute consensus cells
  eval_df_ba <- eval_df[eval_df$souporcell_best_assignment == eval_df$demuxlet_best_assignment &
                          eval_df$souporcell_best_assignment == eval_df$demuxalot_best_assignment &
                          eval_df$souporcell_best_assignment != "unassigned" ,] 
  vireo_n <- nrow(eval_df_ba)
  
  ## Compute balanced accuracy
  eval_df_ba$vireo_eval[eval_df_ba$souporcell_best_assignment == "doublet" & eval_df_ba$souporcell_best_assignment == eval_df_ba$vireo_best_assignment] <- "TN"
  eval_df_ba$vireo_eval[eval_df_ba$souporcell_best_assignment != "doublet" & eval_df_ba$souporcell_best_assignment == eval_df_ba$vireo_best_assignment] <- "TP"
  eval_df_ba$vireo_eval[eval_df_ba$souporcell_best_assignment != "doublet" & eval_df_ba$souporcell_best_assignment != eval_df_ba$vireo_best_assignment] <- "FP"
  eval_df_ba$vireo_eval[eval_df_ba$souporcell_best_assignment == "doublet" & eval_df_ba$souporcell_best_assignment != eval_df_ba$vireo_best_assignment] <- "FN"
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
  vireo_df <- data.frame(BA)
  vireo_df$tool <- "Vireo"
  
  ### Demuxlet
  ## Compute consensus cells
  eval_df_ba <- eval_df[eval_df$souporcell_best_assignment == eval_df$vireo_best_assignment &
                          eval_df$souporcell_best_assignment == eval_df$demuxalot_best_assignment &
                          eval_df$souporcell_best_assignment != "unassigned" ,] 
  demuxlet_n <- nrow(eval_df_ba)
  
  ## Compute balanced accuracy
  eval_df_ba$freemuxlet_eval[eval_df_ba$souporcell_best_assignment == "doublet" & eval_df_ba$souporcell_best_assignment == eval_df_ba$demuxlet_best_assignment] <- "TN"
  eval_df_ba$freemuxlet_eval[eval_df_ba$souporcell_best_assignment != "doublet"  & eval_df_ba$souporcell_best_assignment == eval_df_ba$demuxlet_best_assignment] <- "TP"
  eval_df_ba$freemuxlet_eval[eval_df_ba$souporcell_best_assignment != "doublet" & eval_df_ba$souporcell_best_assignment != eval_df_ba$demuxlet_best_assignment] <- "FP"
  eval_df_ba$freemuxlet_eval[eval_df_ba$souporcell_best_assignment == "doublet" & eval_df_ba$souporcell_best_assignment != eval_df_ba$demuxlet_best_assignment] <- "FN"
  df <- eval_df_ba
  df_summary <- data.frame(table(df$freemuxlet_eval))
  df_summary
  Var1 <- c("FN", "TN", "TP", "FP")
  Freq<- c(0,0,0,0)
  filler_frame <- data.frame(Var1, Freq)
  df_summary <- rbind(df_summary, filler_frame)
  df_summary <- df_summary %>% dplyr::group_by(Var1) %>% 
    dplyr::summarise(Freq = sum(Freq)) %>% as.data.frame()
  
  BA <- ((df_summary$Freq[df_summary$Var1 == "TP"]/(df_summary$Freq[df_summary$Var1 == "TP"] + df_summary$Freq[df_summary$Var1 == "FN"])) + 
           (df_summary$Freq[df_summary$Var1 == "TN"]/(df_summary$Freq[df_summary$Var1 == "TN"] + df_summary$Freq[df_summary$Var1 == "FP"])))/2
  
  demuxlet_BA <- BA
  freemuxlet_df <- data.frame(BA)
  freemuxlet_df$tool <- "Demuxlet"
  
  ### Demuxalot
  ## Compute consensus cells
  eval_df_ba <- eval_df[eval_df$souporcell_best_assignment == eval_df$vireo_best_assignment &
                          eval_df$souporcell_best_assignment == eval_df$demuxlet_best_assignment &
                          eval_df$souporcell_best_assignment != "unassigned" ,] 
  demuxalot_n <- nrow(eval_df_ba)
  
  ## Compute balanced accuracy
  eval_df_ba$demuxalot_eval[eval_df_ba$souporcell_best_assignment == "doublet" & eval_df_ba$souporcell_best_assignment == eval_df_ba$demuxalot_best_assignment] <- "TN"
  eval_df_ba$demuxalot_eval[eval_df_ba$souporcell_best_assignment != "doublet" & eval_df_ba$souporcell_best_assignment == eval_df_ba$demuxalot_best_assignment] <- "TP"
  eval_df_ba$demuxalot_eval[eval_df_ba$souporcell_best_assignment != "doublet" & eval_df_ba$souporcell_best_assignment != eval_df_ba$demuxalot_best_assignment] <- "FP"
  eval_df_ba$demuxalot_eval[eval_df_ba$souporcell_best_assignment == "doublet"  & eval_df_ba$souporcell_best_assignment != eval_df_ba$demuxalot_best_assignment] <- "FN"
  df <- eval_df_ba
  df_summary <- data.frame(table(df$demuxalot_eval))
  df_summary
  Var1 <- c("FN", "TN", "TP", "FP")
  Freq<- c(0,0,0,0)
  filler_frame <- data.frame(Var1, Freq)
  df_summary <- rbind(df_summary, filler_frame)
  df_summary <- df_summary %>% dplyr::group_by(Var1) %>% 
    dplyr::summarise(Freq = sum(Freq)) %>% as.data.frame()
  
  BA <- ((df_summary$Freq[df_summary$Var1 == "TP"]/(df_summary$Freq[df_summary$Var1 == "TP"] + df_summary$Freq[df_summary$Var1 == "FN"])) + 
           (df_summary$Freq[df_summary$Var1 == "TN"]/(df_summary$Freq[df_summary$Var1 == "TN"] + df_summary$Freq[df_summary$Var1 == "FP"])))/2
  
  demuxalot_BA <- BA
  demuxalot_df <- data.frame(BA)
  demuxalot_df$tool <- "Demuxalot"
  
  ### Souporcell 
  ## Compute consensus cells
  eval_df_ba <- eval_df[eval_df$demuxalot_best_assignment == eval_df$vireo_best_assignment &
                          eval_df$demuxalot_best_assignment == eval_df$demuxlet_best_assignment &
                          eval_df$demuxalot_best_assignment != "unassigned" ,] 
  souporcell_n <- nrow(eval_df_ba)
  
  ## Compute balanced accuracy
  eval_df_ba$souporcell_eval[eval_df_ba$demuxalot_best_assignment == "doublet" & eval_df_ba$demuxalot_best_assignment == eval_df_ba$souporcell_best_assignment] <- "TN"
  eval_df_ba$souporcell_eval[eval_df_ba$demuxalot_best_assignment != "doublet" & eval_df_ba$demuxalot_best_assignment == eval_df_ba$souporcell_best_assignment] <- "TP"
  eval_df_ba$souporcell_eval[eval_df_ba$demuxalot_best_assignment != "doublet" & eval_df_ba$demuxalot_best_assignment != eval_df_ba$souporcell_best_assignment] <- "FP"
  eval_df_ba$souporcell_eval[eval_df_ba$demuxalot_best_assignment == "doublet" & eval_df_ba$demuxalot_best_assignment != eval_df_ba$souporcell_best_assignment] <- "FN"
  df <- eval_df_ba
  df_summary <- data.frame(table(df$souporcell_eval))
  df_summary
  Var1 <- c("FN", "TN", "TP", "FP")
  Freq<- c(0,0,0,0)
  filler_frame <- data.frame(Var1, Freq)
  df_summary <- rbind(df_summary, filler_frame)
  df_summary <- df_summary %>% dplyr::group_by(Var1) %>% 
    dplyr::summarise(Freq = sum(Freq)) %>% as.data.frame()
  
  BA <- ((df_summary$Freq[df_summary$Var1 == "TP"]/(df_summary$Freq[df_summary$Var1 == "TP"] + df_summary$Freq[df_summary$Var1 == "FN"])) + 
           (df_summary$Freq[df_summary$Var1 == "TN"]/(df_summary$Freq[df_summary$Var1 == "TN"] + df_summary$Freq[df_summary$Var1 == "FP"])))/2
  
  souporcell_BA <- BA
  souporcell_df <- data.frame(BA)
  souporcell_df$tool <- "Souporcell"
  
  #### Output summary information
  ## Produce table with PWE information
  Tool <- c("Vireo", "Demuxlet", "Demuxalot", "Souporcell")
  Balanced_accuracy <- c(vireo_BA, demuxlet_BA, demuxalot_BA, souporcell_BA)
  n_consensus_droplets <- c(vireo_n, demuxlet_n, demuxalot_n, souporcell_n)
  PWE_summary_df <- data.frame(Tool,Balanced_accuracy,n_consensus_droplets )
  write.csv(PWE_summary_df, paste(par_output_dir,"/step1",'/Balanced_accuracy_summary.csv', sep=""))
  
  ## Plot estimated balanced accuracy
  ggplot(PWE_summary_df, aes(x = Tool, y = Balanced_accuracy, label = round(Balanced_accuracy, digits = 4), fill = Tool)) + 
    geom_bar(stat = "identity") +theme_classic() +
    geom_text() +
    xlab("Demultiplexing tool") +
    ylab("Estimated Balanced Accuracy") + 
    scale_fill_manual(values = c("#d95f02", "#e6ab02", "#7570b3", "#66a61e")) +
    theme(legend.position = "right")
  ggsave(paste(par_output_dir,"/step1","/BA_demultiplexing_tools.pdf", sep=""))
  
  ### combine balanced accuracies into one data frame ####
  BA_df <- rbind(vireo_df,
                 freemuxlet_df,
                 demuxalot_df,
                 souporcell_df)
  
  ### Place Vireo and Souporcell probabilties into proper format for downstream analsyes
  ## Vireo
  # Merge Vireo singlet and doublet probabilities -- take max 
  # Columns 8 and 9 are vireo_max and vireo_doublet, respectively
  eval_df[, "vireo_max_probability"] <- apply(eval_df[, 12:13], 1, max)  
  
  ## Souporcell -- convert log(prob) to prob
  # Singlets
  eval_df$souporcell_singlet_probability <- 1-(10^(eval_df$souporcell_log_probability_singlet))
  # Doublets
  eval_df$souporcell_doublet_probability <- 1-(10^(eval_df$souporcell_log_probability_doublet))
  # Take souporcell max probability between singlets and doublets
  # Columns 25 and 26 are souporcell_singlet and souporcell_doublet, respectively
  eval_df[, "souporcell_max_probability"] <- apply(eval_df[, 26:27], 1, max) 
  
  ### Multiply the assignment probabilities from each of the constituent demultiplexing tools by their respective estimated balanced accuracy for the dataset
  ## Vireo
  eval_df$vireo_weighted_probability <- eval_df$vireo_max_probability*BA_df$BA[BA_df$tool == "Vireo"]
  ## Demuxlet
  eval_df$demuxlet_weighted_probability <- eval_df$demuxlet_max_probability*BA_df$BA[BA_df$tool == "Demuxlet"]
  ## Demuxalot
  eval_df$demuxalot_weighted_probability <- eval_df$demuxalot_max_probability*BA_df$BA[BA_df$tool == "Demuxalot"]
  ## Souporcell
  eval_df$souporcell_weighted_probability <- eval_df$souporcell_max_probability*BA_df$BA[BA_df$tool == "Souporcell"]
  
  ############################################################################################################################################################
  ## Get weighted consensus assignment ##
  ############################################################################################################################################################
  ## Rename dataframe and remove first column
  practice_df <- eval_df
  practice_df <- practice_df[,-1]
  
  ### Create a sample list, not including unassigned and doublet
  ## Vireo
  Vireo_sample_list <- unique(practice_df$vireo_best_assignment)
  ## Demuxalot
  Demuxalot_sample_list <- unique(practice_df$demuxalot_best_assignment)
  ## Demuxlet
  Demuxlet_sample_list <- unique(practice_df$demuxlet_best_assignment)
  ## Souporcell
  Souporcell_sample_list <- unique(practice_df$souporcell_best_assignment)
  
  ## Identify all unique samples identified by each demultiplexing tool
  sample_list <- unlist(append(Vireo_sample_list,Demuxalot_sample_list))
  sample_list <- unlist(append(sample_list,Demuxlet_sample_list)) 
  sample_list <- unlist(append(sample_list,Souporcell_sample_list))   
  remove_sample <- c("doublet", "unassigned")
  sample_list_2 <- sample_list[!sample_list %in% remove_sample]
  sample_list_2 <- unique(sample_list_2)
  if (length(sample_list_2) == par_sample_size){
    message(paste0("Generating weighted-probabilistic ensemble assignments from ", length(sample_list_2), " Sample IDs."))
  } else {
    message(paste0("WARNING: Generating weighted-probabilistic ensemble assignments from ", length(sample_list_2), " Sample IDs. This is not the number of pooled samples defined by the user."))
  }
  
  ## Compute weighted probability for each sample
  for (i in sample_list_2) {
    # Doublets
    practice_df$doublet <- ifelse(practice_df$vireo_best_assignment == "doublet", practice_df$vireo_weighted_probability, 0)
    practice_df$doublet <- ifelse(practice_df$souporcell_best_assignment == "doublet", practice_df$doublet + practice_df$souporcell_weighted_probability, practice_df$doublet+ 0)
    practice_df$doublet <- ifelse(practice_df$demuxalot_best_assignment == "doublet", practice_df$doublet + practice_df$demuxalot_weighted_probability, practice_df$doublet+ 0)
    practice_df$doublet <- ifelse(practice_df$demuxlet_best_assignment == "doublet", practice_df$doublet + practice_df$demuxlet_weighted_probability, practice_df$doublet+ 0)
    
    # Singlets
    practice_df$sample <- ifelse(practice_df$vireo_best_assignment == i, practice_df$vireo_weighted_probability, 0)
    practice_df$sample <- ifelse(practice_df$souporcell_best_assignment == i, practice_df$sample + practice_df$souporcell_weighted_probability, practice_df$sample+ 0)
    practice_df$sample <- ifelse(practice_df$demuxalot_best_assignment == i, practice_df$sample + practice_df$demuxalot_weighted_probability, practice_df$sample+ 0)
    practice_df$sample <- ifelse(practice_df$demuxlet_best_assignment == i, practice_df$sample + practice_df$demuxlet_weighted_probability, practice_df$sample+ 0)
    colnames(practice_df)[ncol(practice_df)] <- i  
  }
  
  ## Select sample assignment with maximum weighted probability for each droplet
  data_long <- gather(practice_df, key="ensemblex_assignment", value="stat", 32:ncol(practice_df)) 
  result <- data_long %>% 
    dplyr::group_by(barcode) %>%
    dplyr::filter(stat == max(stat))
  
  ## Get remaining probabilities for non-assigned samples
  data_long <- gather(practice_df, key="ensemblex_assignment", value="stat", 32:ncol(practice_df))
  result_sum <- data_long %>% 
    dplyr::group_by(barcode) %>%
    dplyr::summarise(total = sum(stat))
  result_test <- merge(result, result_sum, by = "barcode")
  
  ## Calculate ensemblex probability 
  result_test$ensemblex_probability <- result_test$stat/result_test$total
  
  ## Clean up dataframe
  result_test <- dplyr::select(result_test, c("barcode", "ensemblex_assignment", "ensemblex_probability", "vireo_assignment", "souporcell_assignment", "demuxlet_assignment", "demuxalot_assignment", "general_consensus",                 
                                              "vireo_best_assignment", "souporcell_best_assignment", "demuxlet_best_assignment", "demuxalot_best_assignment", "vireo_singlet_probability", "vireo_doublet_probability",         
                                              "vireo_n_vars", "vireo_best_doublet", "vireo_doublet_logLikRatio", "souporcell_log_probability_singlet", "souporcell_log_probability_doublet", "demuxlet_n_snps", 
                                              "demuxlet_n_reads", "demuxlet_max_probability", "demuxlet_DIFF_LLK_SNG_DBL", "demuxalot_max_probability", "demuxalot_doublet_probability", "vireo_max_probability",
                                             "vireo_weighted_probability", "demuxlet_weighted_probability", 
                                              "demuxalot_weighted_probability", "souporcell_weighted_probability"))

  ## save PWE assignment dataframe
  write.csv(result_test, paste(par_output_dir,"/step1",'/step1_cell_assignment.csv', sep=""))
  result_test
}


###########################################################################################################################
# GRAPH-BASED DOUBLET DETECTION
###########################################################################################################################
graph_based_doublet_detection_par_sweep <- function(result_test, par_expected_doublet_rate, par_output_dir){
    
    ## Set seed
    set.seed(1234)
    
    ## create an output directory
    dir.create(paste(par_output_dir,"/step2",sep=''))
    
    ## load Balanced-accuracy dataset
    result_2 <- result_test
    
    ### Perform principal component analysis with select variables
    result_2_lim <- result_2
    result_2_lim <- result_2_lim[,c(1, 14, 19, 20, 21, 17, 23, 25, 13    )] #barcode, vireo_doublet_probability,  souporcell_log_probability_doublet, demuxlet_n_snps, demuxlet_n_reads, vireo_doublet_logLikRatio, demuxlet_DIFF_LLK_SNG_DBL, demuxalot_doublet_probability, vireo_singlet_probability (we dont use this for PCA)
    result_2_lim <- result_2_lim[complete.cases(result_2_lim), ]
    res.pca <- prcomp(result_2_lim[,c(2:8)],scale = T)
    
    ## scree plot
    fviz_eig(res.pca)
    ggsave(paste(par_output_dir,"/step2","/PCA_scree_plot.pdf", sep=""))
    
    ## PCA
    fviz_pca_ind(res.pca,
                col.ind = "black", 
                geom="point", 
                pointsize = 0.5
    )
    ggsave(paste(par_output_dir,"/step2","/PCA_plot.pdf", sep=""))
    
    ### variable contribution to variation
    ## Plot contributions of variables to PC1
    fviz_contrib(res.pca, choice = "var", axes = 1, top = 10) + theme_classic() + coord_flip() + xlab("Variable") + ggtitle("PC1")
    ggsave(paste(par_output_dir,"/step2","/PC1_var_contrib.pdf", sep=""))

    ## Plot contributions of variables to PC2
    fviz_contrib(res.pca, choice = "var", axes = 2, top = 10) + theme_classic() + coord_flip() + xlab("Variable") + ggtitle("PC2")
    ggsave(paste(par_output_dir,"/step2","/PC2_var_contrib.pdf", sep=""))
        
    ## Compute euclidean distance between points
    rownames(result_2_lim) <- result_2_lim$barcode
    res.pca <- prcomp(result_2_lim[,c(2,3,4,5,6,7,8)],scale = T) 
    df_1 <- data.frame(res.pca$x[,1])
    df_1$barcode <- rownames(df_1)
    df_2 <- data.frame(res.pca$x[,2])
    df_2$barcode <- rownames(df_2)
    df_merge_test <- merge(df_1, df_2, by = "barcode")
    colnames(df_merge_test) <- c("barcode", "PC1", "PC2")
    distances <- dist(df_merge_test[c("PC1", "PC2")], diag = TRUE, upper = TRUE)
    distances <- as.matrix(distances)
    colnames(distances) <- df_merge_test$barcode
    rownames(distances) <- df_merge_test$barcode
    
    ### Organize parameters to identify most likely doublets
    ## Organize Vireo doublet log_lik into ordered frame
    vireo_doublet_df <- result_2 %>% select("barcode", "vireo_doublet_logLikRatio")
    vireo_doublet_df <- vireo_doublet_df %>% arrange(desc(vireo_doublet_logLikRatio))
    
    ## Organize Vireo doublet probability
    vireo_doublet_2_df <- result_2 %>% select("barcode", "vireo_doublet_probability")
    vireo_doublet_2_df <- vireo_doublet_2_df %>% arrange(desc(vireo_doublet_probability))
    
    ## Organize Demuxlet "freemuxlet_DIFF_LLK_SNG_DBL" in ordered frame
    freemuxlet_doublet_df <- result_2 %>% select("barcode", "demuxlet_DIFF_LLK_SNG_DBL")
    freemuxlet_doublet_df <- freemuxlet_doublet_df %>% arrange(demuxlet_DIFF_LLK_SNG_DBL)
    
    ## Organize Souporcell log prob doublet 
    souporcell_doublet_df <- result_2 %>% select("barcode", "souporcell_log_probability_doublet")
    souporcell_doublet_df <- souporcell_doublet_df %>% arrange(souporcell_log_probability_doublet)
    
    ## Organize Demuxlet num snp
    freemuxlet_snp_df <- result_2 %>% select("barcode", "demuxlet_n_snps")
    freemuxlet_snp_df <- freemuxlet_snp_df %>% arrange(desc(demuxlet_n_snps))
    
    ## Organize Demuxlet num reads
    freemuxlet_reads_df <- result_2 %>% select("barcode", "demuxlet_n_reads")
    freemuxlet_reads_df <- freemuxlet_reads_df %>% arrange(desc(demuxlet_n_reads))
    
    ## Organize Demuxalot probability
    demuxalot_doublet_df <- result_2 %>% select("barcode", "demuxalot_doublet_probability")
    demuxalot_doublet_df <- demuxalot_doublet_df %>% arrange(desc(demuxalot_doublet_probability))
    
    ### Organize metrics by percentile
    ## Vireo diff 
    vireo_doublet_df <- vireo_doublet_df %>% 
        mutate(percentile  = percent_rank(vireo_doublet_logLikRatio))
    
    ## Vireo doublet probability
    vireo_doublet_2_df <- vireo_doublet_2_df %>% 
        mutate(percentile  = percent_rank(vireo_doublet_probability))
    
    ## Demuxlet
    freemuxlet_doublet_df <- freemuxlet_doublet_df %>% 
        mutate(percentile  = percent_rank(demuxlet_DIFF_LLK_SNG_DBL))
    freemuxlet_doublet_df$percentile <- 1 - freemuxlet_doublet_df$percentile
    
    ## Demuxlet num reads
    freemuxlet_reads_df <- freemuxlet_reads_df %>% 
        mutate(percentile  = percent_rank(demuxlet_n_reads))
    freemuxlet_reads_df$percentile <-  freemuxlet_reads_df$percentile
    
    ## Demuxlet num snps
    freemuxlet_snp_df <- freemuxlet_snp_df %>% 
        mutate(percentile  = percent_rank(demuxlet_n_snps))
    freemuxlet_snp_df$percentile <-  freemuxlet_snp_df$percentile
    
    ## Souporcell doublet prob
    souporcell_doublet_df <- souporcell_doublet_df %>% 
        mutate(percentile  = percent_rank(souporcell_log_probability_doublet))
    souporcell_doublet_df$percentile <- 1 - souporcell_doublet_df$percentile
    
    ## Demuxalot doublet prob
    demuxalot_doublet_df <- demuxalot_doublet_df %>% 
        mutate(percentile  = percent_rank(demuxalot_doublet_probability))
    demuxalot_doublet_df$percentile <- demuxalot_doublet_df$percentile
    
    ################################################################################################################################
    #### parameter sweep #### 
    ################################################################################################################################
    ### Compute varying number of confident doublets (nCD)
    ## 50 nCD 
    suspected_doublets_50 <- list()
    s1 <- rev(seq(0, 1, by = 0.01))
    for (j in s1){
        vireo_doublet_percentile <- subset(vireo_doublet_df, percentile >= j)
        vireo_doublet_2_percentile <- subset(vireo_doublet_2_df, percentile >= j)
        freemuxlet_doublet_percentile <- subset(freemuxlet_doublet_df, percentile >= j)
        souporcell_doublet_percentile <- subset(souporcell_doublet_df, percentile >= j)
        freemuxlet_snp_percentile <- subset(freemuxlet_snp_df, percentile >= j)
        freemuxlet_reads_percentile <- subset(freemuxlet_reads_df, percentile >= j)
        demuxalot_doublet_percentile <- subset(demuxalot_doublet_df, percentile >= j)
        
        for (i in result_2$barcode){
        if (i %in% suspected_doublets_50){
            print("Barcode already doublet")
        } else if (i %in% vireo_doublet_percentile$barcode &
                    i %in% vireo_doublet_2_percentile$barcode &
                    i %in% freemuxlet_doublet_percentile$barcode &
                    i %in% souporcell_doublet_percentile$barcode &
                    i %in% freemuxlet_snp_percentile$barcode &
                    i %in% freemuxlet_reads_percentile$barcode &
                    i %in% demuxalot_doublet_percentile$barcode) {
            print(paste0("adding barcode to doublet list at ", j, "percentile"))
            suspected_doublets_50 <- c(suspected_doublets_50, i)
        } else {
            print(paste0("barcode not a suspected doublet at", j, "percentile"))
        }
        suspected_doublets_50
        }
        if(length(suspected_doublets_50) >= 50 ) { 
        print(paste0("reached expected doublet rate at", j, "percentile"))
        break
        }
        suspected_doublets_50
    }  
    
    ## 100 nCD
    suspected_doublets_100 <- list()
    s1 <- rev(seq(0, 1, by = 0.01))
    for (j in s1){
        vireo_doublet_percentile <- subset(vireo_doublet_df, percentile >= j)
        vireo_doublet_2_percentile <- subset(vireo_doublet_2_df, percentile >= j)
        freemuxlet_doublet_percentile <- subset(freemuxlet_doublet_df, percentile >= j)
        souporcell_doublet_percentile <- subset(souporcell_doublet_df, percentile >= j)
        freemuxlet_snp_percentile <- subset(freemuxlet_snp_df, percentile >= j)
        freemuxlet_reads_percentile <- subset(freemuxlet_reads_df, percentile >= j)
        demuxalot_doublet_percentile <- subset(demuxalot_doublet_df, percentile >= j)
        
        for (i in result_2$barcode){
        if (i %in% suspected_doublets_100){
            print("Barcode already doublet")
        } else if (i %in% vireo_doublet_percentile$barcode &
                    i %in% vireo_doublet_2_percentile$barcode &
                    i %in% freemuxlet_doublet_percentile$barcode &
                    i %in% souporcell_doublet_percentile$barcode &
                    i %in% freemuxlet_snp_percentile$barcode &
                    i %in% freemuxlet_reads_percentile$barcode &
                    i %in% demuxalot_doublet_percentile$barcode) {
            print(paste0("adding barcode to doublet list at ", j, "percentile"))
            suspected_doublets_100 <- c(suspected_doublets_100, i)
        } else {
            print(paste0("barcode not a suspected doublet at", j, "percentile"))
        }
        suspected_doublets_100
        }
        if(length(suspected_doublets_100) >= 100 ) { 
        print(paste0("reached expected doublet rate at", j, "percentile"))
        break
        }
        suspected_doublets_100
    }
    
    ## 150 nCD
    suspected_doublets_150 <- list()
    s1 <- rev(seq(0, 1, by = 0.01))
    for (j in s1){
        vireo_doublet_percentile <- subset(vireo_doublet_df, percentile >= j)
        vireo_doublet_2_percentile <- subset(vireo_doublet_2_df, percentile >= j)
        freemuxlet_doublet_percentile <- subset(freemuxlet_doublet_df, percentile >= j)
        souporcell_doublet_percentile <- subset(souporcell_doublet_df, percentile >= j)
        freemuxlet_snp_percentile <- subset(freemuxlet_snp_df, percentile >= j)
        freemuxlet_reads_percentile <- subset(freemuxlet_reads_df, percentile >= j)
        demuxalot_doublet_percentile <- subset(demuxalot_doublet_df, percentile >= j)
        
        for (i in result_2$barcode){
        if (i %in% suspected_doublets_150){
            print("Barcode already doublet")
        } else if (i %in% vireo_doublet_percentile$barcode &
                    i %in% vireo_doublet_2_percentile$barcode &
                    i %in% freemuxlet_doublet_percentile$barcode &
                    i %in% souporcell_doublet_percentile$barcode &
                    i %in% freemuxlet_snp_percentile$barcode &
                    i %in% freemuxlet_reads_percentile$barcode &
                    i %in% demuxalot_doublet_percentile$barcode) {
            print(paste0("adding barcode to doublet list at ", j, "percentile"))
            suspected_doublets_150 <- c(suspected_doublets_150, i)
        } else {
            print(paste0("barcode not a suspected doublet at", j, "percentile"))
        }
        suspected_doublets_150
        }
        if(length(suspected_doublets_150) >= 150 ) { #OG = 100
        print(paste0("reached expected doublet rate at", j, "percentile"))
        break
        }
        suspected_doublets_150
    }
    
    ## 200 nCD
    suspected_doublets_200 <- list()
    s1 <- rev(seq(0, 1, by = 0.01))
    for (j in s1){
        vireo_doublet_percentile <- subset(vireo_doublet_df, percentile >= j)
        vireo_doublet_2_percentile <- subset(vireo_doublet_2_df, percentile >= j)
        freemuxlet_doublet_percentile <- subset(freemuxlet_doublet_df, percentile >= j)
        souporcell_doublet_percentile <- subset(souporcell_doublet_df, percentile >= j)
        freemuxlet_snp_percentile <- subset(freemuxlet_snp_df, percentile >= j)
        freemuxlet_reads_percentile <- subset(freemuxlet_reads_df, percentile >= j)
        demuxalot_doublet_percentile <- subset(demuxalot_doublet_df, percentile >= j)
        
        for (i in result_2$barcode){
        if (i %in% suspected_doublets_200){
            print("Barcode already doublet")
        } else if (i %in% vireo_doublet_percentile$barcode &
                    i %in% vireo_doublet_2_percentile$barcode &
                    i %in% freemuxlet_doublet_percentile$barcode &
                    i %in% souporcell_doublet_percentile$barcode &
                    i %in% freemuxlet_snp_percentile$barcode &
                    i %in% freemuxlet_reads_percentile$barcode &
                    i %in% demuxalot_doublet_percentile$barcode) {
            print(paste0("adding barcode to doublet list at ", j, "percentile"))
            suspected_doublets_200 <- c(suspected_doublets_200, i)
        } else {
            print(paste0("barcode not a suspected doublet at", j, "percentile"))
        }
        suspected_doublets_200
        }
        if(length(suspected_doublets_200) >= 200 ) { #OG = 100
        print(paste0("reached expected doublet rate at", j, "percentile"))
        break
        }
        suspected_doublets_200
    }
    
    ## 250 nCD 
    suspected_doublets_250 <- list()
    s1 <- rev(seq(0, 1, by = 0.01))
    for (j in s1){
        vireo_doublet_percentile <- subset(vireo_doublet_df, percentile >= j)
        vireo_doublet_2_percentile <- subset(vireo_doublet_2_df, percentile >= j)
        freemuxlet_doublet_percentile <- subset(freemuxlet_doublet_df, percentile >= j)
        souporcell_doublet_percentile <- subset(souporcell_doublet_df, percentile >= j)
        freemuxlet_snp_percentile <- subset(freemuxlet_snp_df, percentile >= j)
        freemuxlet_reads_percentile <- subset(freemuxlet_reads_df, percentile >= j)
        demuxalot_doublet_percentile <- subset(demuxalot_doublet_df, percentile >= j)
        
        for (i in result_2$barcode){
        if (i %in% suspected_doublets_250){
            print("Barcode already doublet")
        } else if (i %in% vireo_doublet_percentile$barcode &
                    i %in% vireo_doublet_2_percentile$barcode &
                    i %in% freemuxlet_doublet_percentile$barcode &
                    i %in% souporcell_doublet_percentile$barcode &
                    i %in% freemuxlet_snp_percentile$barcode &
                    i %in% freemuxlet_reads_percentile$barcode &
                    i %in% demuxalot_doublet_percentile$barcode) {
            print(paste0("adding barcode to doublet list at ", j, "percentile"))
            suspected_doublets_250 <- c(suspected_doublets_250, i)
        } else {
            print(paste0("barcode not a suspected doublet at", j, "percentile"))
        }
        suspected_doublets_250
        }
        if(length(suspected_doublets_250) >= 250 ) { #OG = 100
        print(paste0("reached expected doublet rate at", j, "percentile"))
        break
        }
        suspected_doublets_250
    }
    
    ## 300 nCD
    suspected_doublets_300 <- list()
    s1 <- rev(seq(0, 1, by = 0.01))
    for (j in s1){
        vireo_doublet_percentile <- subset(vireo_doublet_df, percentile >= j)
        vireo_doublet_2_percentile <- subset(vireo_doublet_2_df, percentile >= j)
        freemuxlet_doublet_percentile <- subset(freemuxlet_doublet_df, percentile >= j)
        souporcell_doublet_percentile <- subset(souporcell_doublet_df, percentile >= j)
        freemuxlet_snp_percentile <- subset(freemuxlet_snp_df, percentile >= j)
        freemuxlet_reads_percentile <- subset(freemuxlet_reads_df, percentile >= j)
        demuxalot_doublet_percentile <- subset(demuxalot_doublet_df, percentile >= j)
        
        for (i in result_2$barcode){
        if (i %in% suspected_doublets_300){
            print("Barcode already doublet")
        } else if (i %in% vireo_doublet_percentile$barcode &
                    i %in% vireo_doublet_2_percentile$barcode &
                    i %in% freemuxlet_doublet_percentile$barcode &
                    i %in% souporcell_doublet_percentile$barcode &
                    i %in% freemuxlet_snp_percentile$barcode &
                    i %in% freemuxlet_reads_percentile$barcode &
                    i %in% demuxalot_doublet_percentile$barcode) {
            print(paste0("adding barcode to doublet list at ", j, "percentile"))
            suspected_doublets_300 <- c(suspected_doublets_300, i)
        } else {
            print(paste0("barcode not a suspected doublet at", j, "percentile"))
        }
        suspected_doublets_300
        }
        if(length(suspected_doublets_300) >= 300 ) { #OG = 100
        print(paste0("reached expected doublet rate at", j, "percentile"))
        break
        }
        suspected_doublets_300
    }
    

    ## Make a list of percentile threshold (pT) values based on expected doublet rate for the pool
    message(paste0("Expected doublet rate is ", par_expected_doublet_rate, "; ", round(par_expected_doublet_rate*nrow(result_2), digits = 0), " droplets"  ))
    interval <- par_expected_doublet_rate/6
    pT_list <- rev(seq(1-par_expected_doublet_rate,1-interval, by = interval))
        
    ## 50 nCD pT sweep 
    doublet <- "none"
    percentile <- 0
    barcode <- "none"
    pT <- 0
    filler_frame_50 <- data.frame(doublet,percentile, barcode,pT )

    for (t in pT_list){
        for (j in suspected_doublets_50){
        distances_test_1 <-distances
        distances_test_1 <- data.frame(distances_test_1[,j])
        colnames(distances_test_1) <- "doublet"

        distances_test_1 <- distances_test_1 %>% 
        mutate(percentile  = percent_rank(doublet))
        distances_test_1$percentile <- 1 - distances_test_1$percentile
        
        distances_test_15 <- distances_test_1[distances_test_1$percentile >= t,] 
        distances_test_15$barcode <- rownames(distances_test_15)
        distances_test_15$doublet <- "temp"
        distances_test_15$pT <- t
        
        filler_frame_50 <- rbind(filler_frame_50,distances_test_15)
        filler_frame_50  
        }
        filler_frame_50
    }
    
    ## 100 nCD pT sweep 
    doublet <- "none"
    percentile <- 0
    barcode <- "none"
    pT <- 0
    filler_frame_100 <- data.frame(doublet,percentile, barcode,pT )
    
    for (t in pT_list){
        for (j in suspected_doublets_100){
        distances_test_1 <-distances
        distances_test_1 <- data.frame(distances_test_1[,j])
        colnames(distances_test_1) <- "doublet"
        
        distances_test_1 <- distances_test_1 %>% 
        mutate(percentile  = percent_rank(doublet))
        distances_test_1$percentile <- 1 - distances_test_1$percentile
        
        distances_test_15 <- distances_test_1[distances_test_1$percentile >= t,] 
        distances_test_15$barcode <- rownames(distances_test_15)
        distances_test_15$doublet <- "temp"
        distances_test_15$pT <- t
        
        filler_frame_100 <- rbind(filler_frame_100,distances_test_15)
        filler_frame_100 
        }
        filler_frame_100
    }
    
    ## 150 nCD pT sweep
    doublet <- "none"
    percentile <- 0
    barcode <- "none"
    pT <- 0
    filler_frame_150 <- data.frame(doublet,percentile, barcode,pT )
    
    for (t in pT_list){
        for (j in suspected_doublets_150){
        distances_test_1 <-distances
        distances_test_1 <- data.frame(distances_test_1[,j])
        colnames(distances_test_1) <- "doublet"
        
        distances_test_1 <- distances_test_1 %>% 
        mutate(percentile  = percent_rank(doublet))
        distances_test_1$percentile <- 1 - distances_test_1$percentile
        
        distances_test_15 <- distances_test_1[distances_test_1$percentile >= t,]
        distances_test_15$barcode <- rownames(distances_test_15)
        distances_test_15$doublet <- "temp"
        distances_test_15$pT <- t
        
        filler_frame_150 <- rbind(filler_frame_150,distances_test_15)
        filler_frame_150  
        }
        filler_frame_150
    }
    
    ## 200 nCD pT sweep 
    doublet <- "none"
    percentile <- 0
    barcode <- "none"
    pT <- 0
    filler_frame_200 <- data.frame(doublet,percentile, barcode,pT )
    
    for (t in pT_list){
        for (j in suspected_doublets_200){
        distances_test_1 <-distances
        distances_test_1 <- data.frame(distances_test_1[,j])
        colnames(distances_test_1) <- "doublet"
        
        distances_test_1 <- distances_test_1 %>% 
        mutate(percentile  = percent_rank(doublet))
        distances_test_1$percentile <- 1 - distances_test_1$percentile
        
        distances_test_15 <- distances_test_1[distances_test_1$percentile >= t,] 
        distances_test_15$barcode <- rownames(distances_test_15)
        distances_test_15$doublet <- "temp"
        distances_test_15$pT <- t
        
        filler_frame_200 <- rbind(filler_frame_200,distances_test_15)
        filler_frame_200
        }
        filler_frame_200
    }
    
    ## 250 nCD pT sweep 
    doublet <- "none"
    percentile <- 0
    barcode <- "none"
    pT <- 0
    filler_frame_250 <- data.frame(doublet,percentile, barcode,pT )
    
    for (t in pT_list){
        for (j in suspected_doublets_250){
        distances_test_1 <-distances
        distances_test_1 <- data.frame(distances_test_1[,j])
        colnames(distances_test_1) <- "doublet"
        
        distances_test_1 <- distances_test_1 %>% 
        mutate(percentile  = percent_rank(doublet))
        distances_test_1$percentile <- 1 - distances_test_1$percentile
        
        distances_test_15 <- distances_test_1[distances_test_1$percentile >= t,]
        distances_test_15$barcode <- rownames(distances_test_15)
        distances_test_15$doublet <- "temp"
        distances_test_15$pT <- t
        
        filler_frame_250 <- rbind(filler_frame_250,distances_test_15)
        filler_frame_250
        }
        filler_frame_250
    }
    
    ## 300 nCD pT sweep 
    doublet <- "none"
    percentile <- 0
    barcode <- "none"
    pT <- 0
    filler_frame_300 <- data.frame(doublet,percentile, barcode,pT )
    
    for (t in pT_list){
        for (j in suspected_doublets_300){
        distances_test_1 <-distances
        distances_test_1 <- data.frame(distances_test_1[,j])
        colnames(distances_test_1) <- "doublet"
        
        distances_test_1 <- distances_test_1 %>% 
        mutate(percentile  = percent_rank(doublet))
        distances_test_1$percentile <- 1 - distances_test_1$percentile
        
        distances_test_15 <- distances_test_1[distances_test_1$percentile >= t,] 
        distances_test_15$barcode <- rownames(distances_test_15)
        distances_test_15$doublet <- "temp"
        distances_test_15$pT <- t
        
        filler_frame_300 <- rbind(filler_frame_300,distances_test_15)
        filler_frame_300
        }
        filler_frame_300
    }
    
    ### Compute nearest neighbour frequency and Kutosis of frequency distributions
    ## 50 nCD 
    counts_50 <- filler_frame_50 %>%
        group_by(barcode,pT) %>%
        dplyr::mutate(n = n()) %>%
        ungroup()
    
    counts_50 <- counts_50 %>% 
        distinct(pT, barcode, .keep_all = TRUE)
    
    ncD <- 50
    pT <- 0
    kurtosis <- 0
    fill_frame_50 <- data.frame(ncD,pT, kurtosis)
    
    for (t in pT_list){
        k <- kurtosis(counts_50$n[counts_50$pT == t])
        nCD <- 50
        pT <- t
        kurtosis <- k
        temp_frame <- data.frame(ncD,pT, kurtosis)
        fill_frame_50 <- rbind(fill_frame_50,temp_frame)
        fill_frame_50 <- subset(fill_frame_50, pT != 0)
        fill_frame_50
    }
    
    ## 100 nCD  
    counts_100 <- filler_frame_100 %>%
        group_by(barcode,pT) %>%
        dplyr::mutate(n = n()) %>%
        ungroup()
    
    counts_100 <- counts_100 %>% 
        distinct(pT, barcode, .keep_all = TRUE)
    
    ncD <- 100
    pT <- 0
    kurtosis <- 0
    fill_frame_100 <- data.frame(ncD,pT, kurtosis)
    
    for (t in pT_list){
        k <- kurtosis(counts_100$n[counts_100$pT == t])
        nCD <- 100
        pT <- t
        kurtosis <- k
        temp_frame <- data.frame(ncD,pT, kurtosis)
        fill_frame_100 <- rbind(fill_frame_100,temp_frame)
        fill_frame_100 <- subset(fill_frame_100, pT != 0)
        fill_frame_100
    }
    
    ## 150 nCD 
    counts_150 <- filler_frame_150 %>%
        group_by(barcode,pT) %>%
        dplyr::mutate(n = n()) %>%
        ungroup()
    
    counts_150 <- counts_150 %>% 
        distinct(pT, barcode, .keep_all = TRUE)
    
    ncD <- 150
    pT <- 0
    kurtosis <- 0
    fill_frame_150 <- data.frame(ncD,pT, kurtosis)
    
    for (t in pT_list){
        k <- kurtosis(counts_150$n[counts_150$pT == t])
        nCD <- 150
        pT <- t
        kurtosis <- k
        temp_frame <- data.frame(ncD,pT, kurtosis)
        fill_frame_150 <- rbind(fill_frame_150,temp_frame)
        fill_frame_150 <- subset(fill_frame_150, pT != 0)
        fill_frame_150
    }
    
    ## 200 nCD  
    counts_200 <- filler_frame_200 %>%
        group_by(barcode,pT) %>%
        dplyr::mutate(n = n()) %>%
        ungroup()
    
    counts_200 <- counts_200 %>% 
        distinct(pT, barcode, .keep_all = TRUE)
    
    ncD <- 200
    pT <- 0
    kurtosis <- 0
    fill_frame_200 <- data.frame(ncD,pT, kurtosis)
    
    for (t in pT_list){
        k <- kurtosis(counts_200$n[counts_200$pT == t])
        nCD <- 200
        pT <- t
        kurtosis <- k
        temp_frame <- data.frame(ncD,pT, kurtosis)
        fill_frame_200 <- rbind(fill_frame_200,temp_frame)
        fill_frame_200 <- subset(fill_frame_200, pT != 0)
        fill_frame_200
    }
    
    ## 250 nCD 
    counts_250 <- filler_frame_250 %>%
        group_by(barcode,pT) %>%
        dplyr::mutate(n = n()) %>%
        ungroup()

    counts_250 <- counts_250 %>% 
        distinct(pT, barcode, .keep_all = TRUE)
    
    ncD <- 250
    pT <- 0
    kurtosis <- 0
    fill_frame_250 <- data.frame(ncD,pT, kurtosis)
    
    for (t in pT_list){
        k <- kurtosis(counts_250$n[counts_250$pT == t])
        nCD <- 250
        pT <- t
        kurtosis <- k
        temp_frame <- data.frame(ncD,pT, kurtosis)
        fill_frame_250 <- rbind(fill_frame_250,temp_frame)
        fill_frame_250 <- subset(fill_frame_250, pT != 0)
        fill_frame_250
    }
    
    ## 300 nCD 
    counts_300 <- filler_frame_300 %>%
        group_by(barcode,pT) %>%
        dplyr::mutate(n = n()) %>%
        ungroup()

    counts_300 <- counts_300 %>% 
        distinct(pT, barcode, .keep_all = TRUE)
    
    ncD <- 300
    pT <- 0
    kurtosis <- 0
    fill_frame_300 <- data.frame(ncD,pT, kurtosis)
    
    for (t in pT_list){
        k <- kurtosis(counts_300$n[counts_300$pT == t])
        nCD <- 300
        pT <- t
        kurtosis <- k
        temp_frame <- data.frame(ncD,pT, kurtosis)
        fill_frame_300 <- rbind(fill_frame_300,temp_frame)
        fill_frame_300 <- subset(fill_frame_300, pT != 0)
        fill_frame_300
    }
    
    ### Compute optimal parameters 
    ## Bind all of the Kurtosis frames 
    kurtosis_bind <- rbind(fill_frame_50,fill_frame_100,fill_frame_150,fill_frame_200,fill_frame_250, fill_frame_300)
    
    ## Compute average Kurtosis for each pT
    kurtosis_bind <- kurtosis_bind %>%
        dplyr::group_by(pT) %>%
        dplyr::summarize(Mean_k = mean(kurtosis, na.rm=TRUE))
    
    optimal_pT <- kurtosis_bind$pT[kurtosis_bind$Mean_k == max(kurtosis_bind$Mean_k)] 
    
    ## Plot kurtosis values
    ggplot(kurtosis_bind, aes(x = pT, y = Mean_k)) + 
        geom_vline(xintercept = optimal_pT, col = "red") +
        geom_point() +  
        geom_line() + 
        theme_classic() +
        ggtitle(paste0("Estimated optimal pT: ", optimal_pT))
    ggsave(paste(par_output_dir,"/step2","/optimal_pT.pdf", sep=""))
    
    ## find optimal nCD
    bind_nCD <- rbind(fill_frame_50,fill_frame_100,fill_frame_150,fill_frame_200,fill_frame_250, fill_frame_300)
    bind_nCD <- subset(bind_nCD, pT == optimal_pT)
    ncD <- data.frame(nCD = c(50, 100, 150, 200, 250, 300))
    ncD_list <- c(50, 100, 150, 200, 250, 300)
    
    ## Smooth line
    fm1 <- smooth.spline(bind_nCD[,"ncD"], bind_nCD[,"kurtosis"], df = 3)
    y2 <- predict(fm1, x = ncD)
    y <- data.frame(y2$y)
    y <- y$nCD

    ## Find elbow of the smoothed curve 
    df <- data.frame(ncd = ncD_list, kurtosis = y)
    optimal_nCD <- find_curve_elbow(df, export_type = "row_num", plot_curve = FALSE)

    ## parse value
    optimal_nCD <- df[optimal_nCD, 1]

    ## If cannot find slope, take point preceedinging the largest slope
    if (is.na(optimal_nCD)) {
        message("Could not determine optimal nCD based on the elbow, taking maximum kurtosis value.")
        max <- max(bind_nCD$kurtosis)
        optimal_nCD <- bind_nCD$ncD[bind_nCD$kurtosis == max]
        optimal_nCD <- optimal_nCD[1]
        #slope= diff(bind_nCD$ncD)/diff(bind_nCD$kurtosis)
        #xx <- bind_nCD[which.max(slope),]
        #optimal_nCD <- max(xx$ncD)
    }else{
        optimal_nCD <- optimal_nCD
        message("Determine optimal nCD based on the elbow.")
    }

    ## plot 
    ggplot(bind_nCD, aes(x = ncD, y = kurtosis)) + 
        geom_vline(xintercept = optimal_nCD, col = "red") +
        geom_point() +  
        geom_line() + 
        theme_classic() +
        ggtitle(paste0("Estimated optimal nCD: ", optimal_nCD))
    ggsave(paste(par_output_dir,"/step2","/optimal_nCD.pdf", sep=""))
    
}
graph_based_doublet_detection_manual_par <- function(result_test, par_expected_doublet_rate, par_output_dir,par_ensemblex_pT,par_ensemblex_nCD){
    
    ## Set seed
    set.seed(1234)
    
    ## create an output directory
    dir.create(paste(par_output_dir,"/step2",sep=''))
    
    ## load manually defined optimal parameters
    optimal_pT <- par_ensemblex_pT
    optimal_nCD <- par_ensemblex_nCD
    
    ## load Balanced-accuracy dataset
    result_2 <- result_test

    ### Perform principal component analysis with select variables
    result_2_lim <- result_2
    result_2_lim <- result_2_lim[,c(1, 14, 19, 20, 21, 17, 23, 25, 13 )] #barcode, vireo_doublet_probability,  souporcell_log_probability_doublet, demuxlet_n_snps, demuxlet_n_reads, vireo_doublet_logLikRatio, demuxlet_DIFF_LLK_SNG_DBL, demuxalot_doublet_probability, vireo_singlet_probability (we dont use this for PCA)
    result_2_lim <- result_2_lim[complete.cases(result_2_lim), ]
    res.pca <- prcomp(result_2_lim[,c(2:8)],scale = T)
    
    ## scree plot
    fviz_eig(res.pca)
    ggsave(paste(par_output_dir,"/step2","/PCA_scree_plot.pdf", sep=""))
    
    ## PCA
    fviz_pca_ind(res.pca,
                col.ind = "black", 
                geom="point", 
                pointsize = 0.5
    )
    ggsave(paste(par_output_dir,"/step2","/PCA_plot.pdf", sep=""))
    
    ### variable contribution to variation
    ## Plot contributions of variables to PC1
    fviz_contrib(res.pca, choice = "var", axes = 1, top = 10) + theme_classic() + coord_flip() + xlab("Variable") + ggtitle("PC1")
    ggsave(paste(par_output_dir,"/step2","/PC1_var_contrib.pdf", sep=""))

    ## Plot contributions of variables to PC2
    fviz_contrib(res.pca, choice = "var", axes = 2, top = 10) + theme_classic() + coord_flip() + xlab("Variable") + ggtitle("PC2")
    ggsave(paste(par_output_dir,"/step2","/PC2_var_contrib.pdf", sep=""))
        
    ## Compute euclidean distance between points
    rownames(result_2_lim) <- result_2_lim$barcode
    res.pca <- prcomp(result_2_lim[,c(2,3,4,5,6,7,8)],scale = T) 
    df_1 <- data.frame(res.pca$x[,1])
    df_1$barcode <- rownames(df_1)
    df_2 <- data.frame(res.pca$x[,2])
    df_2$barcode <- rownames(df_2)
    df_merge_test <- merge(df_1, df_2, by = "barcode")
    colnames(df_merge_test) <- c("barcode", "PC1", "PC2")
    distances <- dist(df_merge_test[c("PC1", "PC2")], diag = TRUE, upper = TRUE)
    distances <- as.matrix(distances)
    colnames(distances) <- df_merge_test$barcode
    rownames(distances) <- df_merge_test$barcode
    
    ### Organize parameters to identify most likely doublets
    ## Organize Vireo doublet log_lik into ordered frame
    vireo_doublet_df <- result_2 %>% select("barcode", "vireo_doublet_logLikRatio")
    vireo_doublet_df <- vireo_doublet_df %>% arrange(desc(vireo_doublet_logLikRatio))
    
    ## Organize Vireo doublet probability
    vireo_doublet_2_df <- result_2 %>% select("barcode", "vireo_doublet_probability")
    vireo_doublet_2_df <- vireo_doublet_2_df %>% arrange(desc(vireo_doublet_probability))
    
    ## Organize Demuxlet "freemuxlet_DIFF_LLK_SNG_DBL" in ordered frame
    freemuxlet_doublet_df <- result_2 %>% select("barcode", "demuxlet_DIFF_LLK_SNG_DBL")
    freemuxlet_doublet_df <- freemuxlet_doublet_df %>% arrange(demuxlet_DIFF_LLK_SNG_DBL)
    
    ## Organize Souporcell log prob doublet 
    souporcell_doublet_df <- result_2 %>% select("barcode", "souporcell_log_probability_doublet")
    souporcell_doublet_df <- souporcell_doublet_df %>% arrange(souporcell_log_probability_doublet)
    
    ## Organize Demuxlet num snp
    freemuxlet_snp_df <- result_2 %>% select("barcode", "demuxlet_n_snps")
    freemuxlet_snp_df <- freemuxlet_snp_df %>% arrange(desc(demuxlet_n_snps))
    
    ## Organize Demuxlet num reads
    freemuxlet_reads_df <- result_2 %>% select("barcode", "demuxlet_n_reads")
    freemuxlet_reads_df <- freemuxlet_reads_df %>% arrange(desc(demuxlet_n_reads))
    
    ## Organize Demuxalot probability
    demuxalot_doublet_df <- result_2 %>% select("barcode", "demuxalot_doublet_probability")
    demuxalot_doublet_df <- demuxalot_doublet_df %>% arrange(desc(demuxalot_doublet_probability))
    
    ### Organize metrics by percentile
    ## Vireo diff 
    vireo_doublet_df <- vireo_doublet_df %>% 
        mutate(percentile  = percent_rank(vireo_doublet_logLikRatio))
    
    ## Vireo doublet probability
    vireo_doublet_2_df <- vireo_doublet_2_df %>% 
        mutate(percentile  = percent_rank(vireo_doublet_probability))
    
    ## Demuxlet
    freemuxlet_doublet_df <- freemuxlet_doublet_df %>% 
        mutate(percentile  = percent_rank(demuxlet_DIFF_LLK_SNG_DBL))
    freemuxlet_doublet_df$percentile <- 1 - freemuxlet_doublet_df$percentile
    
    ## Demuxlet num reads
    freemuxlet_reads_df <- freemuxlet_reads_df %>% 
        mutate(percentile  = percent_rank(demuxlet_n_reads))
    freemuxlet_reads_df$percentile <-  freemuxlet_reads_df$percentile
    
    ## Demuxlet num snps
    freemuxlet_snp_df <- freemuxlet_snp_df %>% 
        mutate(percentile  = percent_rank(demuxlet_n_snps))
    freemuxlet_snp_df$percentile <-  freemuxlet_snp_df$percentile
    
    ## Souporcell doublet prob
    souporcell_doublet_df <- souporcell_doublet_df %>% 
        mutate(percentile  = percent_rank(souporcell_log_probability_doublet))
    souporcell_doublet_df$percentile <- 1 - souporcell_doublet_df$percentile
    
    ## Demuxalot doublet prob
    demuxalot_doublet_df <- demuxalot_doublet_df %>% 
        mutate(percentile  = percent_rank(demuxalot_doublet_probability))
    demuxalot_doublet_df$percentile <- demuxalot_doublet_df$percentile
    
    ################################################################################################################################
    #### Compute graph-based doublet detection with optimal pT and nCD values 
    ################################################################################################################################
    ## Report optimal parameters
    message(paste0("Using ", round(optimal_pT, digits = 4), " as optimal pT value"))
    message(paste0("Using ", optimal_nCD, " as optimal nCD value"))
    
    ### Identify high confidence doublets
    suspected_doublets <- list()
    s1 <- rev(seq(0, 1, by = 0.01))
    for (j in s1){
        vireo_doublet_percentile <- subset(vireo_doublet_df, percentile >= j)
        freemuxlet_doublet_percentile <- subset(freemuxlet_doublet_df, percentile >= j)
        souporcell_doublet_percentile <- subset(souporcell_doublet_df, percentile >= j)
        freemuxlet_snp_percentile <- subset(freemuxlet_snp_df, percentile >= j)
        freemuxlet_reads_percentile <- subset(freemuxlet_reads_df, percentile >= j)
        demuxalot_doublet_percentile <- subset(demuxalot_doublet_df, percentile >= j)
        
        for (i in result_2$barcode){
        if (i %in% suspected_doublets){
            print("Barcode already doublet")
        } else if (i %in% vireo_doublet_percentile$barcode &
                    i %in% freemuxlet_doublet_percentile$barcode &
                    i %in% souporcell_doublet_percentile$barcode &
                    i %in% freemuxlet_snp_percentile$barcode &
                    i %in% freemuxlet_reads_percentile$barcode &
                    i %in% demuxalot_doublet_percentile$barcode) {
            print(paste0("adding barcode to doublet list at ", j, "percentile"))
            suspected_doublets <- c(suspected_doublets, i)
        } else {
            print(paste0("barcode not a suspected doublet at", j, "percentile"))
        }
        suspected_doublets
        }
        if(length(suspected_doublets) >= optimal_nCD ) { 
        print(paste0("reached expected doublet rate at", j, "percentile"))
        break
        }
        suspected_doublets
    }  
    
    ### Identify graph-based suspected doublets
    doublet <- "none"
    percentile <- 0
    barcode <- "none"
    filler_frame <- data.frame(doublet,percentile, barcode )
    
    for (j in suspected_doublets){
        distances_test_1 <-distances
        distances_test_1 <- data.frame(distances_test_1[,j])
        colnames(distances_test_1) <- "doublet"
        
        distances_test_1 <- distances_test_1 %>% 
        mutate(percentile  = percent_rank(doublet))
        distances_test_1$percentile <- 1 - distances_test_1$percentile
        
        distances_test_15 <- distances_test_1[distances_test_1$percentile >= optimal_pT,]
        distances_test_15$barcode <- rownames(distances_test_15)
        distances_test_15$doublet <- "temp"
        
        filler_frame <- rbind(filler_frame,distances_test_15)
        filler_frame 
    }
    
    ## Nearest neighbour frequency of GBD-identified doublets
    total_testerquester_count <- filler_frame %>%
        dplyr::group_by(barcode) %>%
        dplyr::mutate(n = n()) %>%
        ungroup()
    total_testerquester_count <- total_testerquester_count[!(duplicated(total_testerquester_count$barcode)),]
    total_testerquester_count <- total_testerquester_count[order(total_testerquester_count$n, decreasing = TRUE),] 
    total_testerquester_count_lim <- total_testerquester_count[c(1:nrow(total_testerquester_count)),]
    remove_dublets <- subset(total_testerquester_count_lim, barcode != "none")
    
    ## Plot kurtosis and save plot
    den <- density(total_testerquester_count_lim$n)
    k <- kurtosis(total_testerquester_count_lim$n)
    k <- round(k, digits = 3)
 
    ### Plot PCA summary
    ## PCA prior to doublet detection
    result_2_lim <- result_2
    result_2_lim$is_doublet <- "singlet"
    result_2_lim$is_doublet[result_2_lim$ensemblex_assignment == "doublet"] <- "doublet"
    result_2_lim <- result_2_lim[,c(1, 14, 19, 20, 21, 17, 23, 25, 13, 31)]
    result_2_lim <- result_2_lim[complete.cases(result_2_lim), ]
    res.pca <- prcomp(result_2_lim[,c(2:8)],scale = T)
    
    fviz_pca_ind(res.pca,
                            col.ind = result_2_lim$is_doublet, 
                            geom="point", 
                            pointsize = 0.5
    ) + ggtitle("Assignment prior to graph-based doublet detection") +
        scale_colour_manual(values = c( "indianred2", "grey"))
        ggsave(paste(par_output_dir,"/step2","/PCA1_graph_based_doublet_detection.pdf", sep=""))

    ## PCA highlighting confident doublets 
    result_2_lim <- result_2
    result_2_lim$is_confident_doublet <- "no"
    result_2_lim$is_confident_doublet[result_2_lim$barcode %in% suspected_doublets ] <- "yes"
    colnames(result_2_lim)
    result_2_lim <- result_2_lim[,c(1, 14, 19, 20, 21, 17, 23, 25, 13, 31)] 
    result_2_lim <- result_2_lim[complete.cases(result_2_lim), ]
    res.pca <- prcomp(result_2_lim[,c(2:8)],scale = T)
    
    fviz_pca_ind(res.pca,
                            col.ind = result_2_lim$is_confident_doublet, 
                            geom="point", 
                            pointsize = 0.5
    ) + ggtitle(paste0("High confidence doublets (nCD = ", optimal_nCD,")")) +
        scale_colour_manual(values = c("grey", "indianred2"))
        ggsave(paste(par_output_dir,"/step2","/PCA2_graph_based_doublet_detection.pdf", sep=""))

    ## PCA after graph baed doublet detection 
    result_2_lim <- result_2
    result_2_lim$is_confident_doublet <- "no"
    result_2_lim$is_confident_doublet[result_2_lim$barcode %in% total_testerquester_count$barcode |  result_2_lim$ensemblex_assignment == "doublet"] <- "yes"
    result_2_lim <- result_2_lim[,c(1, 14, 19, 20, 21, 17, 23, 25, 13, 31)] 
    result_2_lim <- result_2_lim[complete.cases(result_2_lim), ]
    res.pca <- prcomp(result_2_lim[,c(2:8)],scale = T)
    
    fviz_pca_ind(res.pca,
                            col.ind = result_2_lim$is_confident_doublet, 
                            geom="point", 
                            pointsize = 0.5
    ) + ggtitle("Assignment after graph-based doublet detection") +
        scale_colour_manual(values = c("grey", "indianred2"))
    ggsave(paste(par_output_dir,"/step2","/PCA3_graph_based_doublet_detection.pdf", sep=""))
    
    ## Label graph-based expected doublets as doublets
    result_2 <- result_test
    result_2$ensemblex_assignment[result_2$barcode %in% remove_dublets$barcode] <- "doublet"
    write.csv(result_2, paste(par_output_dir,"/step2",'/Step2_cell_assignment.csv', sep=""))
    
    result_2     
}
graph_based_doublet_detection_estimated_par <- function(result_test, par_expected_doublet_rate, par_output_dir){
    
    ## Set seed
    set.seed(1234)
    
    ## create an output directory
    dir.create(paste(par_output_dir,"/step2",sep=''))
    
    ## load Balanced-accuracy dataset
    result_2 <- result_test
    
    ### Perform principal component analysis with select variables
    result_2_lim <- result_2
    result_2_lim <- result_2_lim[,c(1, 14, 19, 20, 21, 17, 23, 25, 13 )] #barcode, vireo_doublet_probability,  souporcell_log_probability_doublet, demuxlet_n_snps, demuxlet_n_reads, vireo_doublet_logLikRatio, demuxlet_DIFF_LLK_SNG_DBL, demuxalot_doublet_probability, vireo_singlet_probability (we dont use this for PCA)
    result_2_lim <- result_2_lim[complete.cases(result_2_lim), ]
    res.pca <- prcomp(result_2_lim[,c(2:8)],scale = T)
    
    ## scree plot
    fviz_eig(res.pca)
    ggsave(paste(par_output_dir,"/step2","/PCA_scree_plot.pdf", sep=""))
    
    ## PCA
    fviz_pca_ind(res.pca,
                col.ind = "black", 
                geom="point", 
                pointsize = 0.5
    )
    ggsave(paste(par_output_dir,"/step2","/PCA_plot.pdf", sep=""))
    
    ### variable contribution to variation
    ## Plot contributions of variables to PC1
    fviz_contrib(res.pca, choice = "var", axes = 1, top = 10) + theme_classic() + coord_flip() + xlab("Variable") + ggtitle("PC1")
    ggsave(paste(par_output_dir,"/step2","/PC1_var_contrib.pdf", sep=""))

    ## Plot contributions of variables to PC2
    fviz_contrib(res.pca, choice = "var", axes = 2, top = 10) + theme_classic() + coord_flip() + xlab("Variable") + ggtitle("PC2")
    ggsave(paste(par_output_dir,"/step2","/PC2_var_contrib.pdf", sep=""))
        
    ## Compute euclidean distance between points
    rownames(result_2_lim) <- result_2_lim$barcode
    res.pca <- prcomp(result_2_lim[,c(2,3,4,5,6,7,8)],scale = T) 
    df_1 <- data.frame(res.pca$x[,1])
    df_1$barcode <- rownames(df_1)
    df_2 <- data.frame(res.pca$x[,2])
    df_2$barcode <- rownames(df_2)
    df_merge_test <- merge(df_1, df_2, by = "barcode")
    colnames(df_merge_test) <- c("barcode", "PC1", "PC2")
    distances <- dist(df_merge_test[c("PC1", "PC2")], diag = TRUE, upper = TRUE)
    distances <- as.matrix(distances)
    colnames(distances) <- df_merge_test$barcode
    rownames(distances) <- df_merge_test$barcode
    
    ### Organize parameters to identify most likely doublets
    ## Organize Vireo doublet log_lik into ordered frame
    vireo_doublet_df <- result_2 %>% select("barcode", "vireo_doublet_logLikRatio")
    vireo_doublet_df <- vireo_doublet_df %>% arrange(desc(vireo_doublet_logLikRatio))
    
    ## Organize Vireo doublet probability
    vireo_doublet_2_df <- result_2 %>% select("barcode", "vireo_doublet_probability")
    vireo_doublet_2_df <- vireo_doublet_2_df %>% arrange(desc(vireo_doublet_probability))
    
    ## Organize Demuxlet "freemuxlet_DIFF_LLK_SNG_DBL" in ordered frame
    freemuxlet_doublet_df <- result_2 %>% select("barcode", "demuxlet_DIFF_LLK_SNG_DBL")
    freemuxlet_doublet_df <- freemuxlet_doublet_df %>% arrange(demuxlet_DIFF_LLK_SNG_DBL)
    
    ## Organize Souporcell log prob doublet 
    souporcell_doublet_df <- result_2 %>% select("barcode", "souporcell_log_probability_doublet")
    souporcell_doublet_df <- souporcell_doublet_df %>% arrange(souporcell_log_probability_doublet)
    
    ## Organize Demuxlet num snp
    freemuxlet_snp_df <- result_2 %>% select("barcode", "demuxlet_n_snps")
    freemuxlet_snp_df <- freemuxlet_snp_df %>% arrange(desc(demuxlet_n_snps))
    
    ## Organize Demuxlet num reads
    freemuxlet_reads_df <- result_2 %>% select("barcode", "demuxlet_n_reads")
    freemuxlet_reads_df <- freemuxlet_reads_df %>% arrange(desc(demuxlet_n_reads))
    
    ## Organize Demuxalot probability
    demuxalot_doublet_df <- result_2 %>% select("barcode", "demuxalot_doublet_probability")
    demuxalot_doublet_df <- demuxalot_doublet_df %>% arrange(desc(demuxalot_doublet_probability))
    
    ### Organize metrics by percentile
    ## Vireo diff 
    vireo_doublet_df <- vireo_doublet_df %>% 
        mutate(percentile  = percent_rank(vireo_doublet_logLikRatio))
    
    ## Vireo doublet probability
    vireo_doublet_2_df <- vireo_doublet_2_df %>% 
        mutate(percentile  = percent_rank(vireo_doublet_probability))
    
    ## Demuxlet
    freemuxlet_doublet_df <- freemuxlet_doublet_df %>% 
        mutate(percentile  = percent_rank(demuxlet_DIFF_LLK_SNG_DBL))
    freemuxlet_doublet_df$percentile <- 1 - freemuxlet_doublet_df$percentile
    
    ## Demuxlet num reads
    freemuxlet_reads_df <- freemuxlet_reads_df %>% 
        mutate(percentile  = percent_rank(demuxlet_n_reads))
    freemuxlet_reads_df$percentile <-  freemuxlet_reads_df$percentile
    
    ## Demuxlet num snps
    freemuxlet_snp_df <- freemuxlet_snp_df %>% 
        mutate(percentile  = percent_rank(demuxlet_n_snps))
    freemuxlet_snp_df$percentile <-  freemuxlet_snp_df$percentile
    
    ## Souporcell doublet prob
    souporcell_doublet_df <- souporcell_doublet_df %>% 
        mutate(percentile  = percent_rank(souporcell_log_probability_doublet))
    souporcell_doublet_df$percentile <- 1 - souporcell_doublet_df$percentile
    
    ## Demuxalot doublet prob
    demuxalot_doublet_df <- demuxalot_doublet_df %>% 
        mutate(percentile  = percent_rank(demuxalot_doublet_probability))
    demuxalot_doublet_df$percentile <- demuxalot_doublet_df$percentile
    
    ################################################################################################################################
    #### parameter sweep #### 
    ################################################################################################################################
    ### Compute varying number of confident doublets (nCD)
    ## 50 nCD 
    suspected_doublets_50 <- list()
    s1 <- rev(seq(0, 1, by = 0.01))
    for (j in s1){
        vireo_doublet_percentile <- subset(vireo_doublet_df, percentile >= j)
        vireo_doublet_2_percentile <- subset(vireo_doublet_2_df, percentile >= j)
        freemuxlet_doublet_percentile <- subset(freemuxlet_doublet_df, percentile >= j)
        souporcell_doublet_percentile <- subset(souporcell_doublet_df, percentile >= j)
        freemuxlet_snp_percentile <- subset(freemuxlet_snp_df, percentile >= j)
        freemuxlet_reads_percentile <- subset(freemuxlet_reads_df, percentile >= j)
        demuxalot_doublet_percentile <- subset(demuxalot_doublet_df, percentile >= j)
        
        for (i in result_2$barcode){
        if (i %in% suspected_doublets_50){
            print("Barcode already doublet")
        } else if (i %in% vireo_doublet_percentile$barcode &
                    i %in% vireo_doublet_2_percentile$barcode &
                    i %in% freemuxlet_doublet_percentile$barcode &
                    i %in% souporcell_doublet_percentile$barcode &
                    i %in% freemuxlet_snp_percentile$barcode &
                    i %in% freemuxlet_reads_percentile$barcode &
                    i %in% demuxalot_doublet_percentile$barcode) {
            print(paste0("adding barcode to doublet list at ", j, "percentile"))
            suspected_doublets_50 <- c(suspected_doublets_50, i)
        } else {
            print(paste0("barcode not a suspected doublet at", j, "percentile"))
        }
        suspected_doublets_50
        }
        if(length(suspected_doublets_50) >= 50 ) { 
        print(paste0("reached expected doublet rate at", j, "percentile"))
        break
        }
        suspected_doublets_50
    }  
    
    ## 100 nCD
    suspected_doublets_100 <- list()
    s1 <- rev(seq(0, 1, by = 0.01))
    for (j in s1){
        vireo_doublet_percentile <- subset(vireo_doublet_df, percentile >= j)
        vireo_doublet_2_percentile <- subset(vireo_doublet_2_df, percentile >= j)
        freemuxlet_doublet_percentile <- subset(freemuxlet_doublet_df, percentile >= j)
        souporcell_doublet_percentile <- subset(souporcell_doublet_df, percentile >= j)
        freemuxlet_snp_percentile <- subset(freemuxlet_snp_df, percentile >= j)
        freemuxlet_reads_percentile <- subset(freemuxlet_reads_df, percentile >= j)
        demuxalot_doublet_percentile <- subset(demuxalot_doublet_df, percentile >= j)
        
        for (i in result_2$barcode){
        if (i %in% suspected_doublets_100){
            print("Barcode already doublet")
        } else if (i %in% vireo_doublet_percentile$barcode &
                    i %in% vireo_doublet_2_percentile$barcode &
                    i %in% freemuxlet_doublet_percentile$barcode &
                    i %in% souporcell_doublet_percentile$barcode &
                    i %in% freemuxlet_snp_percentile$barcode &
                    i %in% freemuxlet_reads_percentile$barcode &
                    i %in% demuxalot_doublet_percentile$barcode) {
            print(paste0("adding barcode to doublet list at ", j, "percentile"))
            suspected_doublets_100 <- c(suspected_doublets_100, i)
        } else {
            print(paste0("barcode not a suspected doublet at", j, "percentile"))
        }
        suspected_doublets_100
        }
        if(length(suspected_doublets_100) >= 100 ) { 
        print(paste0("reached expected doublet rate at", j, "percentile"))
        break
        }
        suspected_doublets_100
    }
    
    ## 150 nCD
    suspected_doublets_150 <- list()
    s1 <- rev(seq(0, 1, by = 0.01))
    for (j in s1){
        vireo_doublet_percentile <- subset(vireo_doublet_df, percentile >= j)
        vireo_doublet_2_percentile <- subset(vireo_doublet_2_df, percentile >= j)
        freemuxlet_doublet_percentile <- subset(freemuxlet_doublet_df, percentile >= j)
        souporcell_doublet_percentile <- subset(souporcell_doublet_df, percentile >= j)
        freemuxlet_snp_percentile <- subset(freemuxlet_snp_df, percentile >= j)
        freemuxlet_reads_percentile <- subset(freemuxlet_reads_df, percentile >= j)
        demuxalot_doublet_percentile <- subset(demuxalot_doublet_df, percentile >= j)
        
        for (i in result_2$barcode){
        if (i %in% suspected_doublets_150){
            print("Barcode already doublet")
        } else if (i %in% vireo_doublet_percentile$barcode &
                    i %in% vireo_doublet_2_percentile$barcode &
                    i %in% freemuxlet_doublet_percentile$barcode &
                    i %in% souporcell_doublet_percentile$barcode &
                    i %in% freemuxlet_snp_percentile$barcode &
                    i %in% freemuxlet_reads_percentile$barcode &
                    i %in% demuxalot_doublet_percentile$barcode) {
            print(paste0("adding barcode to doublet list at ", j, "percentile"))
            suspected_doublets_150 <- c(suspected_doublets_150, i)
        } else {
            print(paste0("barcode not a suspected doublet at", j, "percentile"))
        }
        suspected_doublets_150
        }
        if(length(suspected_doublets_150) >= 150 ) { #OG = 100
        print(paste0("reached expected doublet rate at", j, "percentile"))
        break
        }
        suspected_doublets_150
    }
    
    ## 200 nCD
    suspected_doublets_200 <- list()
    s1 <- rev(seq(0, 1, by = 0.01))
    for (j in s1){
        vireo_doublet_percentile <- subset(vireo_doublet_df, percentile >= j)
        vireo_doublet_2_percentile <- subset(vireo_doublet_2_df, percentile >= j)
        freemuxlet_doublet_percentile <- subset(freemuxlet_doublet_df, percentile >= j)
        souporcell_doublet_percentile <- subset(souporcell_doublet_df, percentile >= j)
        freemuxlet_snp_percentile <- subset(freemuxlet_snp_df, percentile >= j)
        freemuxlet_reads_percentile <- subset(freemuxlet_reads_df, percentile >= j)
        demuxalot_doublet_percentile <- subset(demuxalot_doublet_df, percentile >= j)
        
        for (i in result_2$barcode){
        if (i %in% suspected_doublets_200){
            print("Barcode already doublet")
        } else if (i %in% vireo_doublet_percentile$barcode &
                    i %in% vireo_doublet_2_percentile$barcode &
                    i %in% freemuxlet_doublet_percentile$barcode &
                    i %in% souporcell_doublet_percentile$barcode &
                    i %in% freemuxlet_snp_percentile$barcode &
                    i %in% freemuxlet_reads_percentile$barcode &
                    i %in% demuxalot_doublet_percentile$barcode) {
            print(paste0("adding barcode to doublet list at ", j, "percentile"))
            suspected_doublets_200 <- c(suspected_doublets_200, i)
        } else {
            print(paste0("barcode not a suspected doublet at", j, "percentile"))
        }
        suspected_doublets_200
        }
        if(length(suspected_doublets_200) >= 200 ) { #OG = 100
        print(paste0("reached expected doublet rate at", j, "percentile"))
        break
        }
        suspected_doublets_200
    }
    
    ## 250 nCD 
    suspected_doublets_250 <- list()
    s1 <- rev(seq(0, 1, by = 0.01))
    for (j in s1){
        vireo_doublet_percentile <- subset(vireo_doublet_df, percentile >= j)
        vireo_doublet_2_percentile <- subset(vireo_doublet_2_df, percentile >= j)
        freemuxlet_doublet_percentile <- subset(freemuxlet_doublet_df, percentile >= j)
        souporcell_doublet_percentile <- subset(souporcell_doublet_df, percentile >= j)
        freemuxlet_snp_percentile <- subset(freemuxlet_snp_df, percentile >= j)
        freemuxlet_reads_percentile <- subset(freemuxlet_reads_df, percentile >= j)
        demuxalot_doublet_percentile <- subset(demuxalot_doublet_df, percentile >= j)
        
        for (i in result_2$barcode){
        if (i %in% suspected_doublets_250){
            print("Barcode already doublet")
        } else if (i %in% vireo_doublet_percentile$barcode &
                    i %in% vireo_doublet_2_percentile$barcode &
                    i %in% freemuxlet_doublet_percentile$barcode &
                    i %in% souporcell_doublet_percentile$barcode &
                    i %in% freemuxlet_snp_percentile$barcode &
                    i %in% freemuxlet_reads_percentile$barcode &
                    i %in% demuxalot_doublet_percentile$barcode) {
            print(paste0("adding barcode to doublet list at ", j, "percentile"))
            suspected_doublets_250 <- c(suspected_doublets_250, i)
        } else {
            print(paste0("barcode not a suspected doublet at", j, "percentile"))
        }
        suspected_doublets_250
        }
        if(length(suspected_doublets_250) >= 250 ) { #OG = 100
        print(paste0("reached expected doublet rate at", j, "percentile"))
        break
        }
        suspected_doublets_250
    }
    
    ## 300 nCD
    suspected_doublets_300 <- list()
    s1 <- rev(seq(0, 1, by = 0.01))
    for (j in s1){
        vireo_doublet_percentile <- subset(vireo_doublet_df, percentile >= j)
        vireo_doublet_2_percentile <- subset(vireo_doublet_2_df, percentile >= j)
        freemuxlet_doublet_percentile <- subset(freemuxlet_doublet_df, percentile >= j)
        souporcell_doublet_percentile <- subset(souporcell_doublet_df, percentile >= j)
        freemuxlet_snp_percentile <- subset(freemuxlet_snp_df, percentile >= j)
        freemuxlet_reads_percentile <- subset(freemuxlet_reads_df, percentile >= j)
        demuxalot_doublet_percentile <- subset(demuxalot_doublet_df, percentile >= j)
        
        for (i in result_2$barcode){
        if (i %in% suspected_doublets_300){
            print("Barcode already doublet")
        } else if (i %in% vireo_doublet_percentile$barcode &
                    i %in% vireo_doublet_2_percentile$barcode &
                    i %in% freemuxlet_doublet_percentile$barcode &
                    i %in% souporcell_doublet_percentile$barcode &
                    i %in% freemuxlet_snp_percentile$barcode &
                    i %in% freemuxlet_reads_percentile$barcode &
                    i %in% demuxalot_doublet_percentile$barcode) {
            print(paste0("adding barcode to doublet list at ", j, "percentile"))
            suspected_doublets_300 <- c(suspected_doublets_300, i)
        } else {
            print(paste0("barcode not a suspected doublet at", j, "percentile"))
        }
        suspected_doublets_300
        }
        if(length(suspected_doublets_300) >= 300 ) { #OG = 100
        print(paste0("reached expected doublet rate at", j, "percentile"))
        break
        }
        suspected_doublets_300
    }
    

    ## Make a list of percentile threshold (pT) values based on expected doublet rate for the pool
    message(paste0("Expected doublet rate is ", par_expected_doublet_rate, "; ", round(par_expected_doublet_rate*nrow(result_2), digits = 0), " droplets"  ))
    interval <- par_expected_doublet_rate/6
    pT_list <- rev(seq(1-par_expected_doublet_rate,1-interval, by = interval))
        
    ## 50 nCD pT sweep 
    doublet <- "none"
    percentile <- 0
    barcode <- "none"
    pT <- 0
    filler_frame_50 <- data.frame(doublet,percentile, barcode,pT )

    for (t in pT_list){
        for (j in suspected_doublets_50){
        distances_test_1 <-distances
        distances_test_1 <- data.frame(distances_test_1[,j])
        colnames(distances_test_1) <- "doublet"

        distances_test_1 <- distances_test_1 %>% 
        mutate(percentile  = percent_rank(doublet))
        distances_test_1$percentile <- 1 - distances_test_1$percentile
        
        distances_test_15 <- distances_test_1[distances_test_1$percentile >= t,] 
        distances_test_15$barcode <- rownames(distances_test_15)
        distances_test_15$doublet <- "temp"
        distances_test_15$pT <- t
        
        filler_frame_50 <- rbind(filler_frame_50,distances_test_15)
        filler_frame_50  
        }
        filler_frame_50
    }
    
    ## 100 nCD pT sweep 
    doublet <- "none"
    percentile <- 0
    barcode <- "none"
    pT <- 0
    filler_frame_100 <- data.frame(doublet,percentile, barcode,pT )
    
    for (t in pT_list){
        for (j in suspected_doublets_100){
        distances_test_1 <-distances
        distances_test_1 <- data.frame(distances_test_1[,j])
        colnames(distances_test_1) <- "doublet"
        
        distances_test_1 <- distances_test_1 %>% 
        mutate(percentile  = percent_rank(doublet))
        distances_test_1$percentile <- 1 - distances_test_1$percentile
        
        distances_test_15 <- distances_test_1[distances_test_1$percentile >= t,] 
        distances_test_15$barcode <- rownames(distances_test_15)
        distances_test_15$doublet <- "temp"
        distances_test_15$pT <- t
        
        filler_frame_100 <- rbind(filler_frame_100,distances_test_15)
        filler_frame_100 
        }
        filler_frame_100
    }
    
    ## 150 nCD pT sweep
    doublet <- "none"
    percentile <- 0
    barcode <- "none"
    pT <- 0
    filler_frame_150 <- data.frame(doublet,percentile, barcode,pT )
    
    for (t in pT_list){
        for (j in suspected_doublets_150){
        distances_test_1 <-distances
        distances_test_1 <- data.frame(distances_test_1[,j])
        colnames(distances_test_1) <- "doublet"
        
        distances_test_1 <- distances_test_1 %>% 
        mutate(percentile  = percent_rank(doublet))
        distances_test_1$percentile <- 1 - distances_test_1$percentile
        
        distances_test_15 <- distances_test_1[distances_test_1$percentile >= t,]
        distances_test_15$barcode <- rownames(distances_test_15)
        distances_test_15$doublet <- "temp"
        distances_test_15$pT <- t
        
        filler_frame_150 <- rbind(filler_frame_150,distances_test_15)
        filler_frame_150  
        }
        filler_frame_150
    }
    
    ## 200 nCD pT sweep 
    doublet <- "none"
    percentile <- 0
    barcode <- "none"
    pT <- 0
    filler_frame_200 <- data.frame(doublet,percentile, barcode,pT )
    
    for (t in pT_list){
        for (j in suspected_doublets_200){
        distances_test_1 <-distances
        distances_test_1 <- data.frame(distances_test_1[,j])
        colnames(distances_test_1) <- "doublet"
        
        distances_test_1 <- distances_test_1 %>% 
        mutate(percentile  = percent_rank(doublet))
        distances_test_1$percentile <- 1 - distances_test_1$percentile
        
        distances_test_15 <- distances_test_1[distances_test_1$percentile >= t,] 
        distances_test_15$barcode <- rownames(distances_test_15)
        distances_test_15$doublet <- "temp"
        distances_test_15$pT <- t
        
        filler_frame_200 <- rbind(filler_frame_200,distances_test_15)
        filler_frame_200
        }
        filler_frame_200
    }
    
    ## 250 nCD pT sweep 
    doublet <- "none"
    percentile <- 0
    barcode <- "none"
    pT <- 0
    filler_frame_250 <- data.frame(doublet,percentile, barcode,pT )
    
    for (t in pT_list){
        for (j in suspected_doublets_250){
        distances_test_1 <-distances
        distances_test_1 <- data.frame(distances_test_1[,j])
        colnames(distances_test_1) <- "doublet"
        
        distances_test_1 <- distances_test_1 %>% 
        mutate(percentile  = percent_rank(doublet))
        distances_test_1$percentile <- 1 - distances_test_1$percentile
        
        distances_test_15 <- distances_test_1[distances_test_1$percentile >= t,]
        distances_test_15$barcode <- rownames(distances_test_15)
        distances_test_15$doublet <- "temp"
        distances_test_15$pT <- t
        
        filler_frame_250 <- rbind(filler_frame_250,distances_test_15)
        filler_frame_250
        }
        filler_frame_250
    }
    
    ## 300 nCD pT sweep 
    doublet <- "none"
    percentile <- 0
    barcode <- "none"
    pT <- 0
    filler_frame_300 <- data.frame(doublet,percentile, barcode,pT )
    
    for (t in pT_list){
        for (j in suspected_doublets_300){
        distances_test_1 <-distances
        distances_test_1 <- data.frame(distances_test_1[,j])
        colnames(distances_test_1) <- "doublet"
        
        distances_test_1 <- distances_test_1 %>% 
        mutate(percentile  = percent_rank(doublet))
        distances_test_1$percentile <- 1 - distances_test_1$percentile
        
        distances_test_15 <- distances_test_1[distances_test_1$percentile >= t,] 
        distances_test_15$barcode <- rownames(distances_test_15)
        distances_test_15$doublet <- "temp"
        distances_test_15$pT <- t
        
        filler_frame_300 <- rbind(filler_frame_300,distances_test_15)
        filler_frame_300
        }
        filler_frame_300
    }
    
    ### Compute nearest neighbour frequency and Kutosis of frequency distributions
    ## 50 nCD 
    counts_50 <- filler_frame_50 %>%
        group_by(barcode,pT) %>%
        dplyr::mutate(n = n()) %>%
        ungroup()
    
    counts_50 <- counts_50 %>% 
        distinct(pT, barcode, .keep_all = TRUE)
    
    ncD <- 50
    pT <- 0
    kurtosis <- 0
    fill_frame_50 <- data.frame(ncD,pT, kurtosis)
    
    for (t in pT_list){
        k <- kurtosis(counts_50$n[counts_50$pT == t])
        nCD <- 50
        pT <- t
        kurtosis <- k
        temp_frame <- data.frame(ncD,pT, kurtosis)
        fill_frame_50 <- rbind(fill_frame_50,temp_frame)
        fill_frame_50 <- subset(fill_frame_50, pT != 0)
        fill_frame_50
    }
    
    ## 100 nCD  
    counts_100 <- filler_frame_100 %>%
        group_by(barcode,pT) %>%
        dplyr::mutate(n = n()) %>%
        ungroup()
    
    counts_100 <- counts_100 %>% 
        distinct(pT, barcode, .keep_all = TRUE)
    
    ncD <- 100
    pT <- 0
    kurtosis <- 0
    fill_frame_100 <- data.frame(ncD,pT, kurtosis)
    
    for (t in pT_list){
        k <- kurtosis(counts_100$n[counts_100$pT == t])
        nCD <- 100
        pT <- t
        kurtosis <- k
        temp_frame <- data.frame(ncD,pT, kurtosis)
        fill_frame_100 <- rbind(fill_frame_100,temp_frame)
        fill_frame_100 <- subset(fill_frame_100, pT != 0)
        fill_frame_100
    }
    
    ## 150 nCD 
    counts_150 <- filler_frame_150 %>%
        group_by(barcode,pT) %>%
        dplyr::mutate(n = n()) %>%
        ungroup()
    
    counts_150 <- counts_150 %>% 
        distinct(pT, barcode, .keep_all = TRUE)
    
    ncD <- 150
    pT <- 0
    kurtosis <- 0
    fill_frame_150 <- data.frame(ncD,pT, kurtosis)
    
    for (t in pT_list){
        k <- kurtosis(counts_150$n[counts_150$pT == t])
        nCD <- 150
        pT <- t
        kurtosis <- k
        temp_frame <- data.frame(ncD,pT, kurtosis)
        fill_frame_150 <- rbind(fill_frame_150,temp_frame)
        fill_frame_150 <- subset(fill_frame_150, pT != 0)
        fill_frame_150
    }
    
    ## 200 nCD  
    counts_200 <- filler_frame_200 %>%
        group_by(barcode,pT) %>%
        dplyr::mutate(n = n()) %>%
        ungroup()
    
    counts_200 <- counts_200 %>% 
        distinct(pT, barcode, .keep_all = TRUE)
    
    ncD <- 200
    pT <- 0
    kurtosis <- 0
    fill_frame_200 <- data.frame(ncD,pT, kurtosis)
    
    for (t in pT_list){
        k <- kurtosis(counts_200$n[counts_200$pT == t])
        nCD <- 200
        pT <- t
        kurtosis <- k
        temp_frame <- data.frame(ncD,pT, kurtosis)
        fill_frame_200 <- rbind(fill_frame_200,temp_frame)
        fill_frame_200 <- subset(fill_frame_200, pT != 0)
        fill_frame_200
    }
    
    ## 250 nCD 
    counts_250 <- filler_frame_250 %>%
        group_by(barcode,pT) %>%
        dplyr::mutate(n = n()) %>%
        ungroup()

    counts_250 <- counts_250 %>% 
        distinct(pT, barcode, .keep_all = TRUE)
    
    ncD <- 250
    pT <- 0
    kurtosis <- 0
    fill_frame_250 <- data.frame(ncD,pT, kurtosis)
    
    for (t in pT_list){
        k <- kurtosis(counts_250$n[counts_250$pT == t])
        nCD <- 250
        pT <- t
        kurtosis <- k
        temp_frame <- data.frame(ncD,pT, kurtosis)
        fill_frame_250 <- rbind(fill_frame_250,temp_frame)
        fill_frame_250 <- subset(fill_frame_250, pT != 0)
        fill_frame_250
    }
    
    ## 300 nCD 
    counts_300 <- filler_frame_300 %>%
        group_by(barcode,pT) %>%
        dplyr::mutate(n = n()) %>%
        ungroup()

    counts_300 <- counts_300 %>% 
        distinct(pT, barcode, .keep_all = TRUE)
    
    ncD <- 300
    pT <- 0
    kurtosis <- 0
    fill_frame_300 <- data.frame(ncD,pT, kurtosis)
    
    for (t in pT_list){
        k <- kurtosis(counts_300$n[counts_300$pT == t])
        nCD <- 300
        pT <- t
        kurtosis <- k
        temp_frame <- data.frame(ncD,pT, kurtosis)
        fill_frame_300 <- rbind(fill_frame_300,temp_frame)
        fill_frame_300 <- subset(fill_frame_300, pT != 0)
        fill_frame_300
    }
    
    ################################################################################################################################
    ### Compute optimal parameters 
    ################################################################################################################################
    ## Bind all of the Kurtosis frames 
    kurtosis_bind <- rbind(fill_frame_50,fill_frame_100,fill_frame_150,fill_frame_200,fill_frame_250, fill_frame_300)
    
    ## Compute average Kurtosis for each pT
    kurtosis_bind <- kurtosis_bind %>%
        dplyr::group_by(pT) %>%
        dplyr::summarize(Mean_k = mean(kurtosis, na.rm=TRUE))
    
    optimal_pT <- kurtosis_bind$pT[kurtosis_bind$Mean_k == max(kurtosis_bind$Mean_k)] 
    
    ## Plot kurtosis values
    ggplot(kurtosis_bind, aes(x = pT, y = Mean_k)) + 
        geom_vline(xintercept = optimal_pT, col = "red") +
        geom_point() +  
        geom_line() + 
        theme_classic() +
        ggtitle(paste0("Estimated optimal pT: ", optimal_pT))
    ggsave(paste(par_output_dir,"/step2","/optimal_pT.pdf", sep=""))
    
    ## find optimal nCD
    bind_nCD <- rbind(fill_frame_50,fill_frame_100,fill_frame_150,fill_frame_200,fill_frame_250, fill_frame_300)
    bind_nCD <- subset(bind_nCD, pT == optimal_pT)
    ncD <- data.frame(nCD = c(50, 100, 150, 200, 250, 300))
    ncD_list <- c(50, 100, 150, 200, 250, 300)
    
    ## Smooth line
    fm1 <- smooth.spline(bind_nCD[,"ncD"], bind_nCD[,"kurtosis"], df = 3)
    y2 <- predict(fm1, x = ncD)
    y <- data.frame(y2$y)
    y <- y$nCD

    ## Find elbow of the smoothed curve 
    df <- data.frame(ncd = ncD_list, kurtosis = y)
    optimal_nCD <- find_curve_elbow(df, export_type = "row_num", plot_curve = FALSE)

    ## parse value
    optimal_nCD <- df[optimal_nCD, 1]

    ## If cannot find slope, take point preceedinging the largest slope
    if (is.na(optimal_nCD)) {
        message("Could not determine optimal nCD based on the elbow, taking maximum kurtosis value.")
        max <- max(bind_nCD$kurtosis)
        optimal_nCD <- bind_nCD$ncD[bind_nCD$kurtosis == max]
        optimal_nCD <- optimal_nCD[1]
        #slope= diff(bind_nCD$ncD)/diff(bind_nCD$kurtosis)
        #xx <- bind_nCD[which.max(slope),]
        #optimal_nCD <- max(xx$ncD)
    }else{
        optimal_nCD <- optimal_nCD
        message("Determine optimal nCD based on the elbow.")
    }

    ## plot 
    ggplot(bind_nCD, aes(x = ncD, y = kurtosis)) + 
        geom_vline(xintercept = optimal_nCD, col = "red") +
        geom_point() +  
        geom_line() + 
        theme_classic() +
        ggtitle(paste0("Estimated optimal nCD: ", optimal_nCD))
    ggsave(paste(par_output_dir,"/step2","/optimal_nCD.pdf", sep=""))

    ################################################################################################################################
    #### Compute graph-based doublet detection with optimal pT and nCD values 
    ################################################################################################################################
    ## Report optimal parameters
    message(paste0("Using ", round(optimal_pT, digits = 4), " as optimal pT value"))
    message(paste0("Using ", optimal_nCD, " as optimal nCD value"))
    
    ### Identify high confidence doublets
    suspected_doublets <- list()
    s1 <- rev(seq(0, 1, by = 0.01))
    for (j in s1){
        vireo_doublet_percentile <- subset(vireo_doublet_df, percentile >= j)
        freemuxlet_doublet_percentile <- subset(freemuxlet_doublet_df, percentile >= j)
        souporcell_doublet_percentile <- subset(souporcell_doublet_df, percentile >= j)
        freemuxlet_snp_percentile <- subset(freemuxlet_snp_df, percentile >= j)
        freemuxlet_reads_percentile <- subset(freemuxlet_reads_df, percentile >= j)
        demuxalot_doublet_percentile <- subset(demuxalot_doublet_df, percentile >= j)
        
        for (i in result_2$barcode){
        if (i %in% suspected_doublets){
            print("Barcode already doublet")
        } else if (i %in% vireo_doublet_percentile$barcode &
                    i %in% freemuxlet_doublet_percentile$barcode &
                    i %in% souporcell_doublet_percentile$barcode &
                    i %in% freemuxlet_snp_percentile$barcode &
                    i %in% freemuxlet_reads_percentile$barcode &
                    i %in% demuxalot_doublet_percentile$barcode) {
            print(paste0("adding barcode to doublet list at ", j, "percentile"))
            suspected_doublets <- c(suspected_doublets, i)
        } else {
            print(paste0("barcode not a suspected doublet at", j, "percentile"))
        }
        suspected_doublets
        }
        if(length(suspected_doublets) >= optimal_nCD ) { 
        print(paste0("reached expected doublet rate at", j, "percentile"))
        break
        }
        suspected_doublets
    }  
    
    ### Identify graph-based suspected doublets
    doublet <- "none"
    percentile <- 0
    barcode <- "none"
    filler_frame <- data.frame(doublet,percentile, barcode )
    
    for (j in suspected_doublets){
        distances_test_1 <-distances
        distances_test_1 <- data.frame(distances_test_1[,j])
        colnames(distances_test_1) <- "doublet"
        
        distances_test_1 <- distances_test_1 %>% 
        mutate(percentile  = percent_rank(doublet))
        distances_test_1$percentile <- 1 - distances_test_1$percentile
        
        distances_test_15 <- distances_test_1[distances_test_1$percentile >= optimal_pT,]
        distances_test_15$barcode <- rownames(distances_test_15)
        distances_test_15$doublet <- "temp"
        
        filler_frame <- rbind(filler_frame,distances_test_15)
        filler_frame 
    }
    
    ## Nearest neighbour frequency of GBD-identified doublets
    total_testerquester_count <- filler_frame %>%
        dplyr::group_by(barcode) %>%
        dplyr::mutate(n = n()) %>%
        ungroup()
    total_testerquester_count <- total_testerquester_count[!(duplicated(total_testerquester_count$barcode)),]
    total_testerquester_count <- total_testerquester_count[order(total_testerquester_count$n, decreasing = TRUE),] 
    total_testerquester_count_lim <- total_testerquester_count[c(1:nrow(total_testerquester_count)),]
    remove_dublets <- subset(total_testerquester_count_lim, barcode != "none")
    
    ## Plot kurtosis and save plot
    den <- density(total_testerquester_count_lim$n)
    k <- kurtosis(total_testerquester_count_lim$n)
    k <- round(k, digits = 3)
    
    ### Plot PCA summary
    ## PCA prior to doublet detection
    result_2_lim <- result_2
    result_2_lim$is_doublet <- "singlet"
    result_2_lim$is_doublet[result_2_lim$ensemblex_assignment == "doublet"] <- "doublet"
    result_2_lim <- result_2_lim[,c(1, 14, 19, 20, 21, 17, 23, 25, 13, 31)]
    result_2_lim <- result_2_lim[complete.cases(result_2_lim), ]
    res.pca <- prcomp(result_2_lim[,c(2:8)],scale = T)
    
    fviz_pca_ind(res.pca,
                            col.ind = result_2_lim$is_doublet, 
                            geom="point", 
                            pointsize = 0.5
    ) + ggtitle("Assignment prior to graph-based doublet detection") +
        scale_colour_manual(values = c( "indianred2", "grey"))
        ggsave(paste(par_output_dir,"/step2","/PCA1_graph_based_doublet_detection.pdf", sep=""))

    ## PCA highlighting confident doublets 
    result_2_lim <- result_2
    result_2_lim$is_confident_doublet <- "no"
    result_2_lim$is_confident_doublet[result_2_lim$barcode %in% suspected_doublets ] <- "yes"
    colnames(result_2_lim)
    result_2_lim <- result_2_lim[,c(1, 14, 19, 20, 21, 17, 23, 25, 13, 31)] 
    result_2_lim <- result_2_lim[complete.cases(result_2_lim), ]
    res.pca <- prcomp(result_2_lim[,c(2:8)],scale = T)
    
    fviz_pca_ind(res.pca,
                            col.ind = result_2_lim$is_confident_doublet, 
                            geom="point", 
                            pointsize = 0.5
    ) + ggtitle(paste0("High confidence doublets (nCD = ", optimal_nCD,")")) +
        scale_colour_manual(values = c("grey", "indianred2"))
        ggsave(paste(par_output_dir,"/step2","/PCA2_graph_based_doublet_detection.pdf", sep=""))

    ## PCA after graph baed doublet detection 
    result_2_lim <- result_2
    result_2_lim$is_confident_doublet <- "no"
    result_2_lim$is_confident_doublet[result_2_lim$barcode %in% total_testerquester_count$barcode |  result_2_lim$ensemblex_assignment == "doublet"] <- "yes"
    result_2_lim <- result_2_lim[,c(1, 14, 19, 20, 21, 17, 23, 25, 13, 31)] 
    result_2_lim <- result_2_lim[complete.cases(result_2_lim), ]
    res.pca <- prcomp(result_2_lim[,c(2:8)],scale = T)
    
    fviz_pca_ind(res.pca,
                            col.ind = result_2_lim$is_confident_doublet, 
                            geom="point", 
                            pointsize = 0.5
    ) + ggtitle("Assignment after graph-based doublet detection") +
        scale_colour_manual(values = c("grey", "indianred2"))
    ggsave(paste(par_output_dir,"/step2","/PCA3_graph_based_doublet_detection.pdf", sep=""))
    
    ## Label graph-based expected doublets as doublets
    result_2 <- result_test
    result_2$ensemblex_assignment[result_2$barcode %in% remove_dublets$barcode] <- "doublet"
    write.csv(result_2, paste(par_output_dir,"/step2",'/Step2_cell_assignment.csv', sep=""))
    
    result_2     
}

###########################################################################################################################
# Ensemble-independent doublet detection
###########################################################################################################################
ensemble_independent_doublet_detections <- function(result_2, par_output_dir){
  ## Set seed
  set.seed(1234)
  
  ## Create an output directory
  dir.create(paste(par_output_dir,"/step3",sep=''))

  expected_doublets <- nrow(result_2)*par_expected_doublet_rate

  ### Proportion agreement bar plot with threshold
  ## Demuxalot
  result_2_demuxalot <- result_2[result_2$demuxalot_assignment == "doublet",]
  result_2_demuxalot <- result_2_demuxalot %>% dplyr::select(vireo_best_assignment, demuxlet_best_assignment, souporcell_best_assignment)
  result_2_demuxalot$one <- 0
  result_2_demuxalot$one[result_2_demuxalot$vireo_best_assignment == "doublet"] <- 1
  result_2_demuxalot$two <- 0
  result_2_demuxalot$two[result_2_demuxalot$demuxlet_best_assignment == "doublet"] <- 1
  result_2_demuxalot$three <- 0
  result_2_demuxalot$three[result_2_demuxalot$souporcell_best_assignment == "doublet"] <- 1
  result_2_demuxalot <- result_2_demuxalot %>% dplyr::select(one, two, three)
  result_2_demuxalot$sum <- rowSums(result_2_demuxalot)
  result_2_demuxalot$tool <- "Demuxalot"
  ## Vireo
  result_2_vireo <- result_2[result_2$vireo_assignment == "doublet",]
  result_2_vireo <- result_2_vireo %>% dplyr::select(demuxalot_best_assignment, demuxlet_best_assignment, souporcell_best_assignment)
  result_2_vireo$one <- 0
  result_2_vireo$one[result_2_vireo$demuxalot_best_assignment == "doublet"] <- 1
  result_2_vireo$two <- 0
  result_2_vireo$two[result_2_vireo$demuxlet_best_assignment == "doublet"] <- 1
  result_2_vireo$three <- 0
  result_2_vireo$three[result_2_vireo$souporcell_best_assignment == "doublet"] <- 1
  result_2_vireo <- result_2_vireo %>% dplyr::select(one, two, three)
  result_2_vireo$sum <- rowSums(result_2_vireo)
  result_2_vireo$tool <- "Vireo"
  ## Souporcell
  result_2_souporcell <- result_2[result_2$souporcell_assignment == "doublet",]
  result_2_souporcell <- result_2_souporcell %>% dplyr::select(demuxalot_best_assignment, demuxlet_best_assignment, vireo_best_assignment)
  result_2_souporcell$one <- 0
  result_2_souporcell$one[result_2_souporcell$demuxalot_best_assignment == "doublet"] <- 1
  result_2_souporcell$two <- 0
  result_2_souporcell$two[result_2_souporcell$demuxlet_best_assignment == "doublet"] <- 1
  result_2_souporcell$three <- 0
  result_2_souporcell$three[result_2_souporcell$vireo_best_assignment == "doublet"] <- 1
  result_2_souporcell <- result_2_souporcell %>% dplyr::select(one, two, three)
  result_2_souporcell$sum <- rowSums(result_2_souporcell)
  result_2_souporcell$tool <- "Souporcell"
  ## Demuxlet
  result_2_freemuxlet <- result_2[result_2$demuxlet_assignment == "doublet",]
  result_2_freemuxlet <- result_2_freemuxlet %>% dplyr::select(demuxalot_best_assignment, souporcell_best_assignment, vireo_best_assignment)
  result_2_freemuxlet$one <- 0
  result_2_freemuxlet$one[result_2_freemuxlet$demuxalot_best_assignment == "doublet"] <- 1
  result_2_freemuxlet$two <- 0
  result_2_freemuxlet$two[result_2_freemuxlet$souporcell_best_assignment == "doublet"] <- 1
  result_2_freemuxlet$three <- 0
  result_2_freemuxlet$three[result_2_freemuxlet$vireo_best_assignment == "doublet"] <- 1
  result_2_freemuxlet <- result_2_freemuxlet %>% dplyr::select(one, two, three)
  result_2_freemuxlet$sum <- rowSums(result_2_freemuxlet)
  result_2_freemuxlet$tool <- "Demuxlet"

  ## merge and plot
  merge_df <- rbind(result_2_demuxalot, result_2_vireo, result_2_souporcell, result_2_freemuxlet)
  df2 <- merge_df %>% group_by(tool, sum) %>% count()
  df2 <- df2 %>%
  group_by(tool) %>%
  mutate(label_y = cumsum(n) - 0.5 * n)
    ggplot(df2, aes(x = tool, y = n, fill = as.character(sum), group = tool)) + 
    geom_bar(stat = "identity", position = "stack") +theme_classic() +
    geom_text(aes( y = label_y, label = n), vjust = 1.5, colour = "black") +
    geom_hline(yintercept = expected_doublets)+
    xlab("Demultiplexing tool") +
    ylab("Number of doublets") + 
    scale_fill_manual(values = c("lightgrey", "#ffb9b9", "#ee7272", "#a31818")) +
    theme(legend.position = "right")
  ggsave(paste(par_output_dir,"/step3","/Doublet_overlap_threshold.pdf", sep=""))

  ### Proportion agreement bar plot without threshold
   ## Demuxalot
  result_2_demuxalot <- result_2[result_2$demuxalot_best_assignment == "doublet",]
  result_2_demuxalot <- result_2_demuxalot %>% dplyr::select(vireo_best_assignment, demuxlet_best_assignment, souporcell_best_assignment)
  result_2_demuxalot$one <- 0
  result_2_demuxalot$one[result_2_demuxalot$vireo_best_assignment == "doublet"] <- 1
  result_2_demuxalot$two <- 0
  result_2_demuxalot$two[result_2_demuxalot$demuxlet_best_assignment == "doublet"] <- 1
  result_2_demuxalot$three <- 0
  result_2_demuxalot$three[result_2_demuxalot$souporcell_best_assignment == "doublet"] <- 1
  result_2_demuxalot <- result_2_demuxalot %>% dplyr::select(one, two, three)
  result_2_demuxalot$sum <- rowSums(result_2_demuxalot)
  result_2_demuxalot$tool <- "Demuxalot"
  ## Vireo
  result_2_vireo <- result_2[result_2$vireo_best_assignment == "doublet",]
  result_2_vireo <- result_2_vireo %>% dplyr::select(demuxalot_best_assignment, demuxlet_best_assignment, souporcell_best_assignment)
  result_2_vireo$one <- 0
  result_2_vireo$one[result_2_vireo$demuxalot_best_assignment == "doublet"] <- 1
  result_2_vireo$two <- 0
  result_2_vireo$two[result_2_vireo$demuxlet_best_assignment == "doublet"] <- 1
  result_2_vireo$three <- 0
  result_2_vireo$three[result_2_vireo$souporcell_best_assignment == "doublet"] <- 1
  result_2_vireo <- result_2_vireo %>% dplyr::select(one, two, three)
  result_2_vireo$sum <- rowSums(result_2_vireo)
  result_2_vireo$tool <- "Vireo"
  ## Souporcell
  result_2_souporcell <- result_2[result_2$souporcell_best_assignment == "doublet",]
  result_2_souporcell <- result_2_souporcell %>% dplyr::select(demuxalot_best_assignment, demuxlet_best_assignment, vireo_best_assignment)
  result_2_souporcell$one <- 0
  result_2_souporcell$one[result_2_souporcell$demuxalot_best_assignment == "doublet"] <- 1
  result_2_souporcell$two <- 0
  result_2_souporcell$two[result_2_souporcell$demuxlet_best_assignment == "doublet"] <- 1
  result_2_souporcell$three <- 0
  result_2_souporcell$three[result_2_souporcell$vireo_best_assignment == "doublet"] <- 1
  result_2_souporcell <- result_2_souporcell %>% dplyr::select(one, two, three)
  result_2_souporcell$sum <- rowSums(result_2_souporcell)
  result_2_souporcell$tool <- "Souporcell"
  ## Demuxlet
  result_2_freemuxlet <- result_2[result_2$demuxlet_best_assignment == "doublet",]
  result_2_freemuxlet <- result_2_freemuxlet %>% dplyr::select(demuxalot_best_assignment, souporcell_best_assignment, vireo_best_assignment)
  result_2_freemuxlet$one <- 0
  result_2_freemuxlet$one[result_2_freemuxlet$demuxalot_best_assignment == "doublet"] <- 1
  result_2_freemuxlet$two <- 0
  result_2_freemuxlet$two[result_2_freemuxlet$souporcell_best_assignment == "doublet"] <- 1
  result_2_freemuxlet$three <- 0
  result_2_freemuxlet$three[result_2_freemuxlet$vireo_best_assignment == "doublet"] <- 1
  result_2_freemuxlet <- result_2_freemuxlet %>% dplyr::select(one, two, three)
  result_2_freemuxlet$sum <- rowSums(result_2_freemuxlet)
  result_2_freemuxlet$tool <- "Demuxlet"

  ## merge and plot
  merge_df <- rbind(result_2_demuxalot, result_2_vireo, result_2_souporcell, result_2_freemuxlet)
  df2 <- merge_df %>% group_by(tool, sum) %>% count()
  df2 <- df2 %>%
  group_by(tool) %>%
  mutate(label_y = cumsum(n) - 0.5 * n)
    ggplot(df2, aes(x = tool, y = n, fill = as.character(sum), group = tool)) + 
    geom_bar(stat = "identity", position = "stack") +theme_classic() +
    geom_text(aes( y = label_y, label = n), vjust = 1.5, colour = "black") +
    geom_hline(yintercept = expected_doublets)+
    xlab("Demultiplexing tool") +
    ylab("Number of doublets (no threshold)") + 
    scale_fill_manual(values = c("lightgrey", "#ffb9b9", "#ee7272", "#a31818")) +
    theme(legend.position = "right")
  ggsave(paste(par_output_dir,"/step3","/Doublet_overlap_no_threshold.pdf", sep=""))

  ### Number of ensemblex droplets with EID of each tool with threshold
  ## Souporcell
  result_2_temp_souporcell <- result_2
  result_2_temp_souporcell$ensemblex_assignment[result_2_temp_souporcell$souporcell_assignment == "doublet"] <- "doublet"
  n_souporcell <- result_2_temp_souporcell[result_2_temp_souporcell$ensemblex_assignment == "doublet",] %>% nrow()

  ## Vireo
  result_2_temp_vireo <- result_2
  result_2_temp_vireo$ensemblex_assignment[result_2_temp_vireo$vireo_assignment == "doublet"] <- "doublet"
  n_vireo <- result_2_temp_vireo[result_2_temp_vireo$ensemblex_assignment == "doublet",] %>% nrow()
  
  ## Demuxlet
  result_2_temp_demuxlet <- result_2
  result_2_temp_demuxlet$ensemblex_assignment[result_2_temp_demuxlet$demuxlet_assignment == "doublet"] <- "doublet"
  n_demuxlet <- result_2_temp_demuxlet[result_2_temp_demuxlet$ensemblex_assignment == "doublet",] %>% nrow()
  
  ## Demuxalot
  result_2_temp_demuxalot <- result_2
  result_2_temp_demuxalot$ensemblex_assignment[result_2_temp_demuxalot$demuxalot_assignment == "doublet"] <- "doublet"
  n_demuxalot <- result_2_temp_demuxalot[result_2_temp_demuxalot$ensemblex_assignment == "doublet",] %>% nrow()

  ## Plot
  df <- data.frame(Tool = c("Demuxalot", "Vireo", "Souporcell", "Demuxlet"),
                    n_doublets = c(n_demuxalot, n_vireo, n_souporcell, n_demuxlet))

   ggplot(df, aes(x = Tool, y = n_doublets, label = n_doublets, fill = Tool)) + 
    geom_bar(stat = "identity") +theme_classic() +
    geom_hline(yintercept = expected_doublets)+
    geom_text() +
    xlab("Demultiplexing tool") +
    ylab("Number of doublets (threshold)") + 
    scale_fill_manual(values = c("#d95f02", "#e6ab02", "#7570b3", "#66a61e")) +
    theme(legend.position = "right")
  ggsave(paste(par_output_dir,"/step3","/Number_ensemblex_doublets_EID_threshold.pdf", sep=""))

  ### number of ensemblex droplets with EID of each tool with out threshold
  ## Souporcell
  result_2_temp_souporcell <- result_2
  result_2_temp_souporcell$ensemblex_assignment[result_2_temp_souporcell$souporcell_best_assignment == "doublet"] <- "doublet"
  n_souporcell <- result_2_temp_souporcell[result_2_temp_souporcell$ensemblex_assignment == "doublet",] %>% nrow()

  ## Vireo
  result_2_temp_vireo <- result_2
  result_2_temp_vireo$ensemblex_assignment[result_2_temp_vireo$vireo_best_assignment == "doublet"] <- "doublet"
  n_vireo <- result_2_temp_vireo[result_2_temp_vireo$ensemblex_assignment == "doublet",] %>% nrow()
  
  ## Demuxlet
  result_2_temp_demuxlet <- result_2
  result_2_temp_demuxlet$ensemblex_assignment[result_2_temp_demuxlet$demuxlet_best_assignment == "doublet"] <- "doublet"
  n_demuxlet <- result_2_temp_demuxlet[result_2_temp_demuxlet$ensemblex_assignment == "doublet",] %>% nrow()
  
  ## Demuxalot
  result_2_temp_demuxalot <- result_2
  result_2_temp_demuxalot$ensemblex_assignment[result_2_temp_demuxalot$demuxalot_best_assignment == "doublet"] <- "doublet"
  n_demuxalot <- result_2_temp_demuxalot[result_2_temp_demuxalot$ensemblex_assignment == "doublet",] %>% nrow()

  ## Plot
  df <- data.frame(Tool = c("Demuxalot", "Vireo", "Souporcell", "Demuxlet"),
                    n_doublets = c(n_demuxalot, n_vireo, n_souporcell, n_demuxlet))

   ggplot(df, aes(x = Tool, y = n_doublets, label = n_doublets, fill = Tool)) + 
    geom_bar(stat = "identity") +theme_classic() +
    geom_hline(yintercept = expected_doublets)+
    geom_text() +
    xlab("Demultiplexing tool") +
    ylab("Number of doublets (threshold)") + 
    scale_fill_manual(values = c("#d95f02", "#e6ab02", "#7570b3", "#66a61e")) +
    theme(legend.position = "right")
  ggsave(paste(par_output_dir,"/step3","/Number_ensemblex_doublets_EID_no_threshold.pdf", sep=""))
  
 ####################################################################################################################################################################################
 ## Ensemble inpendent doublet detection ## 
 ####################################################################################################################################################################################
  ## Threshold
  ## If Souporcell says doublet do doublet
  if ((tolower(par_doublet_Souporcell_threshold))=="yes"){
    message("Labelling all doublets identified by Souporcell (threshold) as doublets.")
  result_2$ensemblex_assignment[result_2$souporcell_assignment == "doublet"] <- "doublet"
  }
  ## If Demuxlet says doublet do doublet
  if ((tolower(par_doublet_Demuxlet_threshold))=="yes"){
    message("Labelling all doublets identified by Demuxlet (threshold) as doublets.")
  result_2$ensemblex_assignment[result_2$freemuxlet_assignment == "doublet"] <- "doublet"
  }
  ## If Vireo says doublet do doublet
  if ((tolower(par_doublet_Vireo_threshold))=="yes"){
    message("Labelling all doublets identified by Vireo (threshold) as doublets.")
  result_2$ensemblex_assignment[result_2$vireo_assignment  == "doublet"] <- "doublet"
  }
  ## If Demuxalot says doublet do doublet
  if ((tolower(par_doublet_Demuxalot_threshold))=="yes") {
    message("Labelling all doublets identified by Demuxalot (threshold) as doublets.")
  result_2$ensemblex_assignment[result_2$demuxalot_assignment  == "doublet"] <- "doublet"
  }

  ### No threshold
  ## If Souporcell says doublet do doublet
  if ((tolower(par_doublet_Souporcell_no_threshold))=="yes") {
     message("Labelling all doublets identified by Souporcell (no threshold) as doublets.")
  result_2$ensemblex_assignment[result_2$souporcell_best_assignment == "doublet"] <- "doublet"
  }
  ## If Demuxlet says doublet do doublet
  if ((tolower(par_doublet_Demuxlet_no_threshold))=="yes") {
     message("Labelling all doublets identified by Demuxlet (no threshold) as doublets.")
  result_2$ensemblex_assignment[result_2$demuxlet_best_assignment == "doublet"] <- "doublet"
  }
  ## If Vireo says doublet do doublet
  if ((tolower(par_doublet_Vireo_no_threshold))=="yes") {
     message("Labelling all doublets identified by Vireo (no threshold) as doublets.")
  result_2$ensemblex_assignment[result_2$vireo_best_assignment == "doublet"] <- "doublet"
  }
  ## If Demuxalot says doublet do doublet
  if ((tolower(par_doublet_Demuxalot_no_threshold))=="yes") {
     message("Labelling all doublets identified by Demuxalot (no threshold) as doublets.")
  result_2$ensemblex_assignment[result_2$demuxalot_best_assignment == "doublet"] <- "doublet"
  }

  ## Write csv file
  write.csv(result_2, paste(par_output_dir,"/step3",'/Step3_cell_assignment.csv', sep=""))

  result_2
}
ensemble_independent_doublet_detections_prelim <- function(result_2, par_output_dir){
  ## Set seed
  set.seed(1234)
  
  ## Create an output directory
  dir.create(paste(par_output_dir,"/step3",sep=''))

  expected_doublets <- nrow(result_2)*par_expected_doublet_rate

  ### Proportion agreement bar plot with threshold
  ## Demuxalot
  result_2_demuxalot <- result_2[result_2$demuxalot_assignment == "doublet",]
  result_2_demuxalot <- result_2_demuxalot %>% dplyr::select(vireo_best_assignment, demuxlet_best_assignment, souporcell_best_assignment)
  result_2_demuxalot$one <- 0
  result_2_demuxalot$one[result_2_demuxalot$vireo_best_assignment == "doublet"] <- 1
  result_2_demuxalot$two <- 0
  result_2_demuxalot$two[result_2_demuxalot$demuxlet_best_assignment == "doublet"] <- 1
  result_2_demuxalot$three <- 0
  result_2_demuxalot$three[result_2_demuxalot$souporcell_best_assignment == "doublet"] <- 1
  result_2_demuxalot <- result_2_demuxalot %>% dplyr::select(one, two, three)
  result_2_demuxalot$sum <- rowSums(result_2_demuxalot)
  result_2_demuxalot$tool <- "Demuxalot"
  ## Vireo
  result_2_vireo <- result_2[result_2$vireo_assignment == "doublet",]
  result_2_vireo <- result_2_vireo %>% dplyr::select(demuxalot_best_assignment, demuxlet_best_assignment, souporcell_best_assignment)
  result_2_vireo$one <- 0
  result_2_vireo$one[result_2_vireo$demuxalot_best_assignment == "doublet"] <- 1
  result_2_vireo$two <- 0
  result_2_vireo$two[result_2_vireo$demuxlet_best_assignment == "doublet"] <- 1
  result_2_vireo$three <- 0
  result_2_vireo$three[result_2_vireo$souporcell_best_assignment == "doublet"] <- 1
  result_2_vireo <- result_2_vireo %>% dplyr::select(one, two, three)
  result_2_vireo$sum <- rowSums(result_2_vireo)
  result_2_vireo$tool <- "Vireo"
  ## Souporcell
  result_2_souporcell <- result_2[result_2$souporcell_assignment == "doublet",]
  result_2_souporcell <- result_2_souporcell %>% dplyr::select(demuxalot_best_assignment, demuxlet_best_assignment, vireo_best_assignment)
  result_2_souporcell$one <- 0
  result_2_souporcell$one[result_2_souporcell$demuxalot_best_assignment == "doublet"] <- 1
  result_2_souporcell$two <- 0
  result_2_souporcell$two[result_2_souporcell$demuxlet_best_assignment == "doublet"] <- 1
  result_2_souporcell$three <- 0
  result_2_souporcell$three[result_2_souporcell$vireo_best_assignment == "doublet"] <- 1
  result_2_souporcell <- result_2_souporcell %>% dplyr::select(one, two, three)
  result_2_souporcell$sum <- rowSums(result_2_souporcell)
  result_2_souporcell$tool <- "Souporcell"
  ## Demuxlet
  result_2_freemuxlet <- result_2[result_2$demuxlet_assignment == "doublet",]
  result_2_freemuxlet <- result_2_freemuxlet %>% dplyr::select(demuxalot_best_assignment, souporcell_best_assignment, vireo_best_assignment)
  result_2_freemuxlet$one <- 0
  result_2_freemuxlet$one[result_2_freemuxlet$demuxalot_best_assignment == "doublet"] <- 1
  result_2_freemuxlet$two <- 0
  result_2_freemuxlet$two[result_2_freemuxlet$souporcell_best_assignment == "doublet"] <- 1
  result_2_freemuxlet$three <- 0
  result_2_freemuxlet$three[result_2_freemuxlet$vireo_best_assignment == "doublet"] <- 1
  result_2_freemuxlet <- result_2_freemuxlet %>% dplyr::select(one, two, three)
  result_2_freemuxlet$sum <- rowSums(result_2_freemuxlet)
  result_2_freemuxlet$tool <- "Demuxlet"

  ## merge and plot
  merge_df <- rbind(result_2_demuxalot, result_2_vireo, result_2_souporcell, result_2_freemuxlet)
  df2 <- merge_df %>% group_by(tool, sum) %>% count()
  df2 <- df2 %>%
  group_by(tool) %>%
  mutate(label_y = cumsum(n) - 0.5 * n)
    ggplot(df2, aes(x = tool, y = n, fill = as.character(sum), group = tool)) + 
    geom_bar(stat = "identity", position = "stack") +theme_classic() +
    geom_text(aes( y = label_y, label = n), vjust = 1.5, colour = "black") +
    geom_hline(yintercept = expected_doublets)+
    xlab("Demultiplexing tool") +
    ylab("Number of doublets") + 
    scale_fill_manual(values = c("lightgrey", "#ffb9b9", "#ee7272", "#a31818")) +
    theme(legend.position = "right")
  ggsave(paste(par_output_dir,"/step3","/Doublet_overlap_threshold.pdf", sep=""))

  ### Proportion agreement bar plot without threshold
   ## Demuxalot
  result_2_demuxalot <- result_2[result_2$demuxalot_best_assignment == "doublet",]
  result_2_demuxalot <- result_2_demuxalot %>% dplyr::select(vireo_best_assignment, demuxlet_best_assignment, souporcell_best_assignment)
  result_2_demuxalot$one <- 0
  result_2_demuxalot$one[result_2_demuxalot$vireo_best_assignment == "doublet"] <- 1
  result_2_demuxalot$two <- 0
  result_2_demuxalot$two[result_2_demuxalot$demuxlet_best_assignment == "doublet"] <- 1
  result_2_demuxalot$three <- 0
  result_2_demuxalot$three[result_2_demuxalot$souporcell_best_assignment == "doublet"] <- 1
  result_2_demuxalot <- result_2_demuxalot %>% dplyr::select(one, two, three)
  result_2_demuxalot$sum <- rowSums(result_2_demuxalot)
  result_2_demuxalot$tool <- "Demuxalot"
  ## Vireo
  result_2_vireo <- result_2[result_2$vireo_best_assignment == "doublet",]
  result_2_vireo <- result_2_vireo %>% dplyr::select(demuxalot_best_assignment, demuxlet_best_assignment, souporcell_best_assignment)
  result_2_vireo$one <- 0
  result_2_vireo$one[result_2_vireo$demuxalot_best_assignment == "doublet"] <- 1
  result_2_vireo$two <- 0
  result_2_vireo$two[result_2_vireo$demuxlet_best_assignment == "doublet"] <- 1
  result_2_vireo$three <- 0
  result_2_vireo$three[result_2_vireo$souporcell_best_assignment == "doublet"] <- 1
  result_2_vireo <- result_2_vireo %>% dplyr::select(one, two, three)
  result_2_vireo$sum <- rowSums(result_2_vireo)
  result_2_vireo$tool <- "Vireo"
  ## Souporcell
  result_2_souporcell <- result_2[result_2$souporcell_best_assignment == "doublet",]
  result_2_souporcell <- result_2_souporcell %>% dplyr::select(demuxalot_best_assignment, demuxlet_best_assignment, vireo_best_assignment)
  result_2_souporcell$one <- 0
  result_2_souporcell$one[result_2_souporcell$demuxalot_best_assignment == "doublet"] <- 1
  result_2_souporcell$two <- 0
  result_2_souporcell$two[result_2_souporcell$demuxlet_best_assignment == "doublet"] <- 1
  result_2_souporcell$three <- 0
  result_2_souporcell$three[result_2_souporcell$vireo_best_assignment == "doublet"] <- 1
  result_2_souporcell <- result_2_souporcell %>% dplyr::select(one, two, three)
  result_2_souporcell$sum <- rowSums(result_2_souporcell)
  result_2_souporcell$tool <- "Souporcell"
  ## Demuxlet
  result_2_freemuxlet <- result_2[result_2$demuxlet_best_assignment == "doublet",]
  result_2_freemuxlet <- result_2_freemuxlet %>% dplyr::select(demuxalot_best_assignment, souporcell_best_assignment, vireo_best_assignment)
  result_2_freemuxlet$one <- 0
  result_2_freemuxlet$one[result_2_freemuxlet$demuxalot_best_assignment == "doublet"] <- 1
  result_2_freemuxlet$two <- 0
  result_2_freemuxlet$two[result_2_freemuxlet$souporcell_best_assignment == "doublet"] <- 1
  result_2_freemuxlet$three <- 0
  result_2_freemuxlet$three[result_2_freemuxlet$vireo_best_assignment == "doublet"] <- 1
  result_2_freemuxlet <- result_2_freemuxlet %>% dplyr::select(one, two, three)
  result_2_freemuxlet$sum <- rowSums(result_2_freemuxlet)
  result_2_freemuxlet$tool <- "Demuxlet"

  ## merge and plot
  merge_df <- rbind(result_2_demuxalot, result_2_vireo, result_2_souporcell, result_2_freemuxlet)
  df2 <- merge_df %>% group_by(tool, sum) %>% count()
  df2 <- df2 %>%
  group_by(tool) %>%
  mutate(label_y = cumsum(n) - 0.5 * n)
    ggplot(df2, aes(x = tool, y = n, fill = as.character(sum), group = tool)) + 
    geom_bar(stat = "identity", position = "stack") +theme_classic() +
    geom_text(aes( y = label_y, label = n), vjust = 1.5, colour = "black") +
    geom_hline(yintercept = expected_doublets)+
    xlab("Demultiplexing tool") +
    ylab("Number of doublets (no threshold)") + 
    scale_fill_manual(values = c("lightgrey", "#ffb9b9", "#ee7272", "#a31818")) +
    theme(legend.position = "right")
  ggsave(paste(par_output_dir,"/step3","/Doublet_overlap_no_threshold.pdf", sep=""))

  ### Number of ensemblex droplets with EID of each tool with threshold
  ## Souporcell
  result_2_temp_souporcell <- result_2
  result_2_temp_souporcell$ensemblex_assignment[result_2_temp_souporcell$souporcell_assignment == "doublet"] <- "doublet"
  n_souporcell <- result_2_temp_souporcell[result_2_temp_souporcell$ensemblex_assignment == "doublet",] %>% nrow()

  ## Vireo
  result_2_temp_vireo <- result_2
  result_2_temp_vireo$ensemblex_assignment[result_2_temp_vireo$vireo_assignment == "doublet"] <- "doublet"
  n_vireo <- result_2_temp_vireo[result_2_temp_vireo$ensemblex_assignment == "doublet",] %>% nrow()
  
  ## Demuxlet
  result_2_temp_demuxlet <- result_2
  result_2_temp_demuxlet$ensemblex_assignment[result_2_temp_demuxlet$demuxlet_assignment == "doublet"] <- "doublet"
  n_demuxlet <- result_2_temp_demuxlet[result_2_temp_demuxlet$ensemblex_assignment == "doublet",] %>% nrow()
  
  ## Demuxalot
  result_2_temp_demuxalot <- result_2
  result_2_temp_demuxalot$ensemblex_assignment[result_2_temp_demuxalot$demuxalot_assignment == "doublet"] <- "doublet"
  n_demuxalot <- result_2_temp_demuxalot[result_2_temp_demuxalot$ensemblex_assignment == "doublet",] %>% nrow()

  ## Plot
  df <- data.frame(Tool = c("Demuxalot", "Vireo", "Souporcell", "Demuxlet"),
                    n_doublets = c(n_demuxalot, n_vireo, n_souporcell, n_demuxlet))

   ggplot(df, aes(x = Tool, y = n_doublets, label = n_doublets, fill = Tool)) + 
    geom_bar(stat = "identity") +theme_classic() +
    geom_hline(yintercept = expected_doublets)+
    geom_text() +
    xlab("Demultiplexing tool") +
    ylab("Number of doublets (threshold)") + 
    scale_fill_manual(values = c("#d95f02", "#e6ab02", "#7570b3", "#66a61e")) +
    theme(legend.position = "right")
  ggsave(paste(par_output_dir,"/step3","/Number_ensemblex_doublets_EID_threshold.pdf", sep=""))

  ### number of ensemblex droplets with EID of each tool with out threshold
  ## Souporcell
  result_2_temp_souporcell <- result_2
  result_2_temp_souporcell$ensemblex_assignment[result_2_temp_souporcell$souporcell_best_assignment == "doublet"] <- "doublet"
  n_souporcell <- result_2_temp_souporcell[result_2_temp_souporcell$ensemblex_assignment == "doublet",] %>% nrow()

  ## Vireo
  result_2_temp_vireo <- result_2
  result_2_temp_vireo$ensemblex_assignment[result_2_temp_vireo$vireo_best_assignment == "doublet"] <- "doublet"
  n_vireo <- result_2_temp_vireo[result_2_temp_vireo$ensemblex_assignment == "doublet",] %>% nrow()
  
  ## Demuxlet
  result_2_temp_demuxlet <- result_2
  result_2_temp_demuxlet$ensemblex_assignment[result_2_temp_demuxlet$demuxlet_best_assignment == "doublet"] <- "doublet"
  n_demuxlet <- result_2_temp_demuxlet[result_2_temp_demuxlet$ensemblex_assignment == "doublet",] %>% nrow()
  
  ## Demuxalot
  result_2_temp_demuxalot <- result_2
  result_2_temp_demuxalot$ensemblex_assignment[result_2_temp_demuxalot$demuxalot_best_assignment == "doublet"] <- "doublet"
  n_demuxalot <- result_2_temp_demuxalot[result_2_temp_demuxalot$ensemblex_assignment == "doublet",] %>% nrow()

  ## Plot
  df <- data.frame(Tool = c("Demuxalot", "Vireo", "Souporcell", "Demuxlet"),
                    n_doublets = c(n_demuxalot, n_vireo, n_souporcell, n_demuxlet))

   ggplot(df, aes(x = Tool, y = n_doublets, label = n_doublets, fill = Tool)) + 
    geom_bar(stat = "identity") +theme_classic() +
    geom_hline(yintercept = expected_doublets)+
    geom_text() +
    xlab("Demultiplexing tool") +
    ylab("Number of doublets (threshold)") + 
    scale_fill_manual(values = c("#d95f02", "#e6ab02", "#7570b3", "#66a61e")) +
    theme(legend.position = "right")
  ggsave(paste(par_output_dir,"/step3","/Number_ensemblex_doublets_EID_no_threshold.pdf", sep=""))
}

###########################################################################################################################
# CONFIDENCE-SCORE
###########################################################################################################################
## FUNCTION
confidence_score <- function(result_2, par_output_dir, par_sample_size){
   
  ## Set seed
  set.seed(1234)
  
  ## Create an output directory
  dir.create(paste(par_output_dir,"/confidence",sep=''))
  
  #### Calculate AUC singlet detection using consensus cells as proxy for ground truth
  ### Vireo 
  eval_df <- result_2
  eval_df_lim <- subset(eval_df, souporcell_best_assignment == demuxlet_best_assignment &
                          souporcell_best_assignment == demuxalot_best_assignment &
                          souporcell_best_assignment != "unassigned" &
                          souporcell_best_assignment != "doublet")
  
  eval_df_lim$consensus_eval_ROC <- "bad"
  eval_df_lim$consensus_eval_ROC[eval_df_lim$vireo_best_assignment == eval_df_lim$souporcell_best_assignment] <- "good"
  
  ## Check if we have sufficient values to compute AUC. 
  temp_neg = subset(eval_df_lim, consensus_eval_ROC == "bad")
  neg <- nrow(temp_neg)
  temp_pos = subset(eval_df_lim, consensus_eval_ROC == "good")
  pos <- nrow(temp_pos)
  
    if (pos != 0 & neg != 0){
    if (neg <=(0.01*pos) | pos <=(0.01*neg) ){
        print("Limited droplets obtained to compute Vireo AUC; results may not be reflective of true AUC.")
    } else {
        print("Sufficient droplets obtained to compute Vireo AUC.")
    }

  roc_empirical <- rocit(score = log(eval_df_lim$vireo_max_probability), class = eval_df_lim$consensus_eval_ROC, 
                         negref = "bad") 
  print(summary(roc_empirical))

  vireo_AUC <- roc_empirical$AUC  
  vireo_AUC_singlet <- vireo_AUC
  
  } else {
    print(paste0("Insufficient droplets to compute Vireo AUC. Observed ", neg, " incorrectly classified droplets and ", pos, " correctly classified droplets. Setting Vireo AUC as 0.5 for confidence score computation."))
    vireo_AUC_singlet <- 0.5
  }

  ### Demuxlet 
  eval_df <- result_2
  eval_df_lim <- subset(eval_df, souporcell_best_assignment == vireo_best_assignment &
                          souporcell_best_assignment == demuxalot_best_assignment &
                          souporcell_best_assignment != "unassigned" &
                          souporcell_best_assignment != "doublet")
  
  eval_df_lim$consensus_eval_ROC <- "bad"
  eval_df_lim$consensus_eval_ROC[eval_df_lim$demuxlet_best_assignment == eval_df_lim$souporcell_best_assignment] <- "good"
  
  ## Check if we have sufficient values to compute AUC. 
  temp_neg = subset(eval_df_lim, consensus_eval_ROC == "bad")
  neg <- nrow(temp_neg)
  temp_pos = subset(eval_df_lim, consensus_eval_ROC == "good")
  pos <- nrow(temp_pos)
  
    if (pos != 0 & neg != 0){
    if (neg <=(0.01*pos) | pos <=(0.01*neg) ){
        print("Limited droplets obtained to compute Demuxlet AUC; results may not be reflective of true AUC.")
    } else {
        print("Sufficient droplets obtained to compute Demuxlet AUC.")
    }
  
  roc_empirical <- rocit(score = log(eval_df_lim$demuxlet_max_probability), class = eval_df_lim$consensus_eval_ROC, 
                         negref = "bad") 
  print(summary(roc_empirical))

  freemuxlet_AUC <- roc_empirical$AUC  
  freemuxlet_AUC_singlet <- freemuxlet_AUC
  
    } else {
    print(paste0("Insufficient droplets to compute Demuxlet AUC. Observed ", neg, " incorrectly classified droplets and ", pos, " correctly classified droplets. Setting Demuxlet AUC as 0.5 for confidence score computation."))
    freemuxlet_AUC_singlet <- 0.5
  }
  
  ### Demuxalot 
  eval_df <- result_2
  eval_df_lim <- subset(eval_df, souporcell_best_assignment == vireo_best_assignment &
                          souporcell_best_assignment == demuxlet_best_assignment &
                          souporcell_best_assignment != "unassigned" &
                          souporcell_best_assignment != "doublet")
  
  eval_df_lim$consensus_eval_ROC <- "bad"
  eval_df_lim$consensus_eval_ROC[eval_df_lim$demuxalot_best_assignment == eval_df_lim$souporcell_best_assignment] <- "good"

  ## Check if we have sufficient values to compute AUC. 
  temp_neg = subset(eval_df_lim, consensus_eval_ROC == "bad")
  neg <- nrow(temp_neg)
  temp_pos = subset(eval_df_lim, consensus_eval_ROC == "good")
  pos <- nrow(temp_pos)
  
    if (pos != 0 & neg != 0){
    if (neg <=(0.01*pos) | pos <=(0.01*neg) ){
        print("Limited droplets obtained to compute Demuxalot AUC; results may not be reflective of true AUC.")
    } else {
        print("Sufficient droplets obtained to compute Demuxalot AUC.")
    }
  
  roc_empirical <- rocit(score = log(eval_df_lim$demuxalot_max_probability), class = eval_df_lim$consensus_eval_ROC,
                         negref = "bad") 
  print(summary(roc_empirical))

  demuxalot_AUC <- roc_empirical$AUC
  demuxalot_AUC_singlet <- demuxalot_AUC
  
      } else {
    print(paste0("Insufficient droplets to compute Demuxalot AUC. Observed ", neg, " incorrectly classified droplets and ", pos, " correctly classified droplets. Setting Demuxalot AUC as 0.5 for confidence score computation."))
    demuxalot_AUC_singlet <- 0.5
  }
  
  ### Souporcell 
  eval_df <- result_2
  eval_df_lim <- subset(eval_df, demuxalot_best_assignment == vireo_best_assignment &
                          demuxalot_best_assignment == demuxlet_best_assignment &
                          demuxalot_best_assignment != "unassigned" &
                          demuxalot_best_assignment != "doublet")
  
  eval_df_lim$consensus_eval_ROC <- "bad"
  eval_df_lim$consensus_eval_ROC[eval_df_lim$demuxalot_best_assignment == eval_df_lim$souporcell_best_assignment] <- "good"
  
  ## Check if we have sufficient values to compute AUC. 
  temp_neg = subset(eval_df_lim, consensus_eval_ROC == "bad")
  neg <- nrow(temp_neg)
  temp_pos = subset(eval_df_lim, consensus_eval_ROC == "good")
  pos <- nrow(temp_pos)
  
    if (pos != 0 & neg != 0){
    if (neg <=(0.01*pos) | pos <=(0.01*neg) ){
        print("Limited droplets obtained to compute Souporcell AUC; results may not be reflective of true AUC.")
    } else {
        print("Sufficient droplets obtained to compute Souporcell AUC.")
    }

  roc_empirical <- rocit(score = log(1-(10^(eval_df_lim$souporcell_log_probability_singlet))), class = eval_df_lim$consensus_eval_ROC, ##change_here
                         negref = "bad") 
  print(summary(roc_empirical))
  souporcell_AUC <- roc_empirical$AUC  
  souporcell_AUC_singlet <- souporcell_AUC
  
  } else {
    print(paste0("Insufficient droplets to compute Demuxalot AUC. Observed ", neg, " incorrectly classified droplets and ", pos, " correctly classified droplets. Setting Demuxalot AUC as 0.5 for confidence score computation."))
    souporcell_AUC_singlet <- 0.5
  }
  
  ### Compute ensemblex singlet confidence 
  eval_df$ensemblex_singlet_confidence <- eval_df$ensemblex_probability
  
  ## Vireo
  eval_df$ensemblex_singlet_confidence[eval_df$vireo_singlet_probability >= 0.9 & eval_df$vireo_best_assignment != "doublet" &
                               eval_df$vireo_best_assignment == eval_df$ensemblex_assignment] <- eval_df$ensemblex_singlet_confidence[eval_df$vireo_singlet_probability >= 0.9 & eval_df$vireo_best_assignment != "doublet" &
                                                                                                                eval_df$vireo_best_assignment == eval_df$ensemblex_assignment] + vireo_AUC_singlet
  
  ## Demuxalot
  eval_df$ensemblex_singlet_confidence[eval_df$demuxalot_max_probability >= 0.9 & eval_df$demuxalot_best_assignment != "doublet" &
                               eval_df$demuxalot_best_assignment == eval_df$ensemblex_assignment] <- eval_df$ensemblex_singlet_confidence[eval_df$demuxalot_max_probability >= 0.9 & eval_df$demuxalot_best_assignment != "doublet" &
                                                                                                                    eval_df$demuxalot_best_assignment == eval_df$ensemblex_assignment] + demuxalot_AUC_singlet
  
  ## Freemuxlet
  eval_df$ensemblex_singlet_confidence[eval_df$demuxlet_assignment != "unassigned" & eval_df$demuxlet_best_assignment != "doublet" &
                               eval_df$demuxlet_best_assignment == eval_df$ensemblex_assignment] <- eval_df$ensemblex_singlet_confidence[eval_df$demuxlet_assignment != "unassigned" & eval_df$demuxlet_best_assignment != "doublet" &
                                                                                                                     eval_df$demuxlet_best_assignment == eval_df$ensemblex_assignment] + freemuxlet_AUC_singlet
  
  ## Souporcell
  eval_df$ensemblex_singlet_confidence[eval_df$souporcell_assignment != "unassigned" & eval_df$souporcell_best_assignment != "doublet" &
                               eval_df$souporcell_best_assignment == eval_df$ensemblex_assignment] <- eval_df$ensemblex_singlet_confidence[eval_df$souporcell_assignment != "unassigned" & eval_df$souporcell_best_assignment != "doublet" &
                                                                                                                     eval_df$souporcell_best_assignment == eval_df$ensemblex_assignment] + souporcell_AUC_singlet
  
  ## Unassignable cells
  eval_df$ensemblex_singlet_confidence[eval_df$vireo_n_vars == 0] <- eval_df$ensemblex_singlet_confidence[eval_df$vireo_n_vars == 0]/par_sample_size
  
  ## Consensus
  eval_df$ensemblex_singlet_confidence[eval_df$vireo_best_assignment == eval_df$demuxalot_best_assignment &
                               eval_df$vireo_best_assignment == eval_df$demuxlet_best_assignment &
                               eval_df$vireo_best_assignment == eval_df$souporcell_best_assignment &
                               eval_df$vireo_best_assignment != "doublet" &
                               eval_df$vireo_best_assignment != "unassigned"] <- eval_df$ensemblex_singlet_confidence[eval_df$vireo_best_assignment == eval_df$demuxalot_best_assignment &
                                                                                           eval_df$vireo_best_assignment == eval_df$demuxlet_best_assignment &
                                                                                           eval_df$vireo_best_assignment == eval_df$souporcell_best_assignment &
                                                                                           eval_df$vireo_best_assignment != "doublet" &
                                                                                           eval_df$vireo_best_assignment != "unassigned"] +1
  
  ## Set ensemblex best assignment
  eval_df$ensemblex_best_assignment <- eval_df$ensemblex_assignment 
  ## Set ensemblex assignment
  eval_df$ensemblex_assignment[eval_df$ensemblex_singlet_confidence < 1 & eval_df$ensemblex_assignment != "doublet" ] <- "unassigned" 


  eval_df <- dplyr::select(eval_df, c("barcode", "ensemblex_assignment","ensemblex_best_assignment", "ensemblex_probability", "ensemblex_singlet_confidence" , "vireo_assignment","souporcell_assignment","demuxlet_assignment","demuxalot_assignment","general_consensus" ,                
                                        "vireo_best_assignment","souporcell_best_assignment","demuxlet_best_assignment","demuxalot_best_assignment","vireo_singlet_probability","vireo_doublet_probability","vireo_n_vars","vireo_best_doublet",                
                                        "vireo_doublet_logLikRatio","souporcell_log_probability_singlet", "souporcell_log_probability_doublet", "demuxlet_n_snps", "demuxlet_n_reads","demuxlet_max_probability", "demuxlet_DIFF_LLK_SNG_DBL","demuxalot_max_probability"  ,       
                                        "demuxalot_doublet_probability", "vireo_max_probability","vireo_weighted_probability","demuxlet_weighted_probability","demuxalot_weighted_probability" ,  
                                        "souporcell_weighted_probability" ))

  write.csv(eval_df, paste(par_output_dir,"/confidence",'/ensemblex_final_cell_assignment.csv', sep=""))
  eval_df
}
