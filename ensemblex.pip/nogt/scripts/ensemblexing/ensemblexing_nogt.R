## load parameters
args = commandArgs(trailingOnly=TRUE)
home_dir=args[1]
output_dir=args[2]

source(paste0(home_dir,'/nogt/scripts/ensemblexing/functions_nogt.R'))

packages<-c('dplyr','tidyr','pdfCluster','data.table','readr','lubridate','tidyverse','moments','mousetrap','usethis','devtools','desc','kneedle','ROCit','ggplot2','factoextra','ggpubr','splines','stats','pathviewr')
invisible(lapply(packages, library, character.only = TRUE))

par_sample_size= as.numeric(args[3])
par_expected_doublet_rate= as.numeric(args[4])
par_merge_constituents= args[5]
par_probabilistic_weighted_ensemble= args[6]
par_preliminary_parameter_sweep= args[7]
par_graph_based_doublet_detection= args[8]
par_preliminary_ensemble_independent_doublet= args[9]
par_ensemble_independent_doublet= args[10]
par_doublet_Demuxalot_threshold= args[11]
par_doublet_Demuxalot_no_threshold= args[12]
par_doublet_Freemuxlet_threshold= args[13]
par_doublet_Freemuxlet_no_threshold= args[14]
par_doublet_Souporcell_threshold= args[15]
par_doublet_Souporcell_no_threshold= args[16]
par_doublet_Vireo_threshold= args[17]
par_doublet_Vireo_no_threshold= args[18]
par_compute_singlet_confidence= args[19]
par_ensemblex_nCD= as.numeric(args[20])
par_ensemblex_pT= as.numeric(args[21])

## Print the defined parameters by the user
print("Loaded parameter from config.ini")
print(paste0("par_sample_size=",par_sample_size))
print(paste0("par_expected_doublet_rate=",par_expected_doublet_rate))
print(paste0("par_merge_constituents=",par_merge_constituents))
print(paste0("par_probabilistic_weighted_ensemble=",par_probabilistic_weighted_ensemble))
print(paste0("par_preliminary_parameter_sweep=",par_preliminary_parameter_sweep))
print(paste0("par_graph_based_doublet_detection=",par_graph_based_doublet_detection))
print(paste0("par_ensemblex_nCD=",par_ensemblex_nCD))
print(paste0("par_ensemblex_pT=",par_ensemblex_pT))
print(paste0("par_preliminary_ensemble_independent_doublet=",par_preliminary_ensemble_independent_doublet))
print(paste0("par_ensemble_independent_doublet=",par_ensemble_independent_doublet))
print(paste0("par_doublet_Demuxalot_threshold=",par_doublet_Demuxalot_threshold))
print(paste0("par_doublet_Demuxalot_no_threshold=",par_doublet_Demuxalot_no_threshold))
print(paste0("par_doublet_Freemuxlet_threshold=",par_doublet_Freemuxlet_threshold))
print(paste0("par_doublet_Freemuxlet_no_threshold=",par_doublet_Freemuxlet_no_threshold))
print(paste0("par_doublet_Souporcell_threshold=",par_doublet_Souporcell_threshold))
print(paste0("par_doublet_Souporcell_no_threshold=",par_doublet_Souporcell_no_threshold))
print(paste0("par_doublet_Vireo_threshold=",par_doublet_Vireo_threshold))
print(paste0("par_doublet_Vireo_no_threshold=",par_doublet_Vireo_no_threshold))
print(paste0("par_compute_singlet_confidence=",par_compute_singlet_confidence))
print(paste0("par_doublet_Vireo_no_threshold=",par_doublet_Vireo_no_threshold))
print(paste0("par_compute_singlet_confidence=",par_compute_singlet_confidence))


par_output_dir<- paste0(output_dir,"/ensemblex")
## Demuxalot
par_demuxalot<- paste0(output_dir,"/demuxalot/Demuxalot_result.csv")
## Freemuxlet
par_freemuxlet<- paste0(output_dir,"/freemuxlet/outs.clust1.samples")
## Souporcell
par_souporcell<- paste0(output_dir,"/souporcell/clusters.tsv")
## Vireo
par_vireo<- paste0(output_dir,"/vireo/donor_ids.tsv")


###########################################################################################################################
# Print Warning message
###########################################################################################################################
message("You are running Ensemblex without prior genotype information. Please note that it is not recommended to run Ensemblex without prior genotype information for pools exceeding 16 multiplexed samples. Furthermore, we discourage using Ensemblex without prior genotype information if there is a strong imbalance of sample representation in the pool.")

###########################################################################################################################
# CONSTITUENT DATA PREPARATION AND MERGING
###########################################################################################################################
if (tolower(par_merge_constituents)=="yes"){
    message("Performing constituent data preparation and merge.")
    ## Import constituent tool output files
    vireo<- read.delim(par_vireo, header = T, sep = "\t")
    dim(vireo) 
    souporcell <- read.delim(par_souporcell, header = T, sep = "\t")
    dim(souporcell) 
    freemuxlet <-read.delim(par_freemuxlet, header = T, sep = "\t")
    dim(freemuxlet) 
    demuxalot <- read.delim(par_demuxalot, header = T, sep = ",", check.names=FALSE)
    colnames(demuxalot)[1] <- "X"
    dim(demuxalot) 
    ## check output files for duplicated barcodes
    message("Checking output files for duplicated barcodes.")
    ## Vireo
    if (length(unique(vireo$cell)) != nrow(vireo)){
    vireo <- vireo[!duplicated(vireo$cell),]
    temp <- vireo[duplicated(vireo$cell),]
        if (nrow(temp) == 0){
            print("Duplicated barcodes in Vireo resolved; recommended to investigate initial reason for barcode duplication.")
            print(paste0("Observed ", nrow(vireo), " unique barcodes in Vireo output."))
        }
    } else {
        print("No duplicated barcodes observed in Vireo output.")
        print(paste0("Observed ", nrow(vireo), " unique barcodes in Vireo output."))
    }
    ## Souporcell
    if (length(unique(souporcell$barcode)) != nrow(souporcell)){
    souporcell <- souporcell[!duplicated(souporcell$barcode),]
    temp <- souporcell[duplicated(souporcell$barcode),]
        if (nrow(temp) == 0){
            print("Duplicated barcodes in Souporcell resolved; recommended to investigate initial reason for barcode duplication.")
            print(paste0("Observed ", nrow(souporcell), " unique barcodes in Souporcell output."))
        }
    } else {
        print("No duplicated barcodes observed in Souporcell output.")
        print(paste0("Observed ", nrow(souporcell), " unique barcodes in Souporcell output."))
    }
    ## Freemuxlet
    if (length(unique(freemuxlet$BARCODE)) != nrow(freemuxlet)){
    freemuxlet <- freemuxlet[!duplicated(freemuxlet$BARCODE),]
    temp <- freemuxlet[duplicated(freemuxlet$BARCODE),]
        if (nrow(temp) == 0){
            print("Duplicated barcodes in Freemuxlet resolved; recommended to investigate initial reason for barcode duplication.")
            print(paste0("Observed ", nrow(freemuxlet), " unique barcodes in Freemuxlet output."))
        }
    } else {
        print("No duplicated barcodes observed in Freemuxlet output.")
        print(paste0("Observed ", nrow(freemuxlet), " unique barcodes in Freemuxlet output."))
    }
    ## Demuxalot
    if (length(unique(demuxalot$X)) != nrow(demuxalot)){
    demuxalot <- demuxalot[!duplicated(demuxalot$X),]
    temp <- demuxalot[duplicated(demuxalot$X),]
        if (nrow(temp) == 0){
            print("Duplicated barcodes in Demuxalot resolved; recommended to investigate initial reason for barcode duplication.")
            print(paste0("Observed ", nrow(demuxalot), " unique barcodes in Demuxalot output."))
        }
    } else {
        print("No duplicated barcodes observed in Demuxalot output.")
        print(paste0("Observed ", nrow(demuxalot), " unique barcodes in Demuxalot output."))
    }
    ## Make sure that all of the output files have the same number of droplets
    if (nrow(vireo) == nrow(freemuxlet) &
        nrow(vireo) == nrow(demuxalot) &
        nrow(vireo) == nrow(souporcell)) {
            message(paste0("All demultiplexing tools have reported the same number of cells: ", nrow(souporcell), "." ))
        } else {
            message(paste0("WARNING: Demultiplexing tools have reported different number of cells; Vireo: ", nrow(vireo), "; Freemuxlet: ", nrow(freemuxlet), "; Demuxalot: ", nrow(demuxalot), "; Souporcell: ", nrow(souporcell), "."))
        }
    #### Make sure that the Sample IDs from each output file have the same structure. We will replace "-" or "." with "_"
    ### Vireo
    ## fix potential "-"
    vireo$donor_id <- gsub(pattern = "-", replacement = "_",x = vireo$donor_id, fixed = TRUE)
    vireo$best_singlet <- gsub(pattern = "-", replacement = "_",x = vireo$best_singlet, fixed = TRUE)
    vireo$best_doublet <- gsub(pattern = "-", replacement = "_",x = vireo$best_doublet, fixed = TRUE)
    ## fix potential "."
    vireo$donor_id <- gsub(pattern = ".", replacement = "_",x = vireo$donor_id, fixed = TRUE)
    vireo$best_singlet <- gsub(pattern = ".", replacement = "_",x = vireo$best_singlet, fixed = TRUE)
    vireo$best_doublet <- gsub(pattern = ".", replacement = "_",x = vireo$best_doublet, fixed = TRUE)
    ### Demuxalot
    ## fix potential "-"
    colnames(demuxalot) <- gsub(pattern = "-", replacement = "_",x = colnames(demuxalot), fixed = TRUE)
    ## fix potential "."
    colnames(demuxalot) <- gsub(pattern = ".", replacement = "_",x = colnames(demuxalot), fixed = TRUE)
    ## Freemuxlet
    ## fix potential "-"
    freemuxlet$BEST.GUESS <- gsub(pattern = "-", replacement = "_",x = freemuxlet$BEST.GUESS, fixed = TRUE)
    freemuxlet$NEXT.GUESS <- gsub(pattern = "-", replacement = "_",x = freemuxlet$NEXT.GUESS, fixed = TRUE)
    freemuxlet$SNG.BEST.GUESS <- gsub(pattern = "-", replacement = "_",x = freemuxlet$SNG.BEST.GUESS, fixed = TRUE)
    freemuxlet$SNG.NEXT.GUESS <- gsub(pattern = "-", replacement = "_",x = freemuxlet$SNG.NEXT.GUESS, fixed = TRUE)
    freemuxlet$DBL.BEST.GUESS <- gsub(pattern = "-", replacement = "_",x = freemuxlet$DBL.BEST.GUESS, fixed = TRUE)
    ## fix potential "."
    freemuxlet$BEST.GUESS <- gsub(pattern = ".", replacement = "_",x = freemuxlet$BEST.GUESS, fixed = TRUE)
    freemuxlet$NEXT.GUESS <- gsub(pattern = ".", replacement = "_",x = freemuxlet$NEXT.GUESS, fixed = TRUE)
    freemuxlet$SNG.BEST.GUESS <- gsub(pattern = ".", replacement = "_",x = freemuxlet$SNG.BEST.GUESS, fixed = TRUE)
    freemuxlet$SNG.NEXT.GUESS <- gsub(pattern = ".", replacement = "_",x = freemuxlet$SNG.NEXT.GUESS, fixed = TRUE)
    freemuxlet$DBL.BEST.GUESS <- gsub(pattern = ".", replacement = "_",x = freemuxlet$DBL.BEST.GUESS, fixed = TRUE)
    ### Prepare outputs from each of the constituent demultiplexing tools for downstream analyses
    ## Demuxalot prep
    demuxalot <-  demuxalot_prep(demuxalot)
    ## Freemuxlet prep
    freemuxlet <- freemuxlet_prep(freemuxlet)
    ## Souporcell prep
    souporcell <- souporcell_prep(souporcell, demuxalot, freemuxlet)
    ## Vireo prep
    vireo <- vireo_prep(vireo,souporcell, demuxalot, freemuxlet)  
    ## Merge results and calculate consensus 
    merge_df <-  merge_concensus(vireo,souporcell, freemuxlet, demuxalot)
    ## Ensure that if outputs do not have matching droplet counts, the missing droplets are unassigned.
    merge_df$vireo_best_assignment[is.na(merge_df$vireo_best_assignment)] <- "unassigned"
    merge_df$souporcell_best_assignment[is.na(merge_df$souporcell_best_assignment)] <- "unassigned"
    merge_df$freemuxlet_best_assignment[is.na(merge_df$freemuxlet_best_assignment)] <- "unassigned"
    merge_df$demuxalot_best_assignment[is.na(merge_df$demuxalot_best_assignment)] <- "unassigned"
    ## Check if the merged csv file has duplicated barcodes
    if (length(unique(merge_df$barcode)) != nrow(merge_df)){
    merge_df <- merge_df[!duplicated(merge_df$barcode),]
    temp <- merge_df[duplicated(merge_df$barcode),]
        if (nrow(temp) == 0){
            print("Duplicated barcodes in Merged file resolved; recommended to investigate initial reason for barcode duplication.")
            print(paste0("Observed ", nrow(merge_df), " unique barcodes in Merged file."))
        }
    } else {
        print("No duplicated barcodes observed in Merged file.")
        print(paste0("Observed ", nrow(merge_df), " unique barcodes in Merged file."))
    }
    ## Write merge CSV file
    write.csv(merge_df, paste(par_output_dir,'/constituent_tool_merge.csv',sep =""), row.names=T)
} else {
    message("Skipping constituent data preparation and merge, and loading exisiting csv file.")
}


###########################################################################################################################
# PROBABILISTIC-WEIGHTED ENSEMBLE 
###########################################################################################################################
## Read in merge CSV file
if (tolower(par_probabilistic_weighted_ensemble)=="yes" & file.exists(paste0(par_output_dir,'/constituent_tool_merge.csv'))){
    message("Loading merged output file.")
    merge_df <- read.delim(paste(par_output_dir,'/constituent_tool_merge.csv',sep =""), header = T, sep = ",") 
} else if (tolower(par_probabilistic_weighted_ensemble)=="yes" & !file.exists(paste0(par_output_dir,'/constituent_tool_merge.csv'))) {
    stop('Exiting ensemblex; constituent_tool_merge.csv file cannot be found.')
} else {
    message("Skipping the probabilistic-weighted ensemble component of the ensemblex framework and loading existing csv file.")
}

## Probabilistic-weighted Ensemble (PWE) function
if (tolower(par_probabilistic_weighted_ensemble)=="yes"){
    message("Performing the probabilistic-weighted ensemble component of the ensemblex framework.")
## Run PWE
result_test <- BA_weight_consensus(merge_df, par_sample_size, par_output_dir)
} else {
    message("Skipping the probabilistic-weighted ensemble component of the ensemblex framework and loading existing csv file.")
}


###########################################################################################################################
# GRAPH-BASED DOUBLET DETECTION
###########################################################################################################################
## Read in inputs
if (tolower(par_preliminary_parameter_sweep)=="yes" | tolower(par_graph_based_doublet_detection)=="yes" ){
    if (file.exists(paste0(par_output_dir,"/step1",'/step1_cell_assignment.csv'))){
        message("Loading probabilistic-weighted ensemble output.")
        result_test <- read.delim(paste(par_output_dir,"/step1",'/step1_cell_assignment.csv', sep=""), header = T, sep = ",") 
        result_test <- result_test[,-1]
    } else {
        stop('Exiting ensemblex; step1_cell_assignment.csv cannot be found.')
    }
} else {
    message("Skipping the graph-based doublet detection component of the ensemblex framework and loading existing csv file.")
}

## perform graph-based doublet detection
if (tolower(par_preliminary_parameter_sweep)=="yes" & tolower(par_graph_based_doublet_detection)=="no"){
    message("Performing a preliminary parameter sweep for graph-based doublet detection.")
    ## apply graph-based doublet detection parameter sweep
    graph_based_doublet_detection_par_sweep(result_test, par_expected_doublet_rate, par_output_dir)

} else if (tolower(par_preliminary_parameter_sweep)=="yes" & tolower(par_graph_based_doublet_detection)=="yes" & exists("par_ensemblex_nCD") & exists("par_ensemblex_pT") ){
    message("Performing graph-based doublet detection with manually defined parameters.") 
    ## apply graph-based doublet detection with manually defined parameters
    graph_based_doublet_detection_manual_par(result_test, par_expected_doublet_rate, par_output_dir, par_ensemblex_pT, par_ensemblex_nCD)

} else if (tolower(par_preliminary_parameter_sweep)=="no" & tolower(par_graph_based_doublet_detection)=="yes"){
    message("Performing graph-based doublet detection with ensemblex-estimated parameters.") 
    ## apply graph-based doublet detection with estimated optimal parameters
    graph_based_doublet_detection_estimated_par(result_test, par_expected_doublet_rate, par_output_dir)
} else {
    message("Skipping graph-based doublet detection.")   
}

###########################################################################################################################
# ENSEMBLE-INDEPENDENT DOUBLET DETECTION
###########################################################################################################################
## Read in inputs
if (tolower(par_preliminary_ensemble_independent_doublet)=="yes" | tolower(par_ensemble_independent_doublet)=="yes"){
    if (file.exists(paste0(par_output_dir,"/step2",'/Step2_cell_assignment.csv'))){
        message("Loading graph-based doublet detection output.")
        result_2 <- read.delim(paste(par_output_dir,"/step2",'/Step2_cell_assignment.csv', sep=""), header = T, sep = ",") 
        result_2 <- result_2[,-1]
    } else if (file.exists(paste0(par_output_dir,"/step1",'/step1_cell_assignment.csv'))){
        message("Loading probabilistic-weighted ensemble output.")
        result_2 <- read.delim(paste(par_output_dir,"/step1",'/step1_cell_assignment.csv', sep=""), header = T, sep = ",")  
        result_2 <- result_2[,-1]
    } else {
        stop('Exiting ensemblex; Step2_cell_assignment.csv and step1_cell_assignment.csv files cannot be found.')
    }
} else {
    message("Skipping the ensemble-independent doublet detection component of the ensemblex framework and loading existing csv file.")
}

## Run ensemble-independet doublet detection
if (tolower(par_preliminary_ensemble_independent_doublet)=="yes" & tolower(par_ensemble_independent_doublet)=="yes"){
    ## Run EID
    result_2 <- ensemble_independent_doublet_detections(result_2, par_output_dir)
} else if (tolower(par_preliminary_ensemble_independent_doublet)=="yes" & tolower(par_ensemble_independent_doublet)=="no"){
    ## Run EID
    ensemble_independent_doublet_detections_prelim(result_2, par_output_dir)
} else if (tolower(par_preliminary_ensemble_independent_doublet)=="no" & tolower(par_ensemble_independent_doublet)=="yes") {
    ## Run EID
    result_2 <- ensemble_independent_doublet_detections(result_2, par_output_dir)
} else {
    message("Skipping ensemble independent doublet detection.")
}


###########################################################################################################################
# CONFIDENCE SCORE
###########################################################################################################################
if (tolower(par_compute_singlet_confidence)=="yes"){
    if (file.exists(paste0(par_output_dir,"/step3",'/Step3_cell_assignment.csv'))){
        message("Loading ensemble-independent doublet detection output.")
        result_2 <- read.delim(paste(par_output_dir,"/step3",'/Step3_cell_assignment.csv', sep=""), header = T, sep = ",") 
        result_2 <- result_2[,-1]
    } else if (file.exists(paste0(par_output_dir,"/step2",'/Step2_cell_assignment.csv'))){
        message("Loading graph-based doublet detection output.")
        result_2 <- read.delim(paste(par_output_dir,"/step2",'/Step2_cell_assignment.csv', sep=""), header = T, sep = ",") 
        result_2 <- result_2[,-1]
    } else if (file.exists(paste0(par_output_dir,"/step1",'/step1_cell_assignment.csv'))){
        message("Loading probabilistic-weighted ensemble output.")
        result_2 <- read.delim(paste(par_output_dir,"/step1",'/step1_cell_assignment.csv', sep=""), header = T, sep = ",")  
        result_2 <- result_2[,-1]
    } else {
        stop('Exiting ensemblex; Step3_cell_assignment.csv, Step2_cell_assignment.csv and step1_cell_assignment.csv files cannot be found.')
    }
} else {
    message("Skipping the ensemble-independent doublet detection component of the ensemblex framework and loading existing csv file.")
}

## Compute confidence score
if (tolower(par_compute_singlet_confidence)=="yes"){
## Compute ensemblex singlet confidence
eval_df <- confidence_score(result_2, par_output_dir,par_sample_size)
}
