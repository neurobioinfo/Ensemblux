#########
# NOTES #
#########
# This code was used to produce simulated pools with known ground truth sample labels using 80 independently-sequenced iPSC lines differentiated towards a dopaminergic neuronal state as part of the Foundational Data Initiative for Parkinson's Disease
# Simulated pools were produced using the synth_pool.py script produced by the developers of Vireo. 24 samples were multiplexed per pool with varying numbers of UMIs per cell
# This code was used as an accompanying script to SimulatedPools_UMIperCell.sh to subset UMIs per cell.
# Below is an example of retaining 3000 UMI per cell

########
# MAIN #
########
UMI <- 3000  
file <- read.delim('/home/fiorini9/scratch/mjff/emsemblux_seq_feature_eval/UMI_cell2/24_3000_1/combo_UB_CB.txt', header = F, sep = "\t")
out_dir <- '/home/fiorini9/scratch/mjff/emsemblux_seq_feature_eval/UMI_cell2/24_3000_1'
dim(file) 

file$row <- 1:nrow(file)
V1 <- "none"
V2 <- "none"
row <- "none"
fill_df <- data.frame(V1, V2, row)


for(i in unique(file$V2)){
keep <- i
temp <- subset(file, V2 %in% keep)
lister <- unique(temp$V1)
if(length(lister)< UMI) {
sample <- sample(lister, UMI, replace = T)
}else{
sample <- sample(lister, UMI)
}
temp <- subset(temp, V1 %in% sample)
fill_df <- rbind(fill_df,temp)
}

write.table(fill_df, file = paste(out_dir,"/keep_rows.tsv", sep=""),sep=",", row.names = FALSE, col.names = FALSE, quote = FALSE) 