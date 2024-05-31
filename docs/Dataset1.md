---
layout: post
title: Ensemblux tutorial
description: Ensemblux tutorial
date: 2024-05-02
author: Michael Fiorini
published: true
tags: scRNA 
categories: 
comments: false
---
## Ensemblux pipeline with prior genotype information

- [Introduction](#introduction)
- [Installation](#installation)
- [Step 1: Set up](#step-1-set-up)
- [Step 2: Preparation of input files](#step-2-preparation-of-input-files)
- [Step 3: Genetic demultiplexing by constituent tools](#step-3-genetic-demultiplexing-by-constituent-tools)
- [Step 4: Application of Ensemblux](#step-4-application-of-ensemblux)
- [Resource requirements](#resource-requirements)

 - - - -

## Introduction 
This guide illustrates how to use the Ensemblux pipeline to demultiplexed pooled scRNAseq samples with prior genotype information. Here, we will leverage a pooled scRNAseq dataset produced by [Jerber et al.](https://www.nature.com/articles/s41588-021-00801-6). This pool contains induced pluripotent cell lines (iPSC) from 9 healthy controls that were differentiated towards a dopaminergic neuron state. The Ensemblux pipeline is illustrated in the diagram below:

 <p align="center">
 <img src="https://github.com/neurobioinfo/ensemblux/assets/97498007/58dac8ce-c83a-44ef-aa0e-9c31e771388f" width="200" height="100">
 </p>

**NOTE**: To download the necessary files for the tutorial please see the [Downloading data](midbrain_download.md) section of the Ensemblux documentation.

 - - - -

## Installation 
[to be completed]

module load StdEnv/2023
module load apptainer/1.2.4 

 - - - -

## Step 1: Set up
In Step 1, we will set up the working directory for the Ensemblux pipeline and decide which version of the pipeline we want to use.

First, create a dedicated folder for the analysis (hereafter referred to as the working directory). Then, define the path to the working directory and the path to ensemblux.pip:

```
## Create and navigate to the working directory
cd Ensemblux_tutorial
mkdir working_directory
cd ~/Ensemblux_tutorial/working_directory

## Define the path to ensemblux.pip
Ensemblux_HOME=~/ensemblux.pip

## Define the path to the working directory
Ensemblux_PWD=~/Ensemblux_tutorial/working_directory
```

Next, we can set up the working directory and choose the Ensemblux pipeline for demultiplexing with prior genotype information (`--step init-GT`) using the following code:

```
bash $Ensemblux_HOME/launch_ensemblux.sh -d $Ensemblux_PWD --step init-GT
```

After running the above code, the working directory should have the following structure:

```
Ensemblux_tutorial
└── working_directory
    ├── demuxalot
    ├── demuxlet
    ├── ensemblux_gt
    ├── input_files
    ├── job_info
    │   ├── configs
    │   │   └── ensemblux_config.ini
    │   ├── logs
    │   └── summary_report.txt
    ├── souporcell
    └── vireo_gt
```

Upon setting up the Ensemblux pipeline, we can proceed to Step 2 where we will prepare the input files for Ensemblux's constituent genetic demultiplexing tools.

 - - - -

## Step 2: Preparation of input files
In Step 2, we will define the necessary files needed for Ensemblux's constituent genetic demultiplexing tools and will place them within the working directory. 

**Note**: For the tutorial we will be using the data downloaded in the [Downloading data](midbrain_download.md) section of the Ensemblex documentation. 

First, define all of the required files:

```
BAM=~/Ensemblux_tutorial/CellRanger/outs/possorted_genome_bam.bam

BAM_INDEX=~/Ensemblux_tutorial/CellRanger/outs/possorted_genome_bam.bam.bai

BARCODES=~/Ensemblux_tutorial/CellRanger/outs/filtered_gene_bc_matrices/refdata-cellranger-GRCh37/barcodes.tsv

SAMPLE_VCF=~/Ensemblux_tutorial/sample_genotype/sample_genotype_merge.vcf

REFERENCE_VCF=~/Ensemblux_tutorial/reference_files/common_SNPs_only.recode.vcf

REFERENCE_FASTA=~/Ensemblux_tutorial/reference_files/genome.fa

REFERENCE_FASTA_INDEX=~/Ensemblux_tutorial/reference_files/genome.fa.fai
```

Next, we will sort the pooled samples and reference .vcf files according to the .bam file and place them within the working directory:

```
## Sort pooled samples .vcf file
bash $Ensemblux_HOME/launch_ensemblux.sh -d $Ensemblux_PWD/input_files/pooled_samples.vcf --step sort --vcf $SAMPLE_VCF --bam $Ensemblux_PWD/input_files/pooled_bam.bam

## Sort reference .vcf file
bash $Ensemblux_HOME/launch_ensemblux.sh -d $Ensemblux_PWD/input_files/reference.vcf --step sort --vcf $SAMPLE_VCF --bam $Ensemblux_PWD/input_files/pooled_bam.bam
```
Next, we will place the remaining necessary files within the working directory:

```
cp $BAM $Ensemblux_PWD/input_files/pooled_bam.bam
cp $BAM_INDEX $Ensemblux_PWD/input_files/pooled_bam.bam.bai 
cp $BARCODES $Ensemblux_PWD/input_files/pooled_barcodes.tsv
cp $REFERENCE_FASTA $Ensemblux_PWD/input_files/reference.fa
cp $REFERENCE_FASTA_INDEX $Ensemblux_PWD/input_files/reference.fa.fai
```

After running the above code,  `$Ensemblux_PWD/input_files` should contain the following files:
```
input_files
├── pooled_bam.bam
├── pooled_bam.bam.bai
├── pooled_barcodes.tsv
├── pooled_samples.vcf
├── reference.fa
├── reference.fa.fai
└── reference.vcf
```
**NOTE**: It is important that the file names match those listed above as they are necessary for the Ensemblux pipeline to recognize them.

 - - - -

## Step 3: Genetic demultiplexing by constituent tools
In Step 3, we will demultiplex the pooled samples with each of Ensemblux's constituent genetic demultiplexing tools:

- [Demuxalot](#demuxalot)
- [Demuxlet](#demuxlet)
- [Souporcell](#souporcell)
- [Vireo-GT](#vireo-gt)

First, we will navigate to the `ensemblux_config.ini` file to adjust the demultiplexing parameters for each of the constituent genetic demultiplexing tools:

```
## Navigate to the .ini file
cd $Ensemblux_PWD/job_info/configs

## Open the .ini file and adjust parameters directly in the terminal
nano ensemblux_config.ini
```

For the tutorial, we set the following parameters for the constituent genetic demultiplexing tools:

|Parameter|Value|
|:--|:--|
|PAR_demuxalot_genotype_names|'HPSI0115i-hecn_6,HPSI0214i-pelm_3,HPSI0314i-sojd_3,HPSI0414i-sebn_3,HPSI0514i-uenn_3,HPSI0714i-pipw_4,HPSI0715i-meue_5,HPSI0914i-vaka_5,HPSI1014i-quls_2'|
|PAR_demuxalot_prior_strength|100|
|PAR_demuxalot_minimum_coverage|200|
|PAR_demuxalot_minimum_alternative_coverage|10|
|PAR_demuxalot_n_best_snps_per_donor|100|
|PAR_demuxalot_genotypes_prior_strength|1|
|PAR_demuxalot_doublet_prior|0.25|
|PAR_demuxlet_field|GT|
|PAR_vireo_N|9|
|PAR_vireo_type|GT|
|PAR_vireo_processes|20|
|PAR_vireo_minMAF|0.1|
|PAR_vireo_minCOUNT|20|
|PAR_vireo_forcelearnGT|T|
|PAR_minimap2|'-ax splice -t 8 -G50k -k 21 -w 11 --sr -A2 -B8 -O12,32 -E2,1 -r200 -p.5 -N20 -f1000,5000 -n2 -m20 -s40 -g2000 -2K50m --secondary=no'|
|PAR_freebayes|'-iXu -C 2 -q 20 -n 3 -E 1 -m 30 --min-coverage 6'|
|PAR_vartrix_umi|TRUE|
|PAR_vartrix_mapq|30|
|PAR_vartrix_threads|8|
|PAR_souporcell_k|9|
|PAR_souporcell_t|8|


Now that the parameters have been defined, we can demultiplex the pools with the constituent genetic demultiplexing tools.

 - - - -

#### Demuxalot
To run Demuxalot use the following code:

```
bash $Ensemblux_HOME/launch_ensemblux.sh -d $Ensemblux_PWD --step demuxalot
```

If Demuxalot completed successfully, the following files should be available in `$Ensemblux_PWD/demuxalot`:

```
demuxalot
    ├── Demuxalot_result.csv
    └── new_snps_single_file.betas
```

  - - - -

#### Demuxlet

To run Demuxlet use the following code:

```
bash $Ensemblux_HOME/launch_ensemblux.sh -d $Ensemblux_PWD --step demuxlet
```

If Demuxlet completed successfully, the following files should be available in `$Ensemblux_PWD/demuxlet`:

```
demuxlet
    ├── outs.best
    ├── pileup.cel.gz
    ├── pileup.plp.gz
    ├── pileup.umi.gz
    └── pileup.var.gz
```

  - - - -

#### Souporcell

To run Souporcell use the following code:

```
bash $Ensemblux_HOME/launch_ensemblux.sh -d $Ensemblux_PWD --step souporcell
```

If Souporcell completed successfully, the following files should be available in `$Ensemblux_PWD/souporcell`:

```
souporcell
    ├── alt.mtx
    ├── cluster_genotypes.vcf
    ├── clusters_tmp.tsv
    ├── clusters.tsv
    ├── fq.fq
    ├── minimap.sam
    ├── minitagged.bam
    ├── minitagged_sorted.bam
    ├── minitagged_sorted.bam.bai
    ├── Pool.vcf
    ├── ref.mtx
    └── soup.txt
```
  - - - -

#### Vireo

To run Vireo-GT use the following code:

```
bash $Ensemblux_HOME/launch_ensemblux.sh -d $Ensemblux_PWD --step vireo
```

If Vireo-GT completed successfully, the following files should be available in `$Ensemblux_PWD/vireo_gt`:

```
vireo_gt
    ├── cellSNP.base.vcf.gz
    ├── cellSNP.cells.vcf.gz
    ├── cellSNP.samples.tsv
    ├── cellSNP.tag.AD.mtx
    ├── cellSNP.tag.DP.mtx
    ├── cellSNP.tag.OTH.mtx
    ├── donor_ids.tsv
    ├── fig_GT_distance_estimated.pdf
    ├── fig_GT_distance_input.pdf
    ├── GT_donors.vireo.vcf.gz
    ├── _log.txt
    ├── prob_doublet.tsv.gz
    ├── prob_singlet.tsv.gz
    └── summary.tsv
```
  - - - -

Upon demultiplexing the pooled samples with each of Ensemblux's constituent genetic demultiplexing tools, we can proceed to Step 4 where we will process the output files of the consituent tools with the Ensemblux algorithm to generate the ensemble sample classifications

**NOTE**: To minimize computation time for the tutorial, we have provided the necessary outpu files from the constituent tools [here](https://github.com/neurobioinfo/ensemblux/tree/main/tutorial). To access the files and place them in the working directory, use the following code:

```
## Demuxalot
cd $Ensemblux_PWD/demuxalot
wget https://github.com/neurobioinfo/ensemblux/blob/caad8c250566bfa9a6d7a78b77d2cc338468a58e/tutorial/Demuxalot_result.csv

## Demuxlet
cd $Ensemblux_PWD/demuxlet
wget https://github.com/neurobioinfo/ensemblux/blob/caad8c250566bfa9a6d7a78b77d2cc338468a58e/tutorial/outs.best

## Souporcell
cd $Ensemblux_PWD/souporcell
wget https://github.com/neurobioinfo/ensemblux/blob/caad8c250566bfa9a6d7a78b77d2cc338468a58e/tutorial/clusters.tsv

## Vireo
cd $Ensemblux_PWD/vireo_gt
wget https://github.com/neurobioinfo/ensemblux/blob/caad8c250566bfa9a6d7a78b77d2cc338468a58e/tutorial/donor_ids.tsv

```

 - - - -
## Step 4: Application of Ensemblux
In Step 4, we will process the output files of the four constituent genetic demultiplexing tools with the three-step Ensemblex algorithm:

- Step 1: Probabilistic-weighted ensemble
- Step 2: Graph-based doublet detection
- Step 3: Step 3: Ensemble-independent doublet detection

First, we will navigate to the `ensemblux_config.ini` file to adjust the demultiplexing parameters for the Ensemblex algorithm:

```
## Navigate to the .ini file
cd $Ensemblux_PWD/job_info/configs

## Open the .ini file and adjust parameters directly in the terminal
nano ensemblux_config.ini
```

For the tutorial, we set the following parameters for the Ensemblex algorithm:

|Parameter|Value|
|:--|:--|
|<ins>**Pool parameters**</ins>|
|PAR_ensemblux_sample_size| 9|
|PAR_ensemblux_expected_doublet_rate| 0.10|
|<ins>**Set up parameters**</ins>|
|PAR_ensemblux_merge_constituents|Yes|
|<ins>**Step 1 parameters: Probabilistic-weighted ensemble**</ins>|
|PAR_ensemblux_probabilistic_weighted_ensemble| Yes|
|<ins>**Step 2 parameters: Graph-based doublet detection**</ins>|
|PAR_ensemblux_preliminary_parameter_sweep| No|
|PAR_ensemblux_nCD| NULL|
|PAR_ensemblux_pT| NULL|
|PAR_ensemblux_graph_based_doublet_detection| Yes|
|<ins>**Step 3 parameters: Ensemble-independent doublet detection**</ins>|
|PAR_ensemblux_preliminary_ensemble_independent_doublet| No|
|PAR_ensemblux_ensemble_independent_doublet| Yes|
|PAR_ensemblux_doublet_Demuxalot_threshold| Yes|
|PAR_ensemblux_doublet_Demuxalot_no_threshold| No|
|PAR_ensemblux_doublet_Demuxlet_threshold| No|
|PAR_ensemblux_doublet_Demuxlet_no_threshold| No|
|PAR_ensemblux_doublet_Souporcell_threshold| No|
|PAR_ensemblux_doublet_Souporcell_no_threshold| No|
|PAR_ensemblux_doublet_Vireo_threshold| Yes|
|PAR_ensemblux_doublet_Vireo_no_threshold| No|
|<ins>**Confidence score parameters**</ins>|
|PAR_ensemblux_compute_singlet_confidence| Yes|

If Ensemblex completed successfully, the following files should be available in `$Ensemblux_PWD/ensemblux_gt`:

```
ensemblux_gt
├── confidence
│   └── Ensemblux_final_cell_assignment.csv
├── constituent_tool_merge.csv
├── step1
│   ├── ARI_demultiplexing_tools.pdf
│   ├── BA_demultiplexing_tools.pdf
│   ├── Balanced_accuracy_summary.csv
│   └── step1_cell_assignment.csv
├── step2
│   ├── optimal_nCD.pdf
│   ├── optimal_pT.pdf
│   ├── PC1_var_contrib.pdf
│   ├── PC2_var_contrib.pdf
│   ├── PCA1_graph_based_doublet_detection.pdf
│   ├── PCA2_graph_based_doublet_detection.pdf
│   ├── PCA3_graph_based_doublet_detection.pdf
│   ├── PCA_plot.pdf
│   ├── PCA_scree_plot.pdf
│   └── Step2_cell_assignment.csv
└── step3
    ├── Doublet_overlap_no_threshold.pdf
    ├── Doublet_overlap_threshold.pdf
    ├── Number_Ensemblux_doublets_EID_no_threshold.pdf
    ├── Number_Ensemblux_doublets_EID_threshold.pdf
    └── Step3_cell_assignment.csv
```

Ensemblex's final assignments are described in the `Ensemblux_final_cell_assignment.csv` file. 

Specifically, the `Ensemblux_assignment` column describes Ensemblex's final assignments after application of the singlet confidence threshold (i.e., singlets that fail to meet a singlet confidence of 1.0 are labelled as unassigned); we recomment that users use this column to label their cells for downstream analyses. The `Ensemblux_best_assignment` column describes Ensemblex's best assignments, independent of the singlets confidence threshold (i.e., singlets that fail to meet a singlet confidence of 1.0 are **NOT** labelled as unassigned).

The cell barcodes listed under the `barcode` column can be used to add the `Ensemblux_final_cell_assignment.csv` information to the metadata of a Seurat object. 

 - - - -
## Resource requirements
The following table describes the computational resources used in this tutorial for genetic demultiplexing by the constituent tools and application of the Ensemblex algorithm.

|Tool|Time|CPU|Memory|
|:--|:--|:--|:--|
|Demuxalot|01:34:59|6|12.95 GB|
|Demuxlet|03:16:03|6|138.32 GB|
|Souporcell|2-14:49:21|1|21.83 GB|
|Vireo|2-01:30:24|6|29.42 GB|
|Ensemblex|02:05:27|1|5.67 GB|

 - - - -
