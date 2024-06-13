# Ensemblex: an accuracy-weighted ensemble genetic demultiplexing framework for single-cell RNA sequencing

[![](https://img.shields.io/badge/Documentation-scrnabox-blue)](https://neurobioinfo.github.io/ensemblex/site/) 
<!-- [![biorxiv](https://img.shields.io/badge/biorxiv-manuscript-blue)](https://www.biorxiv.org/content/biorxiv/early/2023/11/15/2023.11.13.566851.full.pdf)  [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7768427.svg)](https://doi.org/10.5281/zenodo.7768427)  -->

-------------
## Contents
- [Introduction](#introduction)
- [Installation](#installation)
- [Step 1: Set up](#step-1-set-up)
- [Step 2: Preparation of input files](#step-2-preparation-of-input-files)
- [Step 3: Genetic demultiplexing by constituent tools](#step-3-genetic-demultiplexing-by-constituent-tools)
- [Step 4: Application of Ensemblex](#step-4-application-of-ensemblex)
- [Tutorial](#tutorial)
- [Contributing](#contributing)
- [Change log](#change-log)
- [License](#license)
- [Acknowledgement](#acknowledgement)

---
## Introduction
Ensemblex is an accuracy-weighted ensemble framework for genetic demultiplexing of pooled single-cell RNA seqeuncing (scRNAseq) data. Ensemblex can be used to demultiplex pools **with** or **without** prior genotype information. When demultiplexing **with** prior genotype information, Ensemblex leverages the sample assignments of four individual, constituent genetic demultiplexing tools:

1. Demuxalot ([Rogozhnikov et al. ](https://www.biorxiv.org/content/10.1101/2021.05.22.443646v2.abstract))
2. Demuxlet ([Kang et al. ](https://www.nature.com/articles/nbt.4042))
3. Souporcell ([Heaton et al. ](https://www.nature.com/articles/s41592-020-0820-1))
4. Vireo-GT ([Huang et al. ](https://link.springer.com/article/10.1186/s13059-019-1865-2))

When demultiplexing **without** prior genotype information, Ensemblex leverages the sample assignments of four individual, constituent genetic demultiplexing tools:

1. Demuxalot ([Rogozhnikov et al. ](https://www.biorxiv.org/content/10.1101/2021.05.22.443646v2.abstract))
2. Freemuxlet ([Kang et al. ](https://www.nature.com/articles/nbt.4042))
3. Souporcell ([Heaton et al. ](https://www.nature.com/articles/s41592-020-0820-1))
4. Vireo ([Huang et al. ](https://link.springer.com/article/10.1186/s13059-019-1865-2))

Upon demultiplexing pools with each of the four constituent genetic demultiplexing tools, Ensemblex processes the output files in a three-step pipeline to identify the most probable sample label for each cell based on the predictions of the constituent tools:

**Step 1**: Probabilistic-weighted ensemble <br />
**Step 2**: Graph-based doublet detection <br />
**Step 3**: Ensemble-independent doublet detection <br />

 <p align="center">
 <img src="https://github.com/neurobioinfo/ensemblex/assets/97498007/d4161be8-8bc6-4c08-af60-4ee3f486e940" width="650" height="700">
 </p>


As output, Ensemblex returns its own cell-specific sample labels and corresponding assignment probabilities and singlet confidence score, as well as the sample labels and corresponding assignment probabilities for each of its constituents. The demultiplexed sample labels could then be used to perform downstream analyses.

To facilitate the application of Ensemblex, we provide a pipeline that demultiplexes pooled cells by each of the individual constituent genetic demultiplexing tools and processes the outputs with the Ensemblex algorithm. 

The pipelines comprise of four distinct steps:

1. [Selection of Ensemblex pipeline and establishing the working directory (Set up)](#step-1-set-up)
2. [Prepare input files for constituent genetic demultiplexing tools](#step-2-preparation-of-input-files)
3. [Genetic demultiplexing by constituent demultiplexing tools](#step-3-genetic-demultiplexing-by-constituent-tools)
4. [Application of the Ensemblex framework](#step-4-application-of-ensemblex)

 <p align="center">
 <img src="https://github.com/neurobioinfo/ensemblex/assets/97498007/395db4b1-2365-4ca8-8483-097b893dd640" width="550" height="600">
 </p>


Below we provide a quick-start guide for using Ensemblex. For comprehensive documentation, please see the [Ensemblex site](https://neurobioinfo.github.io/ensemblex/site/).
In the Ensemblex documentation, we outline each step of the Ensemblex pipeline, illustrate how to run the pipeline, define best practices, and provide a tutorial with pubicly available datasets. 

---
## Installation

[To be completed.]

---
## Step 1: Set up
### Demultiplexing pooled cells _with_ prior genotype information

Initiate the pipeline:
```
## Create and navigate to the working directory
mkdir working_directory
cd /path/to/working_directory

## Define the path to ensemblex.pip
ensemblex_HOME=/path/to/ensemblex.pip

## Define the path to the working directory
ensemblex_PWD=/path/to/working_directory

## Initiate the pipeline for demultiplexing with prior genotype information
bash $ensemblex_HOME/launch_ensemblex.sh -d $ensemblex_PWD --step init-GT
```

### Demultiplexing pooled cells _without_ prior genotype information

Initiate the pipeline:
```
## Create and navigate to the working directory
mkdir working_directory
cd /path/to/working_directory

## Define the path to ensemblex.pip
ensemblex_HOME=/path/to/ensemblex.pip

## Define the path to the working directory
ensemblex_PWD=/path/to/working_directory

## Initiate the pipeline for demultiplexing without prior genotype information
bash $ensemblex_HOME/launch_ensemblex.sh -d $ensemblex_PWD --step init-noGT
```

---
## Step 2: Preparation of input files
### Demultiplexing pooled cells _with_ prior genotype information

The following files are required:

|File|Description|
|:--|:--|
|**gene_expression.bam**|Gene expression bam file of the pooled samples (e.g., 10X Genomics possorted_genome_bam.bam)|
|**gene_expression.bam.bai**|Gene expression bam index file of the pooled samples (e.g., 10X Genomics possorted_genome_bam.bam.bai)|
|**barcodes.tsv**|Barcodes tsv file of the pooled cells  (e.g., 10X Genomics barcodes.tsv)|
|**pooled_samples.vcf**|vcf file describing the genotypes of the pooled samples|
|**genome_reference.fa**|Genome reference fasta file (e.g., 10X Genomics genome.fa)|
|**genome_reference.fa.fai**|Genome reference fasta index file (e.g., 10X Genomics genome.fa.fai)|
|**genotype_reference.vcf**|Population reference vcf file (e.g., 1000 Genomes Project)|

Prepare the input files:

```
## Define all of the required files
BAM=/path/to/possorted_genome_bam.bam
BAM_INDEX=/path/to/possorted_genome_bam.bam.bai
BARCODES=/path/to/barcodes.tsv
SAMPLE_VCF=/path/to/pooled_samples.vcf
REFERENCE_VCF=/path/to/genotype_reference.vcf
REFERENCE_FASTA=/path/to/genome.fa
REFERENCE_FASTA_INDEX=/path/to/genome.fa.fai

## Define the path to the working directory
ensemblex_PWD=/path/to/working_directory

## Copy the files to the input_files directory in the working directory
cp $BAM  $ensemblex_PWD/input_files/pooled_bam.bam
cp $BAM_INDEX  $ensemblex_PWD/input_files/pooled_bam.bam.bai
cp $BARCODES  $ensemblex_PWD/input_files/pooled_barcodes.tsv
cp $SAMPLE_VCF  $ensemblex_PWD/input_files/pooled_samples.vcf
cp $REFERENCE_VCF  $ensemblex_PWD/input_files/reference.vcf
cp $REFERENCE_FASTA  $ensemblex_PWD/input_files/reference.fa
cp $REFERENCE_FASTA_INDEX  $ensemblex_PWD/input_files/reference.fa.fai
```

### Demultiplexing pooled cells _without_ prior genotype information

The following files are required:

|File|Description|
|:--|:--|
|**gene_expression.bam**|Gene expression bam file of the pooled samples (e.g., 10X Genomics possorted_genome_bam.bam)|
|**gene_expression.bam.bai**|Gene expression bam index file of the pooled samples (e.g., 10X Genomics possorted_genome_bam.bam.bai)|
|**barcodes.tsv**|Barcodes tsv file of the pooled cells  (e.g., 10X Genomics barcodes.tsv)|
|**genome_reference.fa**|Genome reference fasta file (e.g., 10X Genomics genome.fa)|
|**genome_reference.fa.fai**|Genome reference fasta index file (e.g., 10X Genomics genome.fa.fai)|
|**genotype_reference.vcf**|Population reference vcf file (e.g., 1000 Genomes Project)|

Prepare the input files:

```
## Define all of the required files
BAM=/path/to/possorted_genome_bam.bam
BAM_INDEX=/path/to/possorted_genome_bam.bam.bai
BARCODES=/path/to/barcodes.tsv
REFERENCE_VCF=/path/to/genotype_reference.vcf
REFERENCE_FASTA=/path/to/genome.fa
REFERENCE_FASTA_INDEX=/path/to/genome.fa.fai

## Define the path to the working directory
ensemblex_PWD=/path/to/working_directory

## Copy the files to the input_files directory in the working directory
cp $BAM  $ensemblex_PWD/input_files/pooled_bam.bam
cp $BAM_INDEX  $ensemblex_PWD/input_files/pooled_bam.bam.bai
cp $BARCODES  $ensemblex_PWD/input_files/pooled_barcodes.tsv
cp $REFERENCE_VCF  $ensemblex_PWD/input_files/reference.vcf
cp $REFERENCE_FASTA  $ensemblex_PWD/input_files/reference.fa
cp $REFERENCE_FASTA_INDEX  $ensemblex_PWD/input_files/reference.fa.fai
```

---
## Step 3: Genetic demultiplexing by constituent tools
### Demultiplexing pooled cells _with_ prior genotype information

Demultiplex the pooled cells with each of Ensemblex's constituent tools:
```
## Define the paths to Ensemblex and the working directory 
ensemblex_HOME=/path/to/ensemblex.pip
ensemblex_PWD=/path/to/working_directory

## Demuxalot
bash $ensemblex_HOME/launch_ensemblex.sh -d $ensemblex_PWD --step demuxalot

## Demuxlet
bash $ensemblex_HOME/launch_ensemblex.sh -d $ensemblex_PWD --step demuxlet

## Souporcell
bash $ensemblex_HOME/launch_ensemblex.sh -d $ensemblex_PWD --step souporcell

## Vireo
bash $ensemblex_HOME/launch_ensemblex.sh -d $ensemblex_PWD --step vireo

```

### Demultiplexing pooled cells _without_ prior genotype information

Demultiplex the pooled cells with each of Ensemblex's constituent tools:
```
## Define the paths to Ensemblex and the working directory 
ensemblex_HOME=/path/to/ensemblex.pip
ensemblex_PWD=/path/to/working_directory

## Freemuxlet
bash $ensemblex_HOME/launch_ensemblex.sh -d $ensemblex_PWD --step freemuxlet

## Souporcell
bash $ensemblex_HOME/launch_ensemblex.sh -d $ensemblex_PWD --step souporcell

## Vireo
bash $ensemblex_HOME/launch_ensemblex.sh -d $ensemblex_PWD --step vireo

## Demuxalot
bash $ensemblex_HOME/launch_ensemblex.sh -d $ensemblex_PWD --step demuxalot
```

---
## Step 4: Application of Ensemblex

### Demultiplexing pooled cells _with_ prior genotype information
```
## Define the paths to Ensemblex and the working directory 
ensemblex_HOME=/path/to/ensemblex.pip
ensemblex_PWD=/path/to/working_directory

## Compute ensemble classification
bash $ensemblex_HOME/launch_ensemblex.sh -d $ensemblex_PWD --step ensemblexing
```
### Demultiplexing pooled cells _without_ prior genotype information
```
## Define the paths to Ensemblex and the working directory 
ensemblex_HOME=/path/to/ensemblex.pip
ensemblex_PWD=/path/to/working_directory

## Compute ensemble classification
bash $ensemblex_HOME/launch_ensemblex.sh -d $ensemblex_PWD --step ensemblexing
```

---
#### Contributing
Any contributions or suggestions for improving the Ensemblex pipeline are welcomed and appreciated. If you encounter any issues, please open an issue in the [GitHub repository](https://github.com/neurobioinfo/ensemblex). Alternatively, you are welcomed to email the developers directly; for any questions please contact Michael Fiorini: michael.fiorini@mail.mcgill.ca

#### Changelog
Every release is documented on the [GitHub Releases page](https://github.com/neurobioinfo/ensemblex/releases).

#### License
This project is licensed under the MIT License.

#### Acknowledgement
The Ensemblex pipeline was produced for projects funded by the Canadian Institute of Health Research and Michael J. Fox Foundation Parkinson's Progression Markers Initiative (MJFF PPMI) in collaboration with The Neuro's Early Drug Discovery Unit (EDDU), McGill University. It is written by [Michael Fiorini](https://github.com/fiorini9) and [Saeid Amiri](https://github.com/saeidamiri1) with supervision from [Rhalena Thomas](https://github.com/RhalenaThomas) and Sali Farhan at the Montreal Neurological Institute-Hospital. Copyright belongs MNI BIOINFO CORE.



**[⬆ back to top](#contents)**
