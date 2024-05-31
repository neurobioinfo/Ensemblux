# Ensemblex: an accuracy-weighted ensemble genetic demultiplexing framework for single-cell RNA sequencing

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
 <img src="https://github.com/mfiorini9/Ensemblux/assets/97498007/b3853b82-231d-43c5-9b00-9f44510a4e84" width="650" height="700">
 </p>


As output, Ensemblex returns its own cell-specific sample labels and corresponding assignment probabilities and singlet confidence score, as well as the sample labels and corresponding assignment probabilities for each of its constituents. The demultiplexed sample labels could then be used to perform downstream analyses.

To facilitate the application of Ensemblex, we provide a pipeline that demultiplexes pooled cells by each of the individual constituent genetic demultiplexing tools and processes the outputs with the Ensemblex algorithm. 

The pipelines comprise of four distinct steps:

1. [Selection of Ensemblex pipeline and establishing the working directory (Set up)](#step-1-set-up)
2. [Prepare input files for constituent genetic demultiplexing tools](#step-2-preparation-of-input-files)
3. [Genetic demultiplexing by constituent demultiplexing tools](#step-3-genetic-demultiplexing-by-constituent-tools)
4. [Application of the Ensemblex framework](#step-4-application-of-ensemblex)

 <p align="center">
 <img src="https://github.com/mfiorini9/Ensemblux/assets/97498007/c1ba33da-e5d6-4d2d-82a5-24f46d7e84e0" width="550" height="600">
 </p>


Below we provide a quick-start guide for using Ensemblex. For comprehensive documentation, please see the [Ensemblex site](https://neurobioinfo.github.io/ensemblux/site/).
In the Ensemblex documentation, we outline each step of the Ensemblex pipeline, illustrate how to run the pipeline, define best practices, and provide a tutorial with pubicly available datasets. 

---
## Installation

[To be completed.]

---
## Step 1: Set up
**Demultiplexing pooled cells _with_ prior genotype information**
```
## Create and navigate to the working directory
mkdir working_directory
cd /path/to/working_directory

## Define the path to ensemblux.pip
Ensemblux_HOME=/path/to/ensemblux.pip

## Define the path to the working directory
Ensemblux_PWD=/path/to/working_directory

## Initiate the pipeline for demultiplexing with prior genotype information
bash $Ensemblux_HOME/launch_ensemblux.sh -d $Ensemblux_PWD --step init-GT
```

**Demultiplexing pooled cells _without_ prior genotype information**
```
## Create and navigate to the working directory
mkdir working_directory
cd /path/to/working_directory

## Define the path to ensemblux.pip
Ensemblux_HOME=/path/to/ensemblux.pip

## Define the path to the working directory
Ensemblux_PWD=/path/to/working_directory

## Initiate the pipeline for demultiplexing without prior genotype information
bash $Ensemblux_HOME/launch_ensemblux.sh -d $Ensemblux_PWD --step init-noGT
```

---
## Step 2: Preparation of input files
**Demultiplexing pooled cells _with_ prior genotype information**

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
Ensemblux_PWD=/path/to/working_directory

## Copy the files to the input_files directory in the working directory
cp $BAM  $Ensemblux_PWD/input_files/pooled_bam.bam
cp $BAM_INDEX  $Ensemblux_PWD/input_files/pooled_bam.bam.bai
cp $BARCODES  $Ensemblux_PWD/input_files/pooled_barcodes.tsv
cp $SAMPLE_VCF  $Ensemblux_PWD/input_files/pooled_samples.vcf
cp $REFERENCE_VCF  $Ensemblux_PWD/input_files/reference.vcf
cp $REFERENCE_FASTA  $Ensemblux_PWD/input_files/reference.fa
cp $REFERENCE_FASTA_INDEX  $Ensemblux_PWD/input_files/reference.fa.fai
```

**Demultiplexing pooled cells _without_ prior genotype information**

The following files are required:

|File|Description|
|:--|:--|
|**gene_expression.bam**|Gene expression bam file of the pooled samples (e.g., 10X Genomics possorted_genome_bam.bam)|
|**gene_expression.bam.bai**|Gene expression bam index file of the pooled samples (e.g., 10X Genomics possorted_genome_bam.bam.bai)|
|**barcodes.tsv**|Barcodes tsv file of the pooled cells  (e.g., 10X Genomics barcodes.tsv)|
|**genome_reference.fa**|Genome reference fasta file (e.g., 10X Genomics genome.fa)|
|**genome_reference.fa.fai**|Genome reference fasta index file (e.g., 10X Genomics genome.fa.fai)|
|**genotype_reference.vcf**|Population reference vcf file (e.g., 1000 Genomes Project)|

```
## Define all of the required files
BAM=/path/to/possorted_genome_bam.bam
BAM_INDEX=/path/to/possorted_genome_bam.bam.bai
BARCODES=/path/to/barcodes.tsv
REFERENCE_VCF=/path/to/genotype_reference.vcf
REFERENCE_FASTA=/path/to/genome.fa
REFERENCE_FASTA_INDEX=/path/to/genome.fa.fai

## Define the path to the working directory
Ensemblux_PWD=/path/to/working_directory

## Copy the files to the input_files directory in the working directory
cp $BAM  $Ensemblux_PWD/input_files/pooled_bam.bam
cp $BAM_INDEX  $Ensemblux_PWD/input_files/pooled_bam.bam.bai
cp $BARCODES  $Ensemblux_PWD/input_files/pooled_barcodes.tsv
cp $REFERENCE_VCF  $Ensemblux_PWD/input_files/reference.vcf
cp $REFERENCE_FASTA  $Ensemblux_PWD/input_files/reference.fa
cp $REFERENCE_FASTA_INDEX  $Ensemblux_PWD/input_files/reference.fa.fai
```

---
## Step 3: Genetic demultiplexing by constituent tools

---
## Step 4: Application of Ensemblex

---
#### Contributing
Any contributions or suggestions for improving the scRNAbox pipeline are welcomed and appreciated. 

If you encounter any issues, please open an issue in the [GitHub repository](https://github.com/neurobioinfo/scrnabox).

#### Changelog
Every release is documented on the [GitHub Releases page](https://github.com/neurobioinfo/ensemblux/releases).

#### License
This project is licensed under the MIT License - see the [LICENSE.md](https://github.com/neurobioinfo/ensemblux/blob/main/LICENSE) file for details.

#### Acknowledgement
The Ensemblex pipeline was produced for projects funded by Canadian Institute of Health Research and Michael J. Fox Foundation. The Ensemblex pipeline was produced as part Dark Genome Project. It is written by [Michael Fiorini](https://github.com/fiorini9), [Rhalena Thomas](https://github.com/RhalenaThomas), and [Saeid Amiri](https://github.com/saeidamiri1), with associate of Sali Farhan and Edward Fon at the Montreal Neurological Institute-Hospital. Copyright belongs [MNI BIOINFO CORE](https://github.com/neurobioinfo). 

**[⬆ back to top](#contents)**
