## Data Download
- [Introduction](#introduction)
- [Downloading and processing scRNAseq data](#downloading-and-processing-scrnaseq-data)
- [Downloading sample genotype data](#downloading-sample-genotype-data)
- [Downloading reference genotype data](#downloading-reference-genotype-data)
- [Downloading genome reference file](#downloading-genome-reference-file)

 - - - -

## Introduction

For the tutorial, we will leverage a pooled scRNAseq dataset produced by [Jerber et al.](https://www.nature.com/articles/s41588-021-00801-6). This pool contains induced pluripotent cell lines (iPSC) from 9 healthy controls that were differentiated towards a dopaminergic neuron state. 

In this section of the tutorial, we will:

1. Download and process the pooled scRNAseq data with the CellRanger *counts* pipeline
2. Download and process the sample genotype data
3. Download reference genotype data
4. Download a reference genome file

Before we begin, we will create a designated folder for the Ensemblex tutorial:

```
mkdir ensemblex_tutorial
cd ensemblex_tutorial
```

 - - - -

## Downloading and processing scRNAseq data

We will begin by downloading the pooled scRNAseq data from the Sequence Read Archive (SRA):

```
## Create a folder to place pooled scRNAseq data
mkdir pooled_scRNAseq
cd pooled_scRNAseq


## Download pooled scRNAseq FASTQ files
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR470/009/ERR4700019/ERR4700019_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR470/009/ERR4700019/ERR4700019_2.fastq.gz

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR470/000/ERR4700020/ERR4700020_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR470/000/ERR4700020/ERR4700020_2.fastq.gz

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR470/001/ERR4700021/ERR4700021_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR470/001/ERR4700021/ERR4700021_2.fastq.gz

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR470/002/ERR4700022/ERR4700022_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR470/002/ERR4700022/ERR4700022_2.fastq.gz


## Rename pooled scRNAseq FASTQ files
mv ERR4700019_1.fastq.gz ~/ensemblex_tutorial/pooled_scRNAseq/pool_S1_L001_R1_001.fastq.gz 
mv ERR4700019_2.fastq.gz ~/ensemblex_tutorial/pooled_scRNAseq/pool_S1_L001_R2_001.fastq.gz 

mv ERR4700020_1.fastq.gz ~/ensemblex_tutorial/pooled_scRNAseq/pool_S1_L002_R1_001.fastq.gz 
mv ERR4700020_2.fastq.gz ~/ensemblex_tutorial/pooled_scRNAseq/pool_S1_L002_R2_001.fastq.gz

mv ERR4700021_1.fastq.gz ~/ensemblex_tutorial/pooled_scRNAseq/pool_S1_L003_R1_001.fastq.gz
mv ERR4700021_2.fastq.gz ~/ensemblex_tutorial/pooled_scRNAseq/pool_S1_L003_R2_001.fastq.gz
    
mv ERR4700022_1.fastq.gz ~/ensemblex_tutorial/pooled_scRNAseq/pool_S1_L004_R1_001.fastq.gz
mv ERR4700022_2.fastq.gz ~/ensemblex_tutorial/pooled_scRNAseq/pool_S1_L004_R2_001.fastq.gz

```
Next, we will process the pooled scRNAseq data with the CellRanger *counts* pipeline:

```
## Create CellRanger directory
cd ~/ensemblex_tutorial
mkdir CellRanger
cd CellRanger

cellranger count \
--id=pool \
--fastqs=/home/fiorini9/scratch/ensemblex_pipeline_test/ensemblex_tutorial/pooled_scRNAseq \
--sample=pool \
--transcriptome=~/10xGenomics/refdata-cellranger-GRCh37
```

If the CellRanger *counts* pipeline completed successfully, it will have generated the following files that we will use for genetic demultiplexing downstream:

- possorted_genome_bam.bam
- possorted_genome_bam.bam.bai
- barcodes.tsv

**NOTE**: For more information regarding the CellRanger *counts* pipeline, please see the [10X documentation](https://www.10xgenomics.com/support/software/cell-ranger/latest/tutorials/cr-tutorial-ct).

 - - - -

## Downloading sample genotype data

Next, we will download the whole exome .vcf files corresponding to the nine pooled individuals from which the iPSC lines derived. We will download the .vcf files from the European Nucleotide Archive (ENA):

```
## Create a folder to place sample genotype data
cd ~/ensemblex_tutorial
mkdir sample_genotype
cd sample_genotype

## HPSI0115i-hecn_6
wget ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ487/ERZ487971/HPSI0115i-hecn_6.wes.exomeseq.SureSelect_HumanAllExon_v5.mpileup.20170327.genotypes.vcf.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ487/ERZ487971/HPSI0115i-hecn_6.wes.exomeseq.SureSelect_HumanAllExon_v5.mpileup.20170327.genotypes.vcf.gz.tbi

## HPSI0214i-pelm_3 
wget ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ122/ERZ122924/HPSI0214i-pelm_3.wes.exomeseq.SureSelect_HumanAllExon_v5.mpileup.20150415.genotypes.vcf.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ122/ERZ122924/HPSI0214i-pelm_3.wes.exomeseq.SureSelect_HumanAllExon_v5.mpileup.20150415.genotypes.vcf.gz.tbi

## HPSI0314i-sojd_3 
wget ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ266/ERZ266723/HPSI0314i-sojd_3.wes.exomeseq.SureSelect_HumanAllExon_v5.mpileup.20160122.genotypes.vcf.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ266/ERZ266723/HPSI0314i-sojd_3.wes.exomeseq.SureSelect_HumanAllExon_v5.mpileup.20160122.genotypes.vcf.gz.tbi

## HPSI0414i-sebn_3 
wget ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ376/ERZ376769/HPSI0414i-sebn_3.wes.exomeseq.SureSelect_HumanAllExon_v5.mpileup.20161031.genotypes.vcf.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ376/ERZ376769/HPSI0414i-sebn_3.wes.exomeseq.SureSelect_HumanAllExon_v5.mpileup.20161031.genotypes.vcf.gz.tbi

## HPSI0514i-uenn_3 
wget ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ488/ERZ488039/HPSI0514i-uenn_3.wes.exomeseq.SureSelect_HumanAllExon_v5.mpileup.20170327.genotypes.vcf.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ488/ERZ488039/HPSI0514i-uenn_3.wes.exomeseq.SureSelect_HumanAllExon_v5.mpileup.20170327.genotypes.vcf.gz.tbi

## HPSI0714i-pipw_4 
wget ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ376/ERZ376869/HPSI0714i-pipw_4.wes.exomeseq.SureSelect_HumanAllExon_v5.mpileup.20161031.genotypes.vcf.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ376/ERZ376869/HPSI0714i-pipw_4.wes.exomeseq.SureSelect_HumanAllExon_v5.mpileup.20161031.genotypes.vcf.gz.tbi

## HPSI0715i-meue_5 
wget ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ376/ERZ376787/HPSI0715i-meue_5.wes.exomeseq.SureSelect_HumanAllExon_v5.mpileup.20161031.genotypes.vcf.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ376/ERZ376787/HPSI0715i-meue_5.wes.exomeseq.SureSelect_HumanAllExon_v5.mpileup.20161031.genotypes.vcf.gz.tbi

## HPSI0914i-vaka_5 
wget ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ487/ERZ487965/HPSI0914i-vaka_5.wes.exomeseq.SureSelect_HumanAllExon_v5.mpileup.20170327.genotypes.vcf.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ487/ERZ487965/HPSI0914i-vaka_5.wes.exomeseq.SureSelect_HumanAllExon_v5.mpileup.20170327.genotypes.vcf.gz.tbi

## HPSI1014i-quls_2
wget ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ487/ERZ487886/HPSI1014i-quls_2.wes.exomeseq.SureSelect_HumanAllExon_v5.mpileup.20170327.genotypes.vcf.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ487/ERZ487886/HPSI1014i-quls_2.wes.exomeseq.SureSelect_HumanAllExon_v5.mpileup.20170327.genotypes.vcf.gz.tbi
```

Upon downloading the individual genotype data, we will merge the individual files to generate a single .vcf file. 

```
## Merge .vcf files
module load bcftools
bcftools merge *.vcf.gz > sample_genotype_merge.vcf
```

The resulting `sample_genotype_merge.vcf` file will be used as prior genotype information for genetic demultiplexing downstream.

 - - - -

## Downloading reference genotype data
Next, we will download a reference genotype file from the [1000 Genomes Project, Phase 3](https://www.internationalgenome.org/):

```
## Create a folder to place the reference files
cd ~/ensemblex_tutorial
mkdir reference_files
cd reference_files

## Download reference .vcf
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.vcf.gz
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.vcf.gz.tbi

## Unzip .vcf file
gunzip ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.vcf.gz

## Only keep SNPs
module load vcftools
vcftools --vcf ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.vcf --remove-indels --recode --recode-INFO-all --out SNPs_only


## Only keep common variants
module load bcftools
bcftools filter -e 'AF<0.01' SNPs_only.recode.vcf > common_SNPs_only.recode.vcf
```

The resulting `common_SNPs_only.recode.vcf` file will be used as reference genotype data for genetic demultiplexing downstream.

 - - - -

## Downloading genome reference file
Finally, we will prepare a reference genome. For our tutorial we will use the GRCh37 10X reference genome. For information regarding references, see the [10X documentation](https://www.10xgenomics.com/support/software/cell-ranger/downloads/cr-ref-build-steps).

```
## Copy pre-built reference genome to working directory
cp /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/10xGenomics/refdata-cellranger-GRCh37/fasta/genome.fa ~/ensemblex_pipeline_test/ensemblex_tutorial/reference_files
```

We will use the `genome.fa` reference genome for genetic demultiplexing downstream.

 - - - -

 To run the Ensemblex pipeline on the downloaded data please see the [Ensemblex with prior genotype information](Dataset1.md) section of the Ensemblex pipeline.
