# Step 2: Preparing input files for genetic demultiplexing
In Step 2, we will define the necessary files needed for Ensemblux's constituent genetic demultiplexing tools and will place them within the working directory. The necessary files vary depending on the version of the Ensemblux pipeline being used:

- [Demultiplexing with prior genotype information](#demultiplexing-with-prior-genotype-information)
- [Demultiplexing without prior genotype information](#demultiplexing-without-prior-genotype-information)

 - - - -
## Demultiplexing with prior genotype information
#### Required files
To demultiplex the pooled samples **with** prior genotype information, the following files are required:

|File|Description|
|:--|:--|
|**gene_expression.bam**|Gene expression bam file of the pooled samples (e.g., 10X Genomics possorted_genome_bam.bam)|
|**gene_expression.bam.bai**|Gene expression bam index file of the pooled samples (e.g., 10X Genomics possorted_genome_bam.bam.bai)|
|**barcodes.tsv**|Barcodes tsv file of the pooled cells  (e.g., 10X Genomics barcodes.tsv)|
|**pooled_samples.vcf**|vcf file describing the genotypes of the pooled samples|
|**genome_reference.fa**|Genome reference fasta file (e.g., 10X Genomics: /Homo_sapiens.GRCh37/genome/10xGenomics/refdata-cellranger-GRCh37/fasta/genome.fa)|
|**genome_reference.fa.fai**|Genome reference fasta index file (e.g., 10X Genomics: ~/Homo_sapiens.GRCh37/genome/10xGenomics/refdata-cellranger-GRCh37/fasta/genome.fa.fai)|
|**genotype_reference.vcf**|Population reference vcf file (e.g., 1000 Genomes Project)|

**NOTE:** We demonstrate how to download reference vcf and fasta files in the [Tutorial](midbrain_download.md) section of the Ensemblux documentation. 

#### Placing files into the Ensemblux pipeline working directory
First, define all of the required files:
```
BAM=/path/to/possorted_genome_bam.bam
BAM_INDEX=/path/to/possorted_genome_bam.bam.bai
BARCODES=/path/to/barcodes.tsv
SAMPLE_VCF=/path/to/pooled_samples.vcf
REFERENCE_VCF=/path/to/genotype_reference.vcf
REFERENCE_FASTA=/path/to/genome.fa
REFERENCE_FASTA_INDEX=/path/to/genome.fa.fai
```
Then, place the required files in the Ensemblux pipeline working directory:

```
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

If the file transfer was successful, the input_files directory of the Ensemblux pipeline working directory will contain the following files:
```
working_directory
└── input_files
    ├── pooled_bam.bam
    ├── pooled_bam.bam.bai
    ├── pooled_barcodes.tsv
    ├── pooled_samples.vcf
    ├── reference.fa
    ├── reference.fa.fai
    └── reference.vcf
```
**NOTE:** You will notice that the names of the input files have been standardized, it is important that the input files have the corresonding name for the Ensemblux pipeline to work properly. 

 Upon placing the required files in the Ensemblux pipeline, we can proceed to Step 3 where we will demultiplex the pooled samples using Ensemblux's constituent genetic demultiplexing tools: [Genetic demultiplexing by consituent tools](Step2.md)
 - - - -
 
## Demultiplexing without prior genotype information
#### Required files
To demultiplex the pooled samples **without** prior genotype information, the following files are required:

|File|Description|
|:--|:--|
|**gene_expression.bam**|Gene expression bam file of the pooled samples (e.g., 10X Genomics possorted_genome_bam.bam)|
|**gene_expression.bam.bai**|Gene expression bam index file of the pooled samples (e.g., 10X Genomics possorted_genome_bam.bam.bai)|
|**barcodes.tsv**|Barcodes tsv file of the pooled cells  (e.g., 10X Genomics barcodes.tsv)|
|**genome_reference.fa**|Genome reference fasta file (e.g., 10X Genomics: /Homo_sapiens.GRCh37/genome/10xGenomics/refdata-cellranger-GRCh37/fasta/genome.fa)|
|**genome_reference.fa.fai**|Genome reference fasta index file (e.g., 10X Genomics: ~/Homo_sapiens.GRCh37/genome/10xGenomics/refdata-cellranger-GRCh37/fasta/genome.fa.fai)|
|**genotype_reference.vcf**|Population reference vcf file (e.g., 1000 Genomes Project)|

**NOTE:** We demonstrate how to download reference vcf and fasta files in the [Tutorial](midbrain_download.md) section of the Ensemblux documentation. 

#### Placing files into the Ensemblux pipeline working directory
First, define all of the required files:
```
BAM=/path/to/possorted_genome_bam.bam
BAM_INDEX=/path/to/possorted_genome_bam.bam.bai
BARCODES=/path/to/barcodes.tsv
REFERENCE_VCF=/path/to/genotype_reference.vcf
REFERENCE_FASTA=/path/to/genome.fa
REFERENCE_FASTA_INDEX=/path/to/genome.fa.fai
```
Then, place the required files in the Ensemblux pipeline working directory:

```
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

If the file transfer was successful, the input_files directory of the Ensemblux pipeline working directory will contain the following files:
```
working_directory
└── input_files
    ├── pooled_bam.bam
    ├── pooled_bam.bam.bai
    ├── pooled_barcodes.tsv
    ├── reference.fa
    ├── reference.fa.fai
    └── reference.vcf
```
**NOTE:** You will notice that the names of the input files have been standardized, it is important that the input files have the corresonding name for the Ensemblux pipeline to work properly. 

 Upon placing the required files in the Ensemblux pipeline, we can proceed to Step 3 where we will demultiplex the pooled samples using Ensemblux's constituent genetic demultiplexing tools: [Genetic demultiplexing by consituent tools](Step2.md)