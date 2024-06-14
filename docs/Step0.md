# Step 1: Setting up the Ensemblex pipeline

In Step 1, we will set up the working directory for the Ensemblex pipeline and decide which version of the pipeline we want to use:

1. [Demultiplexing with prior genotype information](#demultiplexing-with-prior-genotype-information)
2. [Demultiplexing without prior genotype information](#demultiplexing-without-prior-genotype-information)

 - - - -

## Demultiplexing with prior genotype information
First, create a dedicated folder for the analysis (hereafter referred to as the working directory). Then, define the path to the working directory and the path to ensemblex.pip:

```
## Create and navigate to the working directory
mkdir working_directory
cd /path/to/working_directory

## Define the path to ensemblex.pip
ensemblex_HOME=/path/to/ensemblex.pip

## Define the path to the working directory
ensemblex_PWD=/path/to/working_directory
```
Next, we can set up the working directory for demultiplexing **with** prior genotype information using the following code:

```
bash $ensemblex_HOME/launch_ensemblex.sh -d $ensemblex_PWD --step init-GT
```
After running the above code, the working directory should have the following structure
```
working_directory
├── demuxalot
├── demuxlet
├── ensemblex_gt
├── input_files
├── job_info
│   ├── configs
│   │   └── ensemblex_config.ini
│   ├── logs
│   └── summary_report.txt
├── souporcell
└── vireo_gt
```

 Upon setting up the Ensemblex pipeline, we can proceed to Step 2 where we will prepare the input files for Ensemblex's constituent genetic demultiplexing tools: [Preparation of input files](Step1.md)

 - - - -
## Demultiplexing without prior genotype information
First, create a dedicated folder for the analysis (hereafter referred to as the working directory). Then, define the path to the working directory and the path to ensemblex.pip:

```
## Create and navigate to the working directory
mkdir working_directory
cd /path/to/working_directory

## Define the path to ensemblex.pip
ensemblex_HOME=/path/to/ensemblex.pip

## Define the path to the working directory
ensemblex_PWD=/path/to/working_directory
```
Next, we can set up the working directory for demultiplexing **without** prior genotype information using the following code:

```
bash $ensemblex_HOME/launch_ensemblex.sh -d $ensemblex_PWD --step init-noGT
```
After running the above code, the working directory should have the following structure
```
working_directory
├── demuxalot
├── freemuxlet
├── ensemblex
├── input_files
├── job_info
│   ├── configs
│   │   └── ensemblex_config.ini
│   ├── logs
│   └── summary_report.txt
├── souporcell
└── vireo
```
 Upon setting up the Ensemblex pipeline, we can proceed to Step 2 where we will prepare the input files for Ensemblex's constituent genetic demultiplexing tools: [Preparation of input files](Step1.md)

