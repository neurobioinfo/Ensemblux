# Welcome to the Ensemblux documentation!
Ensemblux is an accuracy-weighted, probabilistic ensemble framework for genetic demultiplexing of population-scale single-cell RNA seqeuncing (scRNAseq) sample pooling. Ensemblux was developed in response to limitations of the individual genetic demultiplexing tools that impede investigators from capitalizing on the full potential of their dataset: 

1. Accurate demultiplexing of ceratin subsets of pooled cells can only be achieved by single tools
2. Correctly classified cells are unecessarily discarded from downstream analyses due to insufficient assignment probabilities
3. Demultiplexing accuracy decreases as the number of pooled samples increases

By addressing the above limitations, we demonstrated that Ensemblux achieves higher demultiplexing accuracy across pools ranging in size from four to 80 multiplexed samples compared to the individual tools, limits the introduction of technical artifacts into scRNAseq analysis, and retains a high proportion of pooled cells for downstream analyses. The ensemble method capitalizes on the added confidence of combining distinct statistical frameworks for genetic demultiplexing, but the modular algorithm can adapt to the overall performance of its constituent tools on the respective dataset; thus, making it resilient against a poorly performing constituent tool.

Ensemblux can be used to demultiplex pools **with** or **without** prior genotype information. When demultiplexing **with** prior genotype information, Ensemblux leverages the sample assignments of four individual, constituent genetic demultiplexing tools:

1. Demuxalot (**citation**)
2. Demuxlet (**citation**)
3. Souporcell (**citation**)
4. Vireo-GT (**citation**)

When demultiplexing **without** prior genotype information, Ensemblux leverages the sample assignments of four individual, constituent genetic demultiplexing tools:

1. Demuxalot (**citation**)
2. Freemuxlet (**citation**)
3. Souporcell (**citation**)
4. Vireo (**citation**)

Upon demultiplexing pools with each of the four constituent genetic demultiplexing tools, Ensemblux processes the output files in a three-step pipeline that was designed to obtain the most probable sample label for each individual cell and capture a high proportion of heterogenic doublets:

**Step 1**: Probabilistic-weighted ensemble <br />
**Step 2**: Graph-based doublet detection <br />
**Step 3**: Ensemble-independent doublet detection <br />

As output, Ensemblux returns its own cell-specific sample labels and corresponding assignment probabilities and singlet confidence score, as well as the sample labels and corresponding assignment probabilities for each of its constituents. The demultiplexed sample labels could then be used to perform downstream analyses.

 <p align="center">
 <img src="https://github.com/mfiorini9/Ensemblux/assets/97498007/b3853b82-231d-43c5-9b00-9f44510a4e84" width="650" height="100">
 </p>

**Figure 1. Overview of the Ensemblux worflow.** **A)** temp **B)** temp  **C)** temp  **D)** temp

To facilitate the application of Ensemblux, we provide a pipeline that demultiplexes pooled cells by each of the individual constituent genetic demultiplexing tools and processes the outputs with the Ensemblux algorithm. In this documentation, we outline each step of the Ensemblux pipeline, illustrate how to run the pipeline, define best practices, and provide a tutorial with pubicly available datasets. 

For a comprehensive descripttion of Ensemblux, ground-truth benchmarking, and application to real-world datasets, see our pre-print manuscript: **Pre-print**

 - - - -

## Contents
- The Ensemblux Framework:
    - [Overview](overview.md)
- The Ensemblux Pipeline:
    - [Overview](overview_pipeline.md)
    - [Installation](installation.md)
    - [Step 1: Set up](Step0.md)
    - [Step 2: FASTQ to expression matrix](Step1.md)
    - [Step 3: Seurat object and ambient RNA](Step2.md)
    - [Step 4: Quality control and filtering](Step3.md)
- Documentation:    
    - [Job configurations](config.md)
    - [Execution parameters](reference.md) 
    - [Outputs](outputs.md) 
- Tutorial:
    - [Downloading data](midbrain_download.md)
    - [Enemblux with prior genotype information](Dataset1.md)
    - [Ensemblux without prior genotype information](pbmc_download.md)              
- About:
    - [License](LICENSE.md)
    - [Help and Feedback](contributing.md)
    - [Acknowledgement](Acknowledgement.md)
