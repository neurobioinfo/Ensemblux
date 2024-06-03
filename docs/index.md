# Welcome to the ensemblex documentation!
ensemblex is an accuracy-weighted ensemble framework for genetic demultiplexing of pooled single-cell RNA seqeuncing (scRNAseq) data. By addressing the limitiations of individual genetic demultiplexing tools, we demonstrated that ensemblex:

- Achieves higher demultiplexing accuracy
- Limits the introduction of technical noise into scRNAseq analysis 
- Retains a high proportion of cells for downstream analyses. 

The ensemble method capitalizes on the added confidence of combining distinct statistical frameworks for genetic demultiplexing, but the modular algorithm can adapt to the overall performance of its constituent tools on the respective dataset, making it resilient against a poorly performing constituent tool. 

ensemblex can be used to demultiplex pools **with** or **without** prior genotype information. When demultiplexing **with** prior genotype information, ensemblex leverages the sample assignments of four individual, constituent genetic demultiplexing tools:

1. Demuxalot ([Rogozhnikov et al. ](https://www.biorxiv.org/content/10.1101/2021.05.22.443646v2.abstract))
2. Demuxlet ([Kang et al. ](https://www.nature.com/articles/nbt.4042))
3. Souporcell ([Heaton et al. ](https://www.nature.com/articles/s41592-020-0820-1))
4. Vireo-GT ([Huang et al. ](https://link.springer.com/article/10.1186/s13059-019-1865-2))

When demultiplexing **without** prior genotype information, ensemblex leverages the sample assignments of four individual, constituent genetic demultiplexing tools:

1. Demuxalot ([Rogozhnikov et al. ](https://www.biorxiv.org/content/10.1101/2021.05.22.443646v2.abstract))
2. Freemuxlet ([Kang et al. ](https://www.nature.com/articles/nbt.4042))
3. Souporcell ([Heaton et al. ](https://www.nature.com/articles/s41592-020-0820-1))
4. Vireo ([Huang et al. ](https://link.springer.com/article/10.1186/s13059-019-1865-2))

Upon demultiplexing pools with each of the four constituent genetic demultiplexing tools, ensemblex processes the output files in a three-step pipeline to identify the most probable sample label for each cell based on the predictions of the constituent tools:

**Step 1**: Probabilistic-weighted ensemble <br />
**Step 2**: Graph-based doublet detection <br />
**Step 3**: Ensemble-independent doublet detection <br />

As output, ensemblex returns its own cell-specific sample labels and corresponding assignment probabilities and singlet confidence score, as well as the sample labels and corresponding assignment probabilities for each of its constituents. The demultiplexed sample labels could then be used to perform downstream analyses.

 <p align="center">
 <img src="https://github.com/neurobioinfo/ensemblex/assets/97498007/d4161be8-8bc6-4c08-af60-4ee3f486e940" width="650" height="100">
 </p>

**Figure 1. Overview of the ensemblex worflow.** **A)** The Ensemblex workflow begins with demultiplexing pooled samples by each of the constituent tools. The outputs from each individual demultiplexing tool are then used as input into the Ensemblex framework.  **B)** The Ensemblex framework comprises three distinct steps that are assembled into a pipeline: 1) accuracy-weighted probabilistic ensemble, 2) graph-based doublet detection, and 3) ensemble-independent doublet detection.  **C)** As output, ensemblex returns its own sample-cell assignments as well as the sample-cell assignments of each of its constituent tools.  **D)** ensemblex's sample-cell assignments can be used to perform downstream analysis on the pooled scRNAseq data. 

To facilitate the application of ensemblex, we provide a pipeline that demultiplexes pooled cells by each of the individual constituent genetic demultiplexing tools and processes the outputs with the ensemblex algorithm. In this documentation, we outline each step of the ensemblex pipeline, illustrate how to run the pipeline, define best practices, and provide a tutorial with pubicly available datasets. 

For a comprehensive descripttion of ensemblex, ground-truth benchmarking, and application to real-world datasets, see our pre-print manuscript: **Pre-print**

 - - - -

## Contents
- The ensemblex Algorithm:
    - [ensemblex algorithm overview](overview.md)
- The ensemblex Pipeline:
    - [ensemblex pipeline overview](overview_pipeline.md)
    - [Installation](installation.md)
    - [Step 1: Set up](Step0.md)
    - [Step 2: FASTQ to expression matrix](Step1.md)
    - [Step 3: Seurat object and ambient RNA](Step2.md)
    - [Step 4: Quality control and filtering](Step3.md)
- Documentation:    
    - [Execution parameters](reference.md) 
    - [ensemblex outputs](outputs.md) 
- Tutorial:
    - [Downloading data](midbrain_download.md)
    - [Enemblux with prior genotype information](Dataset1.md)            
- About:
    - [Help and Feedback](contributing.md)
    - [Acknowledgement](Acknowledgement.md)
    - [License](LICENSE.md)