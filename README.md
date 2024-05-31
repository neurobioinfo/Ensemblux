# Ensemblex: an accuracy-weighted ensemble genetic demultiplexing framework for single-cell RNA sequencing

-------------
## Contents
- [Introduction](#introduction)
- [Installation](#installation)
- [Pipeline steps](#pipeline-steps)
- [Running Ensemblex](#running-scrnabox)
- [Tutorial](#tutorial)

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

1. [Selection of Ensemblex pipeline and establishing the working directory (Set up)](Step0.md)
2. [Prepare input files for constituent genetic demultiplexing tools](Step1.md)
3. [Genetic demultiplexing by constituent demultiplexing tools](Step2.md)
4. [Application of the Ensemblex framework](Step3.md)

 <p align="center">
 <img src="https://github.com/mfiorini9/Ensemblux/assets/97498007/c1ba33da-e5d6-4d2d-82a5-24f46d7e84e0" width="550" height="100">
 </p>


Below we provide a quick-start guide for using Ensemblex. For more detailed documentation, please see the [Ensemblex site](https://neurobioinfo.github.io/ensemblux/site/).
In this documentation, we outline each step of the Ensemblex pipeline, illustrate how to run the pipeline, define best practices, and provide a tutorial with pubicly available datasets. 

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
