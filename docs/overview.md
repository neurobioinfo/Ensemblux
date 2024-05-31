# Ensemblux algorithm overview

- [Workflow](#workflow)
- [Step 1: Accuracy-weighted probabilistic ensemble](#step-1-accuracy-weighted-probabilistic-ensemble)
- [Step 2: Graph-based doublet detection](#step-2-graph-based-doublet-detection)
- [Step 3: Ensemble-independent doublet detection](#step-3-ensemble-independent-doublet-detection)
- [Contribution of each step to overall demultiplexing accuracy](#contribution-of-each-step-to-overall-demultiplexing-accuracy)

 - - - -
### Workflow

The Ensemblux workflow begins by demultiplexing pooled cells with each of its constituent tools: Demuxalot, Demuxlet, Souporcell and Vireo-GT if using prior genotype information or Demuxalot, Freemuxlet, Souporcell and Vireo if prior genotype information is not available.

 <p align="center">
 <img src="https://github.com/mfiorini9/Ensemblux/assets/97498007/651e984d-da37-46e0-91b6-421aaa64dae9" width="650" height="100">
 </p>

 **Figure 1. Input into the Ensemblux framework.**  The Ensemblex workflow begins with demultiplexing pooled samples by each of the constituent tools. The outputs from each individual demultiplexing tool are then used as input into the Ensemblex framework. 

Upon demultiplexing pools with each individual constituent genetic demultiplexing tool, Ensemblux processes the outputs in a three-step pipeline:

- [Step 1: Accuracy-weighted probabilistic-weighted ensemble](#step-1-accuracy-weighted-probabilistic-ensemble)
- [Step 2: Graph-based doublet detection](#step-2-graph-based-doublet-detection)
- [Step 3: Ensemble-independent doublet detection](#step-3-ensemble-independent-doublet-detection)

 <p align="center">
 <img src="https://github.com/mfiorini9/Ensemblux/assets/97498007/0b3dc87f-e95a-4ee7-91b3-488ac818d3ac" width="650" height="100">
 </p>

**Figure 2. Overview of the three-step Ensemblux framework.**  The Ensemblex framework comprises three distinct steps that are assembled into a pipeline: 1) accuracy-weighted probabilistic ensemble, 2) graph-based doublet detection, and 3) ensemble-independent doublet detection.

For demonstration purposes throughout this section, we leveraged simulated pools with known ground-truth sample labels that were generated with 80 independetly-sequenced induced pluripotent stem cell (iPSC) lines from individuals with Parkinson's disease and neurologically healthy controls. The lines were differentiated towards a dopaminergic cell fate as part of the Foundational Data Initiative for Parkinson's disease (FOUNDIN-PD; [Bressan et al.](https://www.cell.com/cell-genomics/pdf/S2666-979X(23)00017-4.pdf))

 - - - -

### Step 1: Accuracy-weighted probabilistic ensemble
The accuracy-weighted probabilistic ensemble component of the Ensemblux  utilizes an unsupervised weighting model to identify the most probable sample label for each cell. Ensemblux weighs each constituent tool’s assignment probability distribution by its estimated balanced accuracy for the dataset in a framework that was largely inspired by the work of [Large et al.](https://link.springer.com/article/10.1007/s10618-019-00638-y). To estimate the balanced accuracy of a particular constituent tool (e.g. Demuxalot) for real-word datasets lacking ground-truth labels, Ensemblux leverages the cells with a consensus assignment across the three remaining tools (e.g. Demuxlet, Souporcell, and Vireo-GT) as a proxy for ground-truth. The weighted assignment probabilities across all four constituent tools are then used to inform the most probable sample label for each cell.

 <p align="center">
 <img src="https://github.com/mfiorini9/Ensemblux/assets/97498007/df167d7c-76c1-41da-9f3a-689cf8b54883" width="350" height="100">
 </p>

**Figure 3. Graphical representation of the accuracy-weighted probabilistic ensemble component of the Ensemblux framework.**  
 - - - -

### Step 2: Graph-based doublet detection
The graph-based doublet detection component of the Ensemblux framework was implemented to identify doublets that are incorrectly labeled as singlets by the accuracy-weighted probablistic ensemble component (Step 1). To demonstrate Step 2 of the Ensemblux framework we leveraged a simulated pool comprising 24 pooled samples, 17,384 cells, and a 15% doublet rate.

 <p align="center">
 <img src="https://github.com/mfiorini9/Ensemblux/assets/97498007/da8c8b16-cdec-4dd2-9b9e-f2df4fb1f530" width="350" height="100">
 </p>

**Figure 4. Graphical representation of the graph-based doublet detection component of the Ensemblux framework.**  

The graph-based doublet detection component begins by leveraging select variables returned from each constituent tool:

1. Demuxalot: doublet probability;
2. Demuxlet/Freemuxlet: singlet log likelihood – doublet log likelihood;
3. Demuxlet/Freemuxlet: number of single nucleotide polymorphisms (SNP) per cell;
4. Demuxlet/Freemuxlet: number of reads per cell;
5. Souporcell: doublet log probability;
6. Vireo: doublet probability;
7. Vireo: doublet log likelihood ratio 

 <p align="center">
 <img src="https://github.com/mfiorini9/Ensemblux/assets/97498007/7eb2aa3f-b0d3-45f8-9763-940cc392ef1f" width="650" height="100">
 </p>

**Figure 5. Select variables returned by the constituent genetic demultiplexing tools used for graph-based doubet detection.**

Using these variables, Ensemblux screens each pooled cell to identify the *n* most confident doublets in the pool and performs a principal component analysis (PCA).

 <p align="center">
 <img src="https://github.com/mfiorini9/Ensemblux/assets/97498007/001919fd-bc3b-4b35-bb30-8ee34e084b3e" width="650" height="100">
 </p>

**Figure 6. PCA of pooled cells using select variables returned by the constituent genetic demultiplexing tools.** **A)** PCA highlighting ground truth cell labels: singlet or doublet. **B)** PCA highlighting the *n* most confident doublets identified by Ensemblux.

The PCA embedding is then converted into a Euclidean distance matrix and each cell is assigned a percentile rank based on their distance to each confident doublet. After performing an automated parameter sweep, Ensemblux identifies the droplets that appear most frequently amongst the nearest neighbours of confident doublets as doublets.

 <p align="center">
 <img src="https://github.com/mfiorini9/Ensemblux/assets/97498007/cd21ff69-eee0-4234-a058-1cbe8ead52c4" width="650" height="100">
 </p>

**Figure 7. PCA of pooled cells labeled according to Ensemblux labels prior to and after graph-based doublet detection.** **A)** PCA highlighting ground truth cell labels: singlet or doublet. **B)** PCA highlighting  Ensemblux's labels prior to graph-based doublet detection. **C)** PCA highlighting  Ensemblux's labels after graph-based doublet detection. 

 - - - -

### Step 3: Ensemble-independent doublet detection
The ensemble-independent doublet detection component of the Ensemblux framework was implemented to further improve Ensemblux's ability to identify doublets. Benchmarking on simulated pools with known ground-truth sample labels revealed that certain genetic demultiplexing tools, namely Demuxalot and Vireo, showed high doublet detection specificity.

 <p align="center">
 <img src="https://github.com/mfiorini9/Ensemblux/assets/97498007/e064357a-53bb-41fa-9bd3-f6609c2eefe4" width="350" height="100">
 </p>

**Figure 8. Constituent genetic demultiplexing tool doublet specificity on computationally multiplexed pools with ground truth sample labels.** Doublet specificity was evaluated on pools ranging in size from 4 to 80 multiplexed samples. 

 However, Steps 1 and 2 of the Ensemblux workflow failed to correctly label a subset of doublet calls by these tools. To mitigate this issue and maximize the rate of doublet identification, Ensemblux labels the cells that are identified as doublets by Vireo or Demuxalot as doublets, by default; however, users can nominate different tools for the ensemble-independent doublet detection component depending on the desired doublet detection stringency. 

 <p align="center">
 <img src="https://github.com/mfiorini9/Ensemblux/assets/97498007/9d9bb4fd-f46d-47fe-9aea-a65fc423fbea" width="250" height="100">
 </p>

**Figure 9. Graphical representation of the ensemble-independent doublet detection component of the Ensemblux framework.**  

 - - - -
### Contribution of each step to overall demultiplexing accuracy 
We sequentially applied each step of the Ensemblux framework to 96 computationally multiplexed pools with known ground truth sample labels ranging in size from 4 to 80 samples. The proportion of correctly classified singlets and doublets identified by Ensemblux after each step of the framework is shown in Figure 10. 

 <p align="center">
 <img src="https://github.com/mfiorini9/Ensemblux/assets/97498007/073ed0f8-3535-4e7d-9668-bec0b9450bde" width="750" height="100">
 </p>

**Figure 10. Contribution of each component of the Ensemblux framework to demultiplexing accuracy.** The average proportion of correctly classified **A)** singlets and **B)** doublets across replicates at a given pool size is shown after sequentially applying each step of the Ensemblex framework. The right panels show the average proportion of correct classifications across all 96 pools. The blue points show the proportion of cells that were correctly classified by at least one tool: Demuxalot, Demuxlet, Souporcell, or Vireo.   


 - - - -
For detailed methodology please see our **pre-print manuscript**.







