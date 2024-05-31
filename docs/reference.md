## Adjustable execution parameters for the Ensemblux pipeline

- [Introduction](#introduction)
- [How to modify the parameter files](#how-to-modify-the-parameter-files)
- [Constituent genetic demultiplexing tools with prior genotype information](#step-parameters)
    - [Demuxalot](#step-1-fastq-to-gene-expression-matrix-standard-track)
    - [Demuxlet](#step-1-fastq-to-gene-expression-matrix-hto-track)
    - [Souporcell](#step-2-create-seurat-object-and-remove-ambient-rna)
    - [Vireo](#step-3-quality-control-and-generation-of-filtered-data-objects)
- [Constituent genetic demultiplexing tools without prior genotype information](#step-parameters)
    - [Demuxalot](#step-1-fastq-to-gene-expression-matrix-standard-track)
    - [Freemuxlet](#step-1-fastq-to-gene-expression-matrix-hto-track)
    - [Souporcell](#step-2-create-seurat-object-and-remove-ambient-rna)
    - [Vireo](#step-3-quality-control-and-generation-of-filtered-data-objects)
- [Ensemblux algorithm](#step-parameters) 

- - - -

## Introduction
Prior to running the Ensemblux pipeline, users should modify the execution parameters for the constituent genetic demultiplexing tools and the Ensemblux algorithm. Upon running [Step 1: Set up](Step0.md), a `/job_info` folder will be created in the wording directory. Within the `/job_info` folder is a `/configs` folder which contains the `ensemblux_config.ini`; this .ini file contains all of the adjustable parameters for the Ensemblux pipeline. 

```
working_directory
└── job_info
    ├── configs
    │   └── ensemblux_config.ini
    ├── logs
    └── summary_report.txt
```
To ensure replicability, the execution parameters are documented in `~/working_directory/job_info/summary_report.txt`.

- - - -
## How to modify the parameter files
The following section illustrates how to modify the `ensemblux_config.ini` parameter file directly from the terminal. To begin, navigate to the `/configs` folder and view its contents:

```
cd ~/working_directory/job_info/configs
ls
```
The following file will be available: `ensemblux_config.ini`

To modify the `ensemblux_config.ini` parameter file directly in the terminal we will use [Nano](https://help.ubuntu.com/community/Nano):

```
nano ensemblux_config.ini
```
This will open `ensemblux_config.ini` in the terminal and allow users to modify the parameters. To save the modifications and exit the parameter file, type `ctrl+o` followed by `ctrl+x`.
- - - -

## Constituent genetic demultiplexing tools with prior genotype information

#### Demuxalot

The following parameters are adjustable for Demuxalot:

|Parameter|Default|Description|
|:--|:--|:--|
|PAR_demuxalot_genotype_names|NULL| List of Sample ID's in the sample VCF file (e.g., 'Sample_1,Sample_2,Sample_3').|
|PAR_demuxalot_minimum_coverage|200| Minimum read coverage.|
|PAR_demuxalot_minimum_alternative_coverage|10| Minimum alternative read coverage.|
|PAR_demuxalot_n_best_snps_per_donor|100| Number of best snps for each donor to use for demultiplexing.|
|PAR_demuxalot_genotypes_prior_strength|1| Genotype prior strength.|
|PAR_demuxalot_doublet_prior|0.25| Doublet prior strength. |

- - - -

#### Demuxlet

The following parameters are adjustable for Demuxlet:

|Parameter|Default|Description|
|:--|:--|:--|
|PAR_demuxlet_field|GT| Field to extract the genotypes (GT), genotype likelihood (PL), or posterior probability (GP) from the sample .vcf file.|

**NOTE**: We are currently working on expanding the execution parameters for Demuxlet. 

- - - -

#### Vireo

The following parameters are adjustable for Vireo:

|Parameter|Default|Description|
|:--|:--|:--|
|PAR_vireo_N|NULL| Number of pooled samples.|
|PAR_vireo_type|GT| Field to extract the genotypes (GT), genotype likelihood (PL), or posterior probability (GP) from the sample .vcf file.|
|PAR_vireo_processes|20| Number of subprocesses for computing.|
|PAR_vireo_minMAF|0.1| Minimum minor allele frequency. |
|PAR_vireo_minCOUNT|20| Minimum aggregated count. |
|PAR_vireo_forcelearnGT|T| Whether or not to treat donor GT as prior only.|

**NOTE**: We are currently working on expanding the execution parameters for Vireo. 

- - - -

#### Souporcell

The following parameters are adjustable for Souporcell:

|Parameter|Default|Description|
|:--|:--|:--|
|PAR_minimap2|-ax splice -t 8 -G50k -k 21 -w 11 --sr -A2 -B8 -O12,32 -E2,1 -r200 -p.5 -N20 -f1000,5000 -n2 -m20 -s40 -g2000 -2K50m --secondary=no| For information regarding the minimap2 parameters, please see the [documentation](https://github.com/lh3/minimap2).  |
|PAR_freebayes|-iXu -C 2 -q 20 -n 3 -E 1 -m 30 --min-coverage 6| For information regarding the freebayes parameters, please see the [documentation](https://github.com/freebayes/freebayes). |
|PAR_vartrix_umi|TRUE| Whether or no to consider UMI information when populating coverage matrices.|
|PAR_vartrix_mapq|30| Minimum read mapping quality. |
|PAR_vartrix_threads|8| Number of threads for computing.|
|PAR_souporcell_k|NULL| Number of pooled samples. |
|PAR_souporcell_t|8| Number of threads for computing. |


**NOTE**: We are currently working on expanding the execution parameters for Souporcell. 

- - - -

## Constituent genetic demultiplexing tools without prior genotype information

#### Demuxalot

The following parameters are adjustable for Demuxalot:

|Parameter|Default|Description|
|:--|:--|:--|
|PAR_demuxalot_genotype_names|NULL| List of Sample ID's in the sample VCF file generated by Freemuxlet: outs.clust1.vcf (e.g., 'CLUST0,CLUST1,CLUST2').|
|PAR_demuxalot_minimum_coverage|200| Minimum read coverage.|
|PAR_demuxalot_minimum_alternative_coverage|10| Minimum alternative read coverage.|
|PAR_demuxalot_n_best_snps_per_donor|100| Number of best snps for each donor to use for demultiplexing.|
|PAR_demuxalot_genotypes_prior_strength|1| Genotype prior strength.|
|PAR_demuxalot_doublet_prior|0.25| Doublet prior strength. |

- - - -

#### Freemuxlet

The following parameters are adjustable for Freemuxlet:

|Parameter|Default|Description|
|:--|:--|:--|
|PAR_freemuxlet_nsample|NULL| Number of pooled samples.|

**NOTE**: We are currently working on expanding the execution parameters for Freemuxlet. 

- - - -

#### Vireo

The following parameters are adjustable for Vireo:

|Parameter|Default|Description|
|:--|:--|:--|
|PAR_vireo_N|NULL| Number of pooled samples.|
|PAR_vireo_processes|20| Number of subprocesses for computing.|
|PAR_vireo_minMAF|0.1| Minimum minor allele frequency. |
|PAR_vireo_minCOUNT|20| Minimum aggregated count. |

**NOTE**: We are currently working on expanding the execution parameters for Vireo. 

- - - -

#### Souporcell

The following parameters are adjustable for Souporcell:

|Parameter|Default|Description|
|:--|:--|:--|
|PAR_minimap2|-ax splice -t 8 -G50k -k 21 -w 11 --sr -A2 -B8 -O12,32 -E2,1 -r200 -p.5 -N20 -f1000,5000 -n2 -m20 -s40 -g2000 -2K50m --secondary=no| For information regarding the minimap2 parameters, please see the [documentation](https://github.com/lh3/minimap2).  |
|PAR_freebayes|-iXu -C 2 -q 20 -n 3 -E 1 -m 30 --min-coverage 6| For information regarding the freebayes parameters, please see the [documentation](https://github.com/freebayes/freebayes). |
|PAR_vartrix_umi|TRUE| Whether or no to consider UMI information when populating coverage matrices.|
|PAR_vartrix_mapq|30| Minimum read mapping quality. |
|PAR_vartrix_threads|8| Number of threads for computing.|
|PAR_souporcell_k|NULL| Number of pooled samples. |
|PAR_souporcell_t|8| Number of threads for computing. |


**NOTE**: We are currently working on expanding the execution parameters for Souporcell. 

- - - -

## Ensemblux

The following parameters are adjustable for the Ensemblux algorithm:

|Parameter|Default|Description|
|:--|:--|:--|
|<ins>**Pool parameters**</ins>|
|PAR_ensemblux_sample_size| NULL| Number of samples multiplexed in the pool.|
|PAR_ensemblux_expected_doublet_rate| NULL|Expected doublet rate for the pool. If using 10X Genomics, the expected doublet rate can be estimated based on the number of recovered cells. For more information see [10X Genomics Documentation](https://kb.10xgenomics.com/hc/en-us/articles/360059124751-Why-is-the-multiplet-rate-different-for-the-Next-GEM-Single-Cell-3-LT-v3-1-assay-compared-to-other-single-cell-applications).|
|<ins>**Set up parameters**</ins>|
|PAR_ensemblux_merge_constituents|Yes|Whether or not to merge the output files of the constituent demultiplexing tools. If running Ensemblux on a pool for the first time, this parameter should be set to "Yes". Subsequent runs of Ensemblux (e.g., parameter optimization) can have this parameter set to "No" as the pipeline will automatically detect the previously generated merged file.|
|<ins>**Step 1 parameters: Probabilistic-weighted ensemble**</ins>|
|PAR_ensemblux_probabilistic_weighted_ensemble| Yes|Whether or not to perform Step 1: Probabilistic-weighted ensemble. If running Ensemblux on a pool for the first time, this parameter should be set to "Yes". Subsequent runs of Ensemblux (e.g., parameter optimization) can have this parameter set to "No" as the pipeline will automatically detect the previously generated Step 1 output file.|
|<ins>**Step 2 parameters: Graph-based doublet detection**</ins>|
|PAR_ensemblux_preliminary_parameter_sweep| No|Whether or not to perform a preliminary parameter sweep for Step 2: Graph-based doublet detection. Users should utilize the preliminary parameter sweep if they wish to manually define the number of confident doublets in the pool (nCD) and the percentile threshold of the nearest neighour frequency (pT), which can be defined in the following two parameters, respectively. |
|PAR_ensemblux_nCD| NULL|Manually defined number of confident doublets in the pool (nCD). Value can be informed by the output files generated by setting PAR_ensemblux_preliminary_parameter_sweep to "Yes". |
|PAR_ensemblux_pT| NULL|Manually defined percentile threshold of the nearest neighour frequency (pT). Value can be informed by the output files generated by setting PAR_ensemblux_preliminary_parameter_sweep to "Yes".|
|PAR_ensemblux_graph_based_doublet_detection| Yes|Whether or not to perform Step 2: Graph-based doublet detection. If PAR_ensemblux_nCD and PAR_ensemblux_pT are not defined by the user (NULL), Ensemblux will automatically determine the optimal parameter values using an unsupervised parameter sweep. If PAR_ensemblux_nCD and PAR_ensemblux_pT are defined by the user, graph-based doublet detection will be performed with the user-defined values.   |
|<ins>**Step 3 parameters: Ensemble-independent doublet detection**</ins>|
|PAR_ensemblux_preliminary_ensemble_independent_doublet| No|Whether or not to perform a preliminary parameter sweep for Step 3: Ensemble-independent doublet detection. Users should utilize the preliminary parameter sweep if they wish to manually define which constituent tools to utilize for ensemble-independent doublet detection. Users can define which tools to utilize for ensemble-independent doublet detection in the following parameters.|
|PAR_ensemblux_ensemble_independent_doublet| Yes|Whether or not to perform Step 3: Ensemble-independent doublet detection.|
|PAR_ensemblux_doublet_Demuxalot_threshold| Yes|Whether or not to label doublets identified by Demuxalot as doublets. Only doublets with assignment probabilities exceeding Demuxalot's recommended probability threshold will be labeled as doublets by Ensemblux.|
|PAR_ensemblux_doublet_Demuxalot_no_threshold| No|Whether or not to label doublets identified by Demuxalot as doublets, regardless of the corresponding assignment probability.|
|PAR_ensemblux_doublet_Demuxlet_threshold| No|Whether or not to label doublets identified by Demuxlet as doublets. Only doublets with assignment probabilities exceeding Demuxlet's recommended probability threshold will be labeled as doublets by Ensemblux.|
|PAR_ensemblux_doublet_Demuxlet_no_threshold| No|Whether or not to label doublets identified by Demuxlet as doublets, regardless of the corresponding assignment probability.|
|PAR_ensemblux_doublet_Souporcell_threshold| No|Whether or not to label doublets identified by Souporcell as doublets. Only doublets with assignment probabilities exceeding Souporcell's recommended probability threshold will be labeled as doublets by Ensemblux.|
|PAR_ensemblux_doublet_Souporcell_no_threshold| No|Whether or not to label doublets identified by Souporcell as doublets, regardless of the corresponding assignment probability.|
|PAR_ensemblux_doublet_Vireo_threshold| Yes|Whether or not to label doublets identified by Vireo as doublets. Only doublets with assignment probabilities exceeding Vireo's recommended probability threshold will be labeled as doublets by Ensemblux.|
|PAR_ensemblux_doublet_Vireo_no_threshold| No|Whether or not to label doublets identified by Vireo as doublets, regardless of the corresponding assignment probability.|
|<ins>**Confidence score parameters**</ins>|
|PAR_ensemblux_compute_singlet_confidence| Yes|Whether or not to compute Ensemblux's singlet confidence score. This will define low confidence assignments which should be removed from downstream analyses. |




