## Ensemblux algorithm outputs
- [Introduction](#introduction)
- [Outputs](#outputs)
    - [Merging constituent output files](#merging-constituent-output-files)
    - [Step 1: Accuracy-weighted probabilistic ensemble](#step-1-accuracy-weighted-probabilistic-ensemble)
    - [Step 2: Graph-based doublet detection](#step-2-graph-based-doublet-detection)
    - [Step 3: Ensemble-independent doublet detection](#step-3-ensemble-independent-doublet-detection)
    - [Singlet confidence score](#singlet-confidence-score)
- - - -

## Introduction
After applying the Ensemblux algorithm to the output files of the constituent genetic demultiplexing tools in Step 4, the `~/working_directory/ensemblux` folder will have the following structure:

```
working_directory
└── ensemblux
    ├── constituent_tool_merge.csv
    ├── PWE
    ├── GBD
    ├── EID
    └── confidence
```
 - `constituent_tool_merge.csv` is the merged outputs from each constituent genetic demultiplexing tool.
 - `PWE/` contains the outputs from Step 1: probabilistic-weighted ensemble (PWE).
 - `GBD/` contains the outputs from Step 2: graph-based doublet detection (GBD).
 - `EID/` contains the outputs from Step 3: ensemble-independent doublet detection (EID).
 - `confidence/` contains the final Ensemblux output file, whose sample labels have been annotate with the Ensemblux signlet confidence score.


 **Note:** If users re-run a step of the Ensemblux workflow, the outputs from the previous run will automatically be overwritten. If you do not want to lose the outputs from a previous run, it is important to copy the materials to a separate directory. 
- - - -
## Outputs
#### Merging constituent output files
Ensemblux begins by merging the output files of the constituent genetic demultiplexing tools by cell barcode, which produces the `constituent_tool_merge.csv` file. In this file, each constituent genetic demultiplexing tool has two columns corresponding to their sample labels:

 - `demuxalot_assignment`
 - `demuxalot_best_assignment`
 - `demuxlet_assignment`
 - `demuxlet_best_assignment`
 - `souporcell_assignment`
 - `souporcell_best_assignment`
 - `vireo_assignment`
 - `vireo_best_assignment`

Taking Vireo as an example, `vireo_assignment` shows Vireo's sample labels after applying its recommended probability threshold; thus, cells that do not meet Vireo's recommended probability threshold will be labeled as "unassigned". In turn, `vireo_best_assignment` shows Vireo's best guess assignments with out applying the recommended probability threshold; thus, cells that do not meet Vireo's recommended probability threshold will still show the best sample label and will not be labelled as "unassigned". 

The `constituent_tool_merge.csv` file also contains a `general_consensus` column. **This is not Ensemblux's sample labels**. The `general_consensus` column simply shows the sample labels that result from a majority vote classifier; split decisions are labeled as unassigned. 

- - - -

#### Step 1: Accuracy-weighted probabilistic ensemble
After running Step 1 of the Ensemblux algorithm, the `/PWE` folder will contain the following files:

```
working_directory
└── ensemblux
    └── PWE
        ├── ARI_demultiplexing_tools.pdf
        ├── BA_demultiplexing_tools.pdf
        ├── Balanced_accuracy_summary.csv
        └── Step1_cell_assignment.csv
```

|Output type |Name|Description|
|:--|:--|:--|
|Figure|ARI_demultiplexing_tools.pdf | Heatmap showing the Adjusted Rand Index (ARI) between the sample labels of the constituent genetic demultiplexing tools. |
|Figure|BA_demultiplexing_tools.pdf | Barplot showing the estimated balanced accuracy for each constituent genetic demultiplexing tool. |
|File |Balanced_accuracy_summary.csv | Summary file describing the estimated balanced accuracy computation for each constituent genetic demultiplexing tool. |
|File| Step1_cell_assignment.csv|Data file containing Ensemblux's sample labels after Step 1: accuracy-weighted probabilistic ensemble.|

The `Step1_cell_assignment.csv` file contains the following important columns:

 - `Ensemblux_assignment`: Ensemblux sample labels after performing accuracy-weighted probabilistic ensemble.
 - `Ensemblux_probability`: Accuracy-weighted ensemble probability corresponding to Ensemblux's sample labels.

**NOTE**: Prior to using Ensemblux's sample labels for downstream analyses, we recommend computing the Ensemblux singlet confidence score to identify low confidence singlet assignments that should be removed from the dataset to mitigate the introduction of technical artificats. 


- - - -

#### Step 2: Graph-based doublet detection
After running Step 2 of the Ensemblux algorithm, the `/GBD` folder will contain the following files:

```
working_directory
└── ensemblux
    └── PWE
        ├── optimal_nCD.pdf
        ├── optimal_pT.pdf
        ├── PC1_var_contrib.pdf
        ├── PC2_var_contrib.pdf
        ├── PCA1_graph_based_doublet_detection.pdf
        ├── PCA2_graph_based_doublet_detection.pdf
        ├── PCA3_graph_based_doublet_detection.pdf
        ├── PCA_plot.pdf
        ├── PCA_scree_plot.pdf
        └── Step2_cell_assignment.csv
```

|Output type |Name|Description|
|:--|:--|:--|
|Figure|optimal_nCD.pdf | Dot plot showing the optimal nCD value.  |
|Figure|optimal_pT.pdf | Dot plot showing the optimal pT value. |
|Figure|PC1_var_contrib.pdf | Bar plot showing the contribution of each variable to the variation across the first principal component.  |
|Figure|PC2_var_contrib.pdf |Bar plot showing the contribution of each variable to the variation across the second principal component.  |
|Figure|PCA1_graph_based_doublet_detection.pdf | PCA showing Ensemblux sample labels (singlet or doublet) prior to performing graph-based doublet detection.   |
|Figure|PCA2_graph_based_doublet_detection.pdf | PCA showing the cells identified as the *n* most confident doublets in the pool. |
|Figure|PCA3_graph_based_doublet_detection.pdf | PCA showing Ensemblux sample labels (singlet or doublet) after performing graph-based doublet detection.   |
|Figure|PCA_plot.pdf | PCA of pooled cells. |
|Figure|PCA_scree_plot.pdf | Bar plot showing the variance explained by each principal component. |
|File|Step2_cell_assignment.csv | Data file containing Ensemblux's sample labels after Step 2: graph-based doublet detection. |



The `Step2_cell_assignment.csv` file contains the following important column:

 - `Ensemblux_assignment`: Ensemblux sample labels after performing graph-based doublet detection.

**NOTE**: Prior to using Ensemblux's sample labels for downstream analyses, we recommend computing the Ensemblux singlet confidence score to identify low confidence singlet assignments that should be removed from the dataset to mitigate the introduction of technical artificats. 

- - - -

#### Step 3: Ensemble-independent doublet detection
After running Step 3 of the Ensemblux algorithm, the `/EID` folder will contain the following files:

```
working_directory
└── ensemblux
    └── EID
        ├── Doublet_overlap_no_threshold.pdf
        ├── Doublet_overlap_threshold.pdf
        ├── Number_Ensemblux_doublets_EID_no_threshold.pdf
        ├── Number_Ensemblux_doublets_EID_threshold.pdf
        └── Step3_cell_assignment.csv

```

|Output type |Name|Description|
|:--|:--|:--|
|Figure|Doublet_overlap_no_threshold.pdf | Proportion of doublet calls overlapping between constituent genetic demultiplexing tools without applying assignment probability thresholds. |
|Figure|Doublet_overlap_threshold.pdf | Proportion of doublet calls overlapping between constituent genetic demultiplexing tools after applying assignment probability thresholds. |
|Figure |Number_Ensemblux_doublets_EID_no_threshold.pdf | Number of cells that would be labelled as doublets by Ensemblux if a constituent tool was nominated for ensemble-independent doublet detection, without applying assignment probability thresholds.  |
|Figure| Number_Ensemblux_doublets_EID_threshold.pdf|Number of cells that would be labelled as doublets by Ensemblux if a constituent tool was nominated for ensemble-independent doublet detection, after applying assignment probability thresholds. |
|File| Step3_cell_assignment.csv|Data file containing Ensemblux's sample labels after Step 3: ensemble-independent doublet detection.|

The `Step3_cell_assignment.csv` file contains the following important column:

 - `Ensemblux_assignment`: Ensemblux sample labels after performing ensemble-independent doublet detection.

**NOTE**: Prior to using Ensemblux's sample labels for downstream analyses, we recommend computing the Ensemblux singlet confidence score to identify low confidence singlet assignments that should be removed from the dataset to mitigate the introduction of technical artificats. 

- - - -

#### Singlet confidence score

After computing the Ensemblux singlet confidence score, the `/confidence` folder will contain the following file:

```
working_directory
└── ensemblux
    └── confidence
        └── Ensemblux_final_cell_assignment.csv


```

|Output type |Name|Description|
|:--|:--|:--|
|File| Ensemblux_final_cell_assignment.csv|Data file containing Ensemblux's final sample labels after computing the singlet confidence score. |

The `Ensemblux_final_cell_assignment.csv` file contains the following important column:

 - `Ensemblux_assignment`: Ensemblux sample labels after applying the recommended singlet confidence score threshold; singlets with a confidence score < 1  are labeled as "unassigned".
 - `Ensemblux_best_assignment`: Ensemblux's best guess assignments with out applying the recommended confidence score threshold; singlets with a confidence score < 1 will still show the best sample label and will not be labelled as "unassigned".
 - `Ensemblux_singlet_confidence`: Ensemblux singlet confidence score.


**NOTE**: We recommend using the sample labels from `Ensemblux_assignment` for downstream analyses. 

