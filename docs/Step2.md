# Step 3: Genetic demultiplexing by constituent demultiplexing tools
In Step 3, we will demultiplex the pooled samples with each of Ensemblux's constituent genetic demultiplexing tools. The constituent genetic demultiplexing tools will vary depending on the version of the Ensemblux pipeline being used:

- [Demultiplexing with prior genotype information](#demultiplexing-with-prior-genotype-information)
- [Demultiplexing without prior genotype information](#demultiplexing-without-prior-genotype-information)

**NOTE**: The analytical parameters for each constiuent tool can be adjusted using the the `ensemblux_config.ini` file located in `~/working_directory/job_info/configs`. For a comprehensive description of how to adjust the analytical parameters of the Ensemblux pipeline please see [Execution parameters](reference.md). 
 - - - -

## Demultiplexing with prior genotype information
When demultiplexing **with** prior genotype information, Ensemblux leverages the sample labels from

- [Demuxalot](#demuxalot)
- [Demuxlet](#demuxlet)
- [Souporcell](#souporcell)
- [Vireo-GT](#vireo-gt)
 - - - -
#### Demuxalot
To run Demuxalot use the following code:

```
Ensemblux_HOME=/path/to/ensemblux.pip
Ensemblux_PWD=/path/to/working_directory

bash $Ensemblux_HOME/launch_ensemblux.sh -d $Ensemblux_PWD --step demuxalot
```
If Demuxalot completed successfully, the following files should be available in `~/working_directory/demuxalot`

```
working_directory
└── demuxalot
    ├── Demuxalot_result.csv
    └── new_snps_single_file.betas
```
 - - - -
#### Demuxlet
To run Demuxlet use the following code:

```
Ensemblux_HOME=/path/to/ensemblux.pip
Ensemblux_PWD=/path/to/working_directory

bash $Ensemblux_HOME/launch_ensemblux.sh -d $Ensemblux_PWD --step demuxlet
```
If Demuxlet completed successfully, the following files should be available in `~/working_directory/demuxlet`

```
working_directory
└── demuxlet
    ├── outs.best
    ├── pileup.cel.gz
    ├── pileup.plp.gz
    ├── pileup.umi.gz
    └── pileup.var.gz
```
 - - - -
#### Souporcell
To run Souporcell use the following code:

```
Ensemblux_HOME=/path/to/ensemblux.pip
Ensemblux_PWD=/path/to/working_directory

bash $Ensemblux_HOME/launch_ensemblux.sh -d $Ensemblux_PWD --step souporcell
```

If Souporcell completed successfully, the following files should be available in `~/working_directory/souporcell`

```
working_directory
└── souporcell
    ├── alt.mtx
    ├── cluster_genotypes.vcf
    ├── clusters_tmp.tsv
    ├── clusters.tsv
    ├── fq.fq
    ├── minimap.sam
    ├── minitagged.bam
    ├── minitagged_sorted.bam
    ├── minitagged_sorted.bam.bai
    ├── Pool.vcf
    ├── ref.mtx
    └── soup.txt
```

 - - - -
#### Vireo-GT
To run Vireo-GT use the following code:

```
Ensemblux_HOME=/path/to/ensemblux.pip
Ensemblux_PWD=/path/to/working_directory

bash $Ensemblux_HOME/launch_ensemblux.sh -d $Ensemblux_PWD --step vireo
```

If Vireo-GT completed successfully, the following files should be available in `~/working_directory/vireo_gt`

```
working_directory
└── vireo_gt
    ├── cellSNP.base.vcf.gz
    ├── cellSNP.cells.vcf.gz
    ├── cellSNP.samples.tsv
    ├── cellSNP.tag.AD.mtx
    ├── cellSNP.tag.DP.mtx
    ├── cellSNP.tag.OTH.mtx
    ├── donor_ids.tsv
    ├── fig_GT_distance_estimated.pdf
    ├── fig_GT_distance_input.pdf
    ├── GT_donors.vireo.vcf.gz
    ├── _log.txt
    ├── prob_doublet.tsv.gz
    ├── prob_singlet.tsv.gz
    └── summary.tsv
```
 - - - -
  Upon demultiplexing the pooled samples with each of Ensemblux's constituent genetic demultiplexing tools, we can proceed to Step 4 where we will process the output files of the consituent tools with the Ensemblux algorithm to generate the ensemble sample classifications: [Application of Ensemblux](Step3.md)
  - - - - 


## Demultiplexing without prior genotype information
When demultiplexing **without** prior genotype information, Ensemblux leverages the sample labels from

- [Freemuxlet](#freemuxlet)
- [Souporcell](#souporcell)
- [Vireo](#vireo)
- [Demuxalot](#demuxalot)

 - - - -
#### Freemuxlet
To run Freemuxlet use the following code:

```
Ensemblux_HOME=/path/to/ensemblux.pip
Ensemblux_PWD=/path/to/working_directory

bash $Ensemblux_HOME/launch_ensemblux.sh -d $Ensemblux_PWD --step freemuxlet
```
If Freemuxlet completed successfully, the following files should be available in `~/working_directory/freemuxlet`

```
working_directory
└── freemuxlet
    ├── outs.clust1.samples.gz
    ├── outs.clust1.vcf
    ├── outs.lmix
    ├── pileup.cel.gz
    ├── pileup.plp.gz
    ├── pileup.umi.gz
    └── pileup.var.gz
```
 - - - -
#### Souporcell
To run Souporcell use the following code:

```
Ensemblux_HOME=/path/to/ensemblux.pip
Ensemblux_PWD=/path/to/working_directory

bash $Ensemblux_HOME/launch_ensemblux.sh -d $Ensemblux_PWD --step souporcell
```

If Souporcell completed successfully, the following files should be available in `~/working_directory/souporcell`

```
working_directory
└── souporcell
    ├── alt.mtx
    ├── cluster_genotypes.vcf
    ├── clusters_tmp.tsv
    ├── clusters.tsv
    ├── fq.fq
    ├── minimap.sam
    ├── minitagged.bam
    ├── minitagged_sorted.bam
    ├── minitagged_sorted.bam.bai
    ├── Pool.vcf
    ├── ref.mtx
    └── soup.txt
```

 - - - -
#### Vireo
To run Vireo use the following code:

```
Ensemblux_HOME=/path/to/ensemblux.pip
Ensemblux_PWD=/path/to/working_directory

bash $Ensemblux_HOME/launch_ensemblux.sh -d $Ensemblux_PWD --step vireo
```

If Vireo completed successfully, the following files should be available in `~/working_directory/vireo`

```
working_directory
└── vireo
    ├── cellSNP.base.vcf.gz
    ├── cellSNP.cells.vcf.gz
    ├── cellSNP.samples.tsv
    ├── cellSNP.tag.AD.mtx
    ├── cellSNP.tag.DP.mtx
    ├── cellSNP.tag.OTH.mtx
    ├── donor_ids.tsv
    ├── fig_GT_distance_estimated.pdf
    ├── GT_donors.vireo.vcf.gz
    ├── _log.txt
    ├── prob_doublet.tsv.gz
    ├── prob_singlet.tsv.gz
    └── summary.tsv
```
 - - - -

#### Demuxalot
**NOTE**: Because the Demuxalot algorithm requires prior genotype information, the Ensemblux pipeline uses the predicted vcf file generated by Freemuxlet as input into Demuxalot when prior genotype information is not available. Therefore, it is important to wait for Freemuxlet to complete before running Demuxalot. To check if the required Freemuxlet-generated vcf file is available prior to running Demuxalot, you can use the following code:

```
if test -f /path/to/working_directory/freemuxlet/outs.clust1.vcf; then
  echo "File exists."
fi
```

Upon confirming that the required Freemuxlet-generated file exists, we can run Demuxalot using the following code:

```
Ensemblux_HOME=/path/to/ensemblux.pip
Ensemblux_PWD=/path/to/working_directory

bash $Ensemblux_HOME/launch_ensemblux.sh -d $Ensemblux_PWD --step demuxalot
```
If Demuxalot completed successfully, the following files should be available in `~/working_directory/demuxalot`

```
working_directory
└── demuxalot
    ├── Demuxalot_result.csv
    └── new_snps_single_file.betas
```
 - - - -

 - - - -
  Upon demultiplexing the pooled samples with each of Ensemblux's constituent genetic demultiplexing tools, we can proceed to Step 4 where we will process the output files of the consituent tools with the Ensemblux algorithm to generate the ensemble sample classifications: [Application of Ensemblux](Step3.md)
