# Installation
The Ensemblex container is freely available under an MIT open-source license at [https://zenodo.org/records/11639103](https://zenodo.org/records/11639103). 

The Ensemblex container can be downloaded using the following code:

```
## Download the Ensemblex container
curl "https://zenodo.org/records/11639103/files/ensemblex.pip.zip?download=1" --output ensemblex.pip.zip

## Unzip the Ensemblex container
unzip ensemblex.pip.zip
```

If installation was successful the following will be available:
```
ensemblex.pip
├── gt
│   ├── configs
│   │   └── ensemblex_config.ini
│   └── scripts
│       ├── demuxalot
│       │   ├── pipeline_demuxalot.sh
│       │   └── pipline_demuxalot.py
│       ├── demuxlet
│       │   └── pipeline_demuxlet.sh
│       ├── ensemblexing
│       │   ├── ensemblexing.R
│       │   ├── functions.R
│       │   └── pipeline_ensemblexing.sh
│       ├── souporcell
│       │   └── pipeline_souporcell_generate.sh
│       └── vireo
│           └── pipeline_vireo.sh
├── launch
│   ├── launch_gt.sh
│   └── launch_nogt.sh
├── launch_ensemblex.sh
├── nogt
│   ├── configs
│   │   └── ensemblex_config.ini
│   └── scripts
│       ├── demuxalot
│       │   ├── pipeline_demuxalot.py
│       │   └── pipeline_demuxalot.sh
│       ├── ensemblexing
│       │   ├── ensemblexing_nogt.R
│       │   ├── functions_nogt.R
│       │   └── pipeline_ensemblexing.sh
│       ├── freemuxlet
│       │   └── pipeline_freemuxlet.sh
│       ├── souporcell
│       │   └── pipeline_souporcell_generate.sh
│       └── vireo
│           └── pipeline_vireo.sh
├── README
├── soft
│   └── ensemblex.sif
└── tools
    ├── sort_vcf_same_as_bam.sh
    └── utils.sh
```

In addition to the Ensemblex container, users must install [Apptainer](https://apptainer.org/). For example:
```
## Load Apptainer
module load apptainer/1.2.4 
```

To test if the Ensemblex container is installed properly, run the following code:
```
## Define the path to ensemblex.pip
ensemblex_HOME=/path/to/ensemblex.pip

## Print help message
bash $ensemblex_HOME/launch_ensemblex.sh -h
```
Which should return the following help message:
```
------------------- 
Usage:  /home/fiorini9/scratch/ensemblex.pip/launch_ensemblex.sh [arguments]
        mandatory arguments:
                -d  (--dir)  = Working directory (where all the outputs will be printed) (give full path)
                --steps  =  Specify the steps to execute. Begin by selecting either init-GT or init-noGT to establish the working directory. 
                       For GT: vireo, demuxalot, demuxlet, souporcell, ensemblexing 
                       For noGT: vireo, demuxalot, freemuxlet, souporcell, ensemblexing 

        optional arguments:
                -h  (--help)  = See helps regarding the pipeline arguments 
                --vcf  = The path of vcf file 
                --bam  = The path of bam file 
                --sortout  = The path snd nsme of vcf generated using sort  
 ------------------- 
 For a comprehensive help, visit  https://neurobioinfo.github.io/ensemblex/site/ for documentation. 
```
----

Upon installing up the Ensemblex container, we can proceed to Step 1 where we will initiate the Ensemblex pipeline for demultiplexing: [Set up](Step0.md)