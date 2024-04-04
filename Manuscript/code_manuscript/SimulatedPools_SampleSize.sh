#!/bin/bash

#########
# NOTES #
#########
# This code was used to produce simulated pools with known ground truth sample labels using 80 independently-sequenced iPSC lines differentiated towards a dopaminergic neuronal state as part of the Foundational Data Initiative for Parkinson's Disease
# Simulated pools were produced using the synth_pool.py script produced by the developers of Vireo, ranging in size from 4 to 80 multiplexed sample
# Below is an example script for producing a simulated pool from four lines, totaling ~20,000 cells and a 15% doublet rate

#################
# SET VARIABLES #
#################
POOL=4pool

S1=SCRN_PPMI3186_2662_da65_v1
S2=SCRN_PPMI3220_6139_da65_v1
S3=SCRN_PPMI3234_1338_da65_v1
S4=SCRN_PPMI3409_0660_da65_v1

SYNTHDIR=/home/fiorini9/scratch/mjff/MF/R1_20/$POOL 
DATA=/lustre03/project/6070393/COMMON/rawdata/FoundIN-PD/download/FOUNDIN/processed/SCRN/SCRN/cellRanger
BAM=outs/possorted_genome_bam.bam
SOFTS=/lustre03/project/6070393/COMMON/samplepooling/MF/synthetic_data/synth_submit
REF=/lustre03/project/6070393/COMMON/samplepooling/MF/synthetic_data/VCF/reheader_rename_PDHCSWEDD_progen_merged.vcf 

BARCODES1=/lustre03/project/6070393/COMMON/rawdata/FoundIN-PD/download/FOUNDIN/processed/SCRN/SCRN/cellRanger
BARCODES2=outs/filtered_feature_bc_matrix/barcodes.tsv.gz
DATADIR=/home/fiorini9/scratch/mjff/MF/R1_20 
BAM1=/lustre03/project/6070393/COMMON/rawdata/FoundIN-PD/download/FOUNDIN/processed/SCRN/SCRN/cellRanger
BAM2=outs/possorted_genome_bam.bam

########
# MAIN #
########
#1. copy over barcodes files
cd $SYNTHDIR
cp $BARCODES1/$S1/$BARCODES2 $DATADIR/$POOL/$S1.tsv.gz
cp $BARCODES1/$S2/$BARCODES2 $DATADIR/$POOL/$S2.tsv.gz
cp $BARCODES1/$S3/$BARCODES2 $DATADIR/$POOL/$S3.tsv.gz
cp $BARCODES1/$S4/$BARCODES2 $DATADIR/$POOL/$S4.tsv.gz

gunzip *.gz

#2. copy over the bam files
cp $BAM1/$S1/$BAM2  $SYNTHDIR/$S1.bam
cp $BAM1/$S2/$BAM2  $SYNTHDIR/$S2.bam
cp $BAM1/$S3/$BAM2  $SYNTHDIR/$S3.bam
cp $BAM1/$S4/$BAM2  $SYNTHDIR/$S4.bam

#3. filter barcodes files to only include th number of cells that we want
cat $DATADIR/$POOL/$S1.tsv | awk 'NR >= 0  && NR <= 5000' > BC1.tsv
cat $DATADIR/$POOL/$S2.tsv | awk 'NR >= 0  && NR <= 5000' > BC2.tsv
cat $DATADIR/$POOL/$S3.tsv | awk 'NR >= 0  && NR <= 5000' > BC3.tsv
cat $DATADIR/$POOL/$S4.tsv | awk 'NR >= 0  && NR <= 5000' > BC4.tsv

#4. add BC tag to our new barcodes file
sed -e 's/^/CB:Z:/' BC1.tsv > tag_BC1.tsv
sed -e 's/^/CB:Z:/' BC2.tsv > tag_BC2.tsv
sed -e 's/^/CB:Z:/' BC3.tsv > tag_BC3.tsv
sed -e 's/^/CB:Z:/' BC4.tsv > tag_BC4.tsv

#5. Save the header lines
samtools view -H $SYNTHDIR/$S1.bam > $S1.header
samtools view -H $SYNTHDIR/$S2.bam > $S2.header
samtools view -H $SYNTHDIR/$S3.bam > $S3.header
samtools view -H $SYNTHDIR/$S4.bam > $S4.header

#6. Filter .bam to only include desired barcodes
samtools view $SYNTHDIR/$S1.bam | LC_ALL=C grep -F -f tag_BC1.tsv > $SYNTHDIR/$S1.sam.body
samtools view $SYNTHDIR/$S2.bam | LC_ALL=C grep -F -f tag_BC2.tsv > $SYNTHDIR/$S2.sam.body
samtools view $SYNTHDIR/$S3.bam | LC_ALL=C grep -F -f tag_BC3.tsv > $SYNTHDIR/$S3.sam.body
samtools view $SYNTHDIR/$S4.bam | LC_ALL=C grep -F -f tag_BC4.tsv > $SYNTHDIR/$S4.sam.body

#7. Combine .sam header and body
cat $S1.header $SYNTHDIR/$S1.sam.body > $SYNTHDIR/$S1.sam
cat $S2.header $SYNTHDIR/$S2.sam.body > $SYNTHDIR/$S2.sam
cat $S3.header $SYNTHDIR/$S3.sam.body > $SYNTHDIR/$S3.sam
cat $S4.header $SYNTHDIR/$S4.sam.body > $SYNTHDIR/$S4.sam

#8. Convert filtered.sam to BAM format
samtools view -b $SYNTHDIR/$S1.sam > $SYNTHDIR/$S1.bam
samtools view -b $SYNTHDIR/$S2.sam > $SYNTHDIR/$S2.bam
samtools view -b $SYNTHDIR/$S3.sam > $SYNTHDIR/$S3.bam
samtools view -b $SYNTHDIR/$S4.sam > $SYNTHDIR/$S4.bam


#9. Index .bam files and remove .sam
parallel  samtools index ::: *.bam
rm *.sam

#10. Produce simulated pools with synth_pool.py
module  load python/3.8.10
    PYTHONENV0=/lustre03/project/6070393/COMMON/samplepooling/software/ENVscRNA
    source $PYTHONENV0/bin/activate
    export PYTHONPATH=$PYTHONENV0/lib/python3.8/site-packages
    
$PYTHONENV0/bin/python3.8 $SOFTS/synth_pool.py \
-s $SYNTHDIR/$S1.bam,$SYNTHDIR/$S2.bam,$SYNTHDIR/$S3.bam,$SYNTHDIR/$S4.bam \
-b $SYNTHDIR/BC1.tsv,$SYNTHDIR/BC2.tsv,$SYNTHDIR/BC3.tsv,$SYNTHDIR/BC4.tsv \
-d 0.15 \
-r $REF \
-o $SYNTHDIR
