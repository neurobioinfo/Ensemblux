#!/bin/bash

#########
# NOTES #
#########
# This code was used to produce simulated pools with known ground truth sample labels using 80 independently-sequenced iPSC lines differentiated towards a dopaminergic neuronal state as part of the Foundational Data Initiative for Parkinson's Disease
# Simulated pools were produced using the synth_pool.py script produced by the developers of Vireo. 24 samples were multiplexed per pool with 8, 16, 24, 32, 40, ad 48 thousand cells
# Below is an example script for producing a simulated pool from  24 lines, totaling ~23,100 cells and a 15% doublet rate

#################
# SET VARIABLES #
#################
POOL=8_24_1

cd /home/fiorini9/scratch/mjff/emsemblux_seq_feature_eval/sample_size/$POOL

S1=SCRN_PPMI3186_2662_da65_v1
S2=SCRN_PPMI3220_6139_da65_v1
S3=SCRN_PPMI3234_1338_da65_v1
S4=SCRN_PPMI3409_0660_da65_v1
S5=SCRN_PPMI3411_0083_da65_v1
S6=SCRN_PPMI3419_0298_da65_v1
S7=SCRN_PPMI3422_1260_da65_v1
S8=SCRN_PPMI3446_0966_da65_v1
S9=SCRN_PPMI3448_2397_da65_v1
S10=SCRN_PPMI3451_3236_da65_v1
S11=SCRN_PPMI3452_7426_da65_v1
S12=SCRN_PPMI3453_7504_da65_v1
S13=SCRN_PPMI3459_1662_da65_v1
S14=SCRN_PPMI3460_1783_da65_v1
S15=SCRN_PPMI3467_3278_da65_v1
S16=SCRN_PPMI3471_8520_da65_v1
S17=SCRN_PPMI3473_0119_da65_v1
S18=SCRN_PPMI3475_1575_da65_v1
S19=SCRN_PPMI3588_6182_da65_v1
S20=SCRN_PPMI3658_4030_da65_v1
S21=SCRN_PPMI3664_2833_da65_v1
S22=SCRN_PPMI3665_4484_da65_v1
S23=SCRN_PPMI3666_3014_da65_v1
S24=SCRN_PPMI3668_3953_da65_v1

SYNTHDIR=/home/fiorini9/scratch/mjff/emsemblux_seq_feature_eval/sample_size/$POOL 
DATA=/home/fiorini9/nearline/rrg-tdurcan/DATA/FoundIN-PD/download/FOUNDIN/processed/SCRN/SCRN/cellRanger
BAM=outs/possorted_genome_bam.bam
SOFTS=/lustre03/project/6070393/COMMON/samplepooling/MF/synthetic_data/synth_submit
REF=/lustre03/project/6070393/COMMON/samplepooling/MF/synthetic_data/VCF/reheader_rename_PDHCSWEDD_progen_merged.vcf 

BARCODES1=/home/fiorini9/nearline/rrg-tdurcan/DATA/FoundIN-PD/download/FOUNDIN/processed/SCRN/SCRN/cellRanger
BARCODES2=outs/filtered_feature_bc_matrix/barcodes.tsv.gz
DATADIR=/home/fiorini9/scratch/mjff/emsemblux_seq_feature_eval/sample_size 
BAM1=/home/fiorini9/nearline/rrg-tdurcan/DATA/FoundIN-PD/download/FOUNDIN/processed/SCRN/SCRN/cellRanger
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
cp $BARCODES1/$S5/$BARCODES2 $DATADIR/$POOL/$S5.tsv.gz
cp $BARCODES1/$S6/$BARCODES2 $DATADIR/$POOL/$S6.tsv.gz
cp $BARCODES1/$S7/$BARCODES2 $DATADIR/$POOL/$S7.tsv.gz
cp $BARCODES1/$S8/$BARCODES2 $DATADIR/$POOL/$S8.tsv.gz
cp $BARCODES1/$S9/$BARCODES2 $DATADIR/$POOL/$S9.tsv.gz
cp $BARCODES1/$S10/$BARCODES2 $DATADIR/$POOL/$S10.tsv.gz
cp $BARCODES1/$S11/$BARCODES2 $DATADIR/$POOL/$S11.tsv.gz
cp $BARCODES1/$S12/$BARCODES2 $DATADIR/$POOL/$S12.tsv.gz
cp $BARCODES1/$S13/$BARCODES2 $DATADIR/$POOL/$S13.tsv.gz
cp $BARCODES1/$S14/$BARCODES2 $DATADIR/$POOL/$S14.tsv.gz
cp $BARCODES1/$S15/$BARCODES2 $DATADIR/$POOL/$S15.tsv.gz
cp $BARCODES1/$S16/$BARCODES2 $DATADIR/$POOL/$S16.tsv.gz
cp $BARCODES1/$S17/$BARCODES2 $DATADIR/$POOL/$S17.tsv.gz
cp $BARCODES1/$S18/$BARCODES2 $DATADIR/$POOL/$S18.tsv.gz
cp $BARCODES1/$S19/$BARCODES2 $DATADIR/$POOL/$S19.tsv.gz
cp $BARCODES1/$S20/$BARCODES2 $DATADIR/$POOL/$S20.tsv.gz
cp $BARCODES1/$S21/$BARCODES2 $DATADIR/$POOL/$S21.tsv.gz
cp $BARCODES1/$S22/$BARCODES2 $DATADIR/$POOL/$S22.tsv.gz
cp $BARCODES1/$S23/$BARCODES2 $DATADIR/$POOL/$S23.tsv.gz
cp $BARCODES1/$S24/$BARCODES2 $DATADIR/$POOL/$S24.tsv.gz

gunzip *.gz

#2. copy over the bam files
cp $BAM1/$S1/$BAM2  $SYNTHDIR/$S1.bam
cp $BAM1/$S2/$BAM2  $SYNTHDIR/$S2.bam
cp $BAM1/$S3/$BAM2  $SYNTHDIR/$S3.bam
cp $BAM1/$S4/$BAM2  $SYNTHDIR/$S4.bam
cp $BAM1/$S5/$BAM2  $SYNTHDIR/$S5.bam
cp $BAM1/$S6/$BAM2  $SYNTHDIR/$S6.bam
cp $BAM1/$S7/$BAM2  $SYNTHDIR/$S7.bam
cp $BAM1/$S8/$BAM2  $SYNTHDIR/$S8.bam
cp $BAM1/$S9/$BAM2  $SYNTHDIR/$S9.bam
cp $BAM1/$S10/$BAM2  $SYNTHDIR/$S10.bam
cp $BAM1/$S11/$BAM2  $SYNTHDIR/$S11.bam
cp $BAM1/$S12/$BAM2  $SYNTHDIR/$S12.bam
cp $BAM1/$S13/$BAM2  $SYNTHDIR/$S13.bam
cp $BAM1/$S14/$BAM2  $SYNTHDIR/$S14.bam
cp $BAM1/$S15/$BAM2  $SYNTHDIR/$S15.bam
cp $BAM1/$S16/$BAM2  $SYNTHDIR/$S16.bam
cp $BAM1/$S17/$BAM2  $SYNTHDIR/$S17.bam
cp $BAM1/$S18/$BAM2  $SYNTHDIR/$S18.bam
cp $BAM1/$S19/$BAM2  $SYNTHDIR/$S19.bam
cp $BAM1/$S20/$BAM2  $SYNTHDIR/$S20.bam
cp $BAM1/$S21/$BAM2  $SYNTHDIR/$S21.bam
cp $BAM1/$S22/$BAM2  $SYNTHDIR/$S22.bam
cp $BAM1/$S23/$BAM2  $SYNTHDIR/$S23.bam
cp $BAM1/$S24/$BAM2  $SYNTHDIR/$S24.bam


#3. filter barcodes files to only include th number of cells that we want
cat $DATADIR/$POOL/$S1.tsv | awk 'NR >= 0  && NR <= 333' > BC1.tsv
cat $DATADIR/$POOL/$S2.tsv | awk 'NR >= 0  && NR <= 333' > BC2.tsv
cat $DATADIR/$POOL/$S3.tsv | awk 'NR >= 0  && NR <= 333' > BC3.tsv
cat $DATADIR/$POOL/$S4.tsv | awk 'NR >= 0  && NR <= 333' > BC4.tsv
cat $DATADIR/$POOL/$S5.tsv | awk 'NR >= 0  && NR <= 333' > BC5.tsv
cat $DATADIR/$POOL/$S6.tsv | awk 'NR >= 0  && NR <= 333' > BC6.tsv
cat $DATADIR/$POOL/$S7.tsv | awk 'NR >= 0  && NR <= 333' > BC7.tsv
cat $DATADIR/$POOL/$S8.tsv | awk 'NR >= 0  && NR <= 333' > BC8.tsv
cat $DATADIR/$POOL/$S9.tsv | awk 'NR >= 0  && NR <= 333' > BC9.tsv
cat $DATADIR/$POOL/$S10.tsv | awk 'NR >= 0  && NR <= 333' > BC10.tsv
cat $DATADIR/$POOL/$S11.tsv | awk 'NR >= 0  && NR <= 333' > BC11.tsv
cat $DATADIR/$POOL/$S12.tsv | awk 'NR >= 0  && NR <= 333' > BC12.tsv
cat $DATADIR/$POOL/$S13.tsv | awk 'NR >= 0  && NR <= 333' > BC13.tsv
cat $DATADIR/$POOL/$S14.tsv | awk 'NR >= 0  && NR <= 333' > BC14.tsv
cat $DATADIR/$POOL/$S15.tsv | awk 'NR >= 0  && NR <= 333' > BC15.tsv
cat $DATADIR/$POOL/$S16.tsv | awk 'NR >= 0  && NR <= 333' > BC16.tsv
cat $DATADIR/$POOL/$S17.tsv | awk 'NR >= 0  && NR <= 333' > BC17.tsv
cat $DATADIR/$POOL/$S18.tsv | awk 'NR >= 0  && NR <= 333' > BC18.tsv
cat $DATADIR/$POOL/$S19.tsv | awk 'NR >= 0  && NR <= 333' > BC19.tsv
cat $DATADIR/$POOL/$S20.tsv | awk 'NR >= 0  && NR <= 333' > BC20.tsv
cat $DATADIR/$POOL/$S21.tsv | awk 'NR >= 0  && NR <= 333' > BC21.tsv
cat $DATADIR/$POOL/$S22.tsv | awk 'NR >= 0  && NR <= 333' > BC22.tsv
cat $DATADIR/$POOL/$S23.tsv | awk 'NR >= 0  && NR <= 333' > BC23.tsv
cat $DATADIR/$POOL/$S24.tsv | awk 'NR >= 0  && NR <= 333' > BC24.tsv

#4. add BC tag to our new barcodes file
sed -e 's/^/CB:Z:/' BC1.tsv > tag_BC1.tsv
sed -e 's/^/CB:Z:/' BC2.tsv > tag_BC2.tsv
sed -e 's/^/CB:Z:/' BC3.tsv > tag_BC3.tsv
sed -e 's/^/CB:Z:/' BC4.tsv > tag_BC4.tsv
sed -e 's/^/CB:Z:/' BC5.tsv > tag_BC5.tsv
sed -e 's/^/CB:Z:/' BC6.tsv > tag_BC6.tsv
sed -e 's/^/CB:Z:/' BC7.tsv > tag_BC7.tsv
sed -e 's/^/CB:Z:/' BC8.tsv > tag_BC8.tsv
sed -e 's/^/CB:Z:/' BC9.tsv > tag_BC9.tsv
sed -e 's/^/CB:Z:/' BC10.tsv > tag_BC10.tsv
sed -e 's/^/CB:Z:/' BC11.tsv > tag_BC11.tsv
sed -e 's/^/CB:Z:/' BC12.tsv > tag_BC12.tsv
sed -e 's/^/CB:Z:/' BC13.tsv > tag_BC13.tsv
sed -e 's/^/CB:Z:/' BC14.tsv > tag_BC14.tsv
sed -e 's/^/CB:Z:/' BC15.tsv > tag_BC15.tsv
sed -e 's/^/CB:Z:/' BC16.tsv > tag_BC16.tsv
sed -e 's/^/CB:Z:/' BC17.tsv > tag_BC17.tsv
sed -e 's/^/CB:Z:/' BC18.tsv > tag_BC18.tsv
sed -e 's/^/CB:Z:/' BC19.tsv > tag_BC19.tsv
sed -e 's/^/CB:Z:/' BC20.tsv > tag_BC20.tsv
sed -e 's/^/CB:Z:/' BC21.tsv > tag_BC21.tsv
sed -e 's/^/CB:Z:/' BC22.tsv > tag_BC22.tsv
sed -e 's/^/CB:Z:/' BC23.tsv > tag_BC23.tsv
sed -e 's/^/CB:Z:/' BC24.tsv > tag_BC24.tsv

#5. Save the header lines
samtools view -H $SYNTHDIR/$S1.bam > $S1.header
samtools view -H $SYNTHDIR/$S2.bam > $S2.header
samtools view -H $SYNTHDIR/$S3.bam > $S3.header
samtools view -H $SYNTHDIR/$S4.bam > $S4.header
samtools view -H $SYNTHDIR/$S5.bam > $S5.header
samtools view -H $SYNTHDIR/$S6.bam > $S6.header
samtools view -H $SYNTHDIR/$S7.bam > $S7.header
samtools view -H $SYNTHDIR/$S8.bam > $S8.header
samtools view -H $SYNTHDIR/$S9.bam > $S9.header
samtools view -H $SYNTHDIR/$S10.bam > $S10.header
samtools view -H $SYNTHDIR/$S11.bam > $S11.header
samtools view -H $SYNTHDIR/$S12.bam > $S12.header
samtools view -H $SYNTHDIR/$S13.bam > $S13.header
samtools view -H $SYNTHDIR/$S14.bam > $S14.header
samtools view -H $SYNTHDIR/$S15.bam > $S15.header
samtools view -H $SYNTHDIR/$S16.bam > $S16.header
samtools view -H $SYNTHDIR/$S17.bam > $S17.header
samtools view -H $SYNTHDIR/$S18.bam > $S18.header
samtools view -H $SYNTHDIR/$S19.bam > $S19.header
samtools view -H $SYNTHDIR/$S20.bam > $S20.header
samtools view -H $SYNTHDIR/$S21.bam > $S21.header
samtools view -H $SYNTHDIR/$S22.bam > $S22.header
samtools view -H $SYNTHDIR/$S23.bam > $S23.header
samtools view -H $SYNTHDIR/$S24.bam > $S24.header

#6. Filter alignments using filter.txt. Use LC_ALL=C to set C locale instead of UTF-8
samtools view $SYNTHDIR/$S1.bam | LC_ALL=C grep -F -f tag_BC1.tsv > $SYNTHDIR/$S1.sam.body
samtools view $SYNTHDIR/$S2.bam | LC_ALL=C grep -F -f tag_BC2.tsv > $SYNTHDIR/$S2.sam.body
samtools view $SYNTHDIR/$S3.bam | LC_ALL=C grep -F -f tag_BC3.tsv > $SYNTHDIR/$S3.sam.body
samtools view $SYNTHDIR/$S4.bam | LC_ALL=C grep -F -f tag_BC4.tsv > $SYNTHDIR/$S4.sam.body
samtools view $SYNTHDIR/$S5.bam | LC_ALL=C grep -F -f tag_BC5.tsv > $SYNTHDIR/$S5.sam.body
samtools view $SYNTHDIR/$S6.bam | LC_ALL=C grep -F -f tag_BC6.tsv > $SYNTHDIR/$S6.sam.body
samtools view $SYNTHDIR/$S7.bam | LC_ALL=C grep -F -f tag_BC7.tsv > $SYNTHDIR/$S7.sam.body
samtools view $SYNTHDIR/$S8.bam | LC_ALL=C grep -F -f tag_BC8.tsv > $SYNTHDIR/$S8.sam.body
samtools view $SYNTHDIR/$S9.bam | LC_ALL=C grep -F -f tag_BC9.tsv > $SYNTHDIR/$S9.sam.body
samtools view $SYNTHDIR/$S10.bam | LC_ALL=C grep -F -f tag_BC10.tsv > $SYNTHDIR/$S10.sam.body
samtools view $SYNTHDIR/$S11.bam | LC_ALL=C grep -F -f tag_BC11.tsv > $SYNTHDIR/$S11.sam.body
samtools view $SYNTHDIR/$S12.bam | LC_ALL=C grep -F -f tag_BC12.tsv > $SYNTHDIR/$S12.sam.body
samtools view $SYNTHDIR/$S13.bam | LC_ALL=C grep -F -f tag_BC13.tsv > $SYNTHDIR/$S13.sam.body
samtools view $SYNTHDIR/$S14.bam | LC_ALL=C grep -F -f tag_BC14.tsv > $SYNTHDIR/$S14.sam.body
samtools view $SYNTHDIR/$S15.bam | LC_ALL=C grep -F -f tag_BC15.tsv > $SYNTHDIR/$S15.sam.body
samtools view $SYNTHDIR/$S16.bam | LC_ALL=C grep -F -f tag_BC16.tsv > $SYNTHDIR/$S16.sam.body
samtools view $SYNTHDIR/$S17.bam | LC_ALL=C grep -F -f tag_BC17.tsv > $SYNTHDIR/$S17.sam.body
samtools view $SYNTHDIR/$S18.bam | LC_ALL=C grep -F -f tag_BC18.tsv > $SYNTHDIR/$S18.sam.body
samtools view $SYNTHDIR/$S19.bam | LC_ALL=C grep -F -f tag_BC19.tsv > $SYNTHDIR/$S19.sam.body
samtools view $SYNTHDIR/$S20.bam | LC_ALL=C grep -F -f tag_BC20.tsv > $SYNTHDIR/$S20.sam.body
samtools view $SYNTHDIR/$S21.bam | LC_ALL=C grep -F -f tag_BC21.tsv > $SYNTHDIR/$S21.sam.body
samtools view $SYNTHDIR/$S22.bam | LC_ALL=C grep -F -f tag_BC22.tsv > $SYNTHDIR/$S22.sam.body
samtools view $SYNTHDIR/$S23.bam | LC_ALL=C grep -F -f tag_BC23.tsv > $SYNTHDIR/$S23.sam.body
samtools view $SYNTHDIR/$S24.bam | LC_ALL=C grep -F -f tag_BC24.tsv > $SYNTHDIR/$S24.sam.body

#7. Combine header and body
cat $S1.header $SYNTHDIR/$S1.sam.body > $SYNTHDIR/$S1.sam
cat $S2.header $SYNTHDIR/$S2.sam.body > $SYNTHDIR/$S2.sam
cat $S3.header $SYNTHDIR/$S3.sam.body > $SYNTHDIR/$S3.sam
cat $S4.header $SYNTHDIR/$S4.sam.body > $SYNTHDIR/$S4.sam
cat $S5.header $SYNTHDIR/$S5.sam.body > $SYNTHDIR/$S5.sam
cat $S6.header $SYNTHDIR/$S6.sam.body > $SYNTHDIR/$S6.sam
cat $S7.header $SYNTHDIR/$S7.sam.body > $SYNTHDIR/$S7.sam
cat $S8.header $SYNTHDIR/$S8.sam.body > $SYNTHDIR/$S8.sam
cat $S9.header $SYNTHDIR/$S9.sam.body > $SYNTHDIR/$S9.sam
cat $S10.header $SYNTHDIR/$S10.sam.body > $SYNTHDIR/$S10.sam
cat $S11.header $SYNTHDIR/$S11.sam.body > $SYNTHDIR/$S11.sam
cat $S12.header $SYNTHDIR/$S12.sam.body > $SYNTHDIR/$S12.sam
cat $S13.header $SYNTHDIR/$S13.sam.body > $SYNTHDIR/$S13.sam
cat $S14.header $SYNTHDIR/$S14.sam.body > $SYNTHDIR/$S14.sam
cat $S15.header $SYNTHDIR/$S15.sam.body > $SYNTHDIR/$S15.sam
cat $S16.header $SYNTHDIR/$S16.sam.body > $SYNTHDIR/$S16.sam
cat $S17.header $SYNTHDIR/$S17.sam.body > $SYNTHDIR/$S17.sam
cat $S18.header $SYNTHDIR/$S18.sam.body > $SYNTHDIR/$S18.sam
cat $S19.header $SYNTHDIR/$S19.sam.body > $SYNTHDIR/$S19.sam
cat $S20.header $SYNTHDIR/$S20.sam.body > $SYNTHDIR/$S20.sam
cat $S21.header $SYNTHDIR/$S21.sam.body > $SYNTHDIR/$S21.sam
cat $S22.header $SYNTHDIR/$S22.sam.body > $SYNTHDIR/$S22.sam
cat $S23.header $SYNTHDIR/$S23.sam.body > $SYNTHDIR/$S23.sam
cat $S24.header $SYNTHDIR/$S24.sam.body > $SYNTHDIR/$S24.sam

#8. Convert filtered.sam to BAM format
samtools view -b $SYNTHDIR/$S1.sam > $SYNTHDIR/$S1.bam
samtools view -b $SYNTHDIR/$S2.sam > $SYNTHDIR/$S2.bam
samtools view -b $SYNTHDIR/$S3.sam > $SYNTHDIR/$S3.bam
samtools view -b $SYNTHDIR/$S4.sam > $SYNTHDIR/$S4.bam
samtools view -b $SYNTHDIR/$S5.sam > $SYNTHDIR/$S5.bam
samtools view -b $SYNTHDIR/$S6.sam > $SYNTHDIR/$S6.bam
samtools view -b $SYNTHDIR/$S7.sam > $SYNTHDIR/$S7.bam
samtools view -b $SYNTHDIR/$S8.sam > $SYNTHDIR/$S8.bam
samtools view -b $SYNTHDIR/$S9.sam > $SYNTHDIR/$S9.bam
samtools view -b $SYNTHDIR/$S10.sam > $SYNTHDIR/$S10.bam
samtools view -b $SYNTHDIR/$S11.sam > $SYNTHDIR/$S11.bam
samtools view -b $SYNTHDIR/$S12.sam > $SYNTHDIR/$S12.bam
samtools view -b $SYNTHDIR/$S13.sam > $SYNTHDIR/$S13.bam
samtools view -b $SYNTHDIR/$S14.sam > $SYNTHDIR/$S14.bam
samtools view -b $SYNTHDIR/$S15.sam > $SYNTHDIR/$S15.bam
samtools view -b $SYNTHDIR/$S16.sam > $SYNTHDIR/$S16.bam
samtools view -b $SYNTHDIR/$S17.sam > $SYNTHDIR/$S17.bam
samtools view -b $SYNTHDIR/$S18.sam > $SYNTHDIR/$S18.bam
samtools view -b $SYNTHDIR/$S19.sam > $SYNTHDIR/$S19.bam
samtools view -b $SYNTHDIR/$S20.sam > $SYNTHDIR/$S20.bam
samtools view -b $SYNTHDIR/$S21.sam > $SYNTHDIR/$S21.bam
samtools view -b $SYNTHDIR/$S22.sam > $SYNTHDIR/$S22.bam
samtools view -b $SYNTHDIR/$S23.sam > $SYNTHDIR/$S23.bam
samtools view -b $SYNTHDIR/$S24.sam > $SYNTHDIR/$S24.bam

#9. Need to index the files
parallel samtools index ::: *.bam
rm *.sam

#10. Produce simulated pools with synth_pool.py
module  load python/3.8.10
    PYTHONENV0=/lustre03/project/6070393/COMMON/samplepooling/software/ENVscRNA
    source $PYTHONENV0/bin/activate
    export PYTHONPATH=$PYTHONENV0/lib/python3.8/site-packages

cd $SYNTHDIR

$PYTHONENV0/bin/python3.8 $SOFTS/synth_pool.py \
-s $SYNTHDIR/$S1.bam,$SYNTHDIR/$S2.bam,$SYNTHDIR/$S3.bam,$SYNTHDIR/$S4.bam,$SYNTHDIR/$S5.bam,$SYNTHDIR/$S6.bam,$SYNTHDIR/$S7.bam,$SYNTHDIR/$S8.bam,$SYNTHDIR/$S9.bam,$SYNTHDIR/$S10.bam,$SYNTHDIR/$S11.bam,$SYNTHDIR/$S12.bam,$SYNTHDIR/$S13.bam,$SYNTHDIR/$S14.bam,$SYNTHDIR/$S15.bam,$SYNTHDIR/$S16.bam,$SYNTHDIR/$S17.bam,$SYNTHDIR/$S18.bam,$SYNTHDIR/$S19.bam,$SYNTHDIR/$S20.bam,$SYNTHDIR/$S21.bam,$SYNTHDIR/$S22.bam,$SYNTHDIR/$S23.bam,$SYNTHDIR/$S24.bam \
-b $SYNTHDIR/BC1.tsv,$SYNTHDIR/BC2.tsv,$SYNTHDIR/BC3.tsv,$SYNTHDIR/BC4.tsv,$SYNTHDIR/BC5.tsv,$SYNTHDIR/BC6.tsv,$SYNTHDIR/BC7.tsv,$SYNTHDIR/BC8.tsv,$SYNTHDIR/BC9.tsv,$SYNTHDIR/BC10.tsv,$SYNTHDIR/BC11.tsv,$SYNTHDIR/BC12.tsv,$SYNTHDIR/BC13.tsv,$SYNTHDIR/BC14.tsv,$SYNTHDIR/BC15.tsv,$SYNTHDIR/BC16.tsv,$SYNTHDIR/BC17.tsv,$SYNTHDIR/BC18.tsv,$SYNTHDIR/BC19.tsv,$SYNTHDIR/BC20.tsv,$SYNTHDIR/BC21.tsv,$SYNTHDIR/BC22.tsv,$SYNTHDIR/BC23.tsv,$SYNTHDIR/BC24.tsv \
-d 0.06 \
-r $REF \
-o $SYNTHDIR