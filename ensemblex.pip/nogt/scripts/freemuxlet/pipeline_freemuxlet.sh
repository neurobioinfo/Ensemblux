#!/bin/bash

umask 002
source $OUTPUT_DIR/job_info/configs/ensemblex_config.ini
source $OUTPUT_DIR/job_info/.tmp/temp_config.ini

if [ -d $OUTPUT_DIR/ensemblex ]; then
   rm -rf  $OUTPUT_DIR/freemuxlet
   mkdir -p $OUTPUT_DIR/freemuxlet
else 
   mkdir -p $OUTPUT_DIR/freemuxlet
fi

#----------------------------------------------------------------#
#                                                     #
# INITIALIZE VARIABLES                                #
#                                                     #
#----------------------------------------------------------------#
echo "-------------------------------------------"
echo "* step Freemuxlet No-GT submitted at `date +%FT%H.%M.%S`"
echo "-------------------------------------------"
echo "* PIPELINE_HOME:           $PIPELINE_HOME"
echo "* OUTPUT_DIR:              $OUTPUT_DIR"
echo "-------------------------------------------"
echo "------Parameters used in this step---------"
echo "* PAR_freemuxlet_field:             $PAR_freemuxlet_field"
echo "-------------------------------------------"
echo -e "------Output of Run------------------------\n\n"
CONTAINER1=$PIPELINE_HOME/soft/ensemblex.sif

#----------------------------------------------------------------#
# START PIPELINE                                      #
#----------------------------------------------------------------#
echo "Start of pileup step"
$CONTAINER_CMD exec  --bind  $OUTPUT_DIR ${CONTAINER1}  /opt/popscle/bin/popscle dsc-pileup \
--sam $OUTPUT_DIR/input_files/pooled_bam.bam \
--vcf $OUTPUT_DIR/input_files/reference.vcf \
--group-list $OUTPUT_DIR/input_files/pooled_barcodes.tsv \
--out $OUTPUT_DIR/freemuxlet/pileup

echo "End of Step pileup"

echo "Start of Freemuxlet step"
$CONTAINER_CMD exec  --bind  $OUTPUT_DIR ${CONTAINER1}  /opt/popscle/bin/popscle freemuxlet  \
--plp $OUTPUT_DIR/freemuxlet/pileup \
--nsample  $PAR_freemuxlet_nsample \
--group-list $OUTPUT_DIR/input_files/pooled_barcodes.tsv \
--out $OUTPUT_DIR/freemuxlet/outs

echo "End of Freemuxlet step"
exit 0
