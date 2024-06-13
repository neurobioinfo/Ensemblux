#!/bin/bash

umask 002
source $OUTPUT_DIR/job_info/configs/ensemblex_config.ini
source $OUTPUT_DIR/job_info/.tmp/temp_config.ini

#----------------------------------------------------------------#
#                                                     #
# INITIALIZE VARIABLES                                #
#                                                     #
#----------------------------------------------------------------#
echo "-------------------------------------------"
echo "* step Souporcell submitted at `date +%FT%H.%M.%S`"
echo "-------------------------------------------"
echo "* PIPELINE_HOME:           $PIPELINE_HOME"
echo "* OUTPUT_DIR:              $OUTPUT_DIR"
echo "-------------------------------------------"
echo "------Parameters used in this step---------"
echo "* PAR_demuxlet_field:             $PAR_demuxlet_field"
echo "-------------------------------------------"
echo -e "------Output of Run------------------------\n\n"
CONTAINER1=$PIPELINE_HOME/soft/ensemblex.sif
# SOFT_SOUP=/opt/souporcell
#----------------------------------------------------------------#
# START PIPELINE                                      #
#----------------------------------------------------------------#
echo "Start of pileup step"
$CONTAINER_CMD exec  --bind  $OUTPUT_DIR ${CONTAINER1}  /opt/popscle/bin/popscle dsc-pileup \
--sam $OUTPUT_DIR/input_files/pooled_bam.bam \
--vcf $OUTPUT_DIR/input_files/pooled_samples.vcf \
--group-list $OUTPUT_DIR/input_files/pooled_barcodes.tsv \
--out $OUTPUT_DIR/demuxlet/pileup
echo "End of pileup step"

echo "Start of demuxlet step"
$CONTAINER_CMD exec  --bind  $OUTPUT_DIR ${CONTAINER1}  /opt/popscle/bin/popscle demuxlet  \
--plp $OUTPUT_DIR/demuxlet/pileup \
--vcf $OUTPUT_DIR/input_files/pooled_samples.vcf \
--field $PAR_demuxlet_field \
--group-list $OUTPUT_DIR/input_files/pooled_barcodes.tsv \
--out $OUTPUT_DIR/demuxlet/outs
echo "End of demuxlet step"

exit 0
