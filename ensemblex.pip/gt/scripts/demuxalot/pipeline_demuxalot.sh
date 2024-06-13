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
echo "* step Swmuxalot GT submitted at `date +%FT%H.%M.%S`"
echo "-------------------------------------------"
echo "* PIPELINE_HOME:           $PIPELINE_HOME"
echo "* OUTPUT_DIR:              $OUTPUT_DIR"
echo "-------------------------------------------"
echo "------Parameters used in this step---------"
echo "* PAR_demuxalot_genotype_names  $PAR_demuxalot_genotype_names"
echo "* PAR_demuxalot_prior_strength  $PAR_demuxalot_prior_strength"
echo "* PAR_demuxalot_minimum_coverage  $PAR_demuxalot_minimum_coverage"
echo "* PAR_demuxalot_minimum_alternative_coverage  $PAR_demuxalot_minimum_alternative_coverage"
echo "* PAR_demuxalot_n_best_snps_per_donor  $PAR_demuxalot_n_best_snps_per_donor"
echo "* PAR_demuxalot_genotypes_prior_strength  $PAR_demuxalot_genotypes_prior_strength"
echo "* PAR_demuxalot_doublet_prior  $PAR_demuxalot_doublet_prior"
echo "-------------------------------------------"
echo -e "------Output of Run------------------------\n\n"
CONTAINER1=$PIPELINE_HOME/soft/ensemblex.sif
#----------------------------------------------------------------#
# START PIPELINE                                      #
#----------------------------------------------------------------#
echo "Start of demuxalot"
$CONTAINER_CMD exec  --bind  $OUTPUT_DIR,$PIPELINE_HOME ${CONTAINER1} python3 $PIPELINE_HOME/gt/scripts/demuxalot/pipline_demuxalot.py -fl $OUTPUT_DIR -p1 $PAR_demuxalot_genotype_names -p2 $PAR_demuxalot_prior_strength -p3 $PAR_demuxalot_minimum_coverage -p4 $PAR_demuxalot_minimum_alternative_coverage -p5 $PAR_demuxalot_n_best_snps_per_donor -p6 $PAR_demuxalot_genotypes_prior_strength -p7 $PAR_demuxalot_doublet_prior
echo "End of demuxalot"
exit 0
