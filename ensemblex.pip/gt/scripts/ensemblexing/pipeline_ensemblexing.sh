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
echo "------Parameters used in this step---------"
echo "* PAR_ensemblex_sample_size   $PAR_ensemblex_sample_size"
echo "* PAR_ensemblex_expected_doublet_rate   $PAR_ensemblex_expected_doublet_rate"
echo "* PAR_ensemblex_merge_constituents   $PAR_ensemblex_merge_constituents"
echo "* PAR_ensemblex_probabilistic_weighted_ensemble   $PAR_ensemblex_probabilistic_weighted_ensemble"
echo "* PAR_ensemblex_preliminary_parameter_sweep   $PAR_ensemblex_preliminary_parameter_sweep"
echo "* PAR_ensemblex_graph_based_doublet_detection   $PAR_ensemblex_graph_based_doublet_detection"
echo "* PAR_ensemblex_preliminary_ensemble_independent_doublet   $PAR_ensemblex_preliminary_ensemble_independent_doublet"
echo "* PAR_ensemblex_ensemble_independent_doublet   $PAR_ensemblex_ensemble_independent_doublet"
echo "* PAR_ensemblex_doublet_Demuxalot_threshold   $PAR_ensemblex_doublet_Demuxalot_threshold"
echo "* PAR_ensemblex_doublet_Demuxalot_no_threshold   $PAR_ensemblex_doublet_Demuxalot_no_threshold"
echo "* PAR_ensemblex_doublet_Demuxlet_threshold   $PAR_ensemblex_doublet_Demuxlet_threshold"
echo "* PAR_ensemblex_doublet_Demuxlet_no_threshold   $PAR_ensemblex_doublet_Demuxlet_no_threshold"
echo "* PAR_ensemblex_doublet_Souporcell_threshold   $PAR_ensemblex_doublet_Souporcell_threshold"
echo "* PAR_ensemblex_doublet_Souporcell_no_threshold   $PAR_ensemblex_doublet_Souporcell_no_threshold"
echo "* PAR_ensemblex_doublet_Vireo_threshold   $PAR_ensemblex_doublet_Vireo_threshold"
echo "* PAR_ensemblex_doublet_Vireo_no_threshold   $PAR_ensemblex_doublet_Vireo_no_threshold"
echo "* PAR_ensemblex_compute_singlet_confidence   $PAR_ensemblex_compute_singlet_confidence"
echo "* PAR_ensemblex_nCD      $PAR_ensemblex_nCD"
echo "* PAR_ensemblex_pT      $PAR_ensemblex_pT"


echo "-------------------------------------------"
echo -e "------Output of Run------------------------\n\n"
CONTAINER1=$PIPELINE_HOME/soft/ensemblex.sif

#----------------------------------------------------------------#
# START PIPELINE                                      #
#----------------------------------------------------------------#
echo "Start of ensemblexing"
$CONTAINER_CMD exec  --bind  $OUTPUT_DIR,$PIPELINE_HOME ${CONTAINER1} Rscript $PIPELINE_HOME/gt/scripts/ensemblexing/ensemblexing.R $PIPELINE_HOME $OUTPUT_DIR $PAR_ensemblex_sample_size $PAR_ensemblex_expected_doublet_rate $PAR_ensemblex_merge_constituents $PAR_ensemblex_probabilistic_weighted_ensemble $PAR_ensemblex_preliminary_parameter_sweep $PAR_ensemblex_graph_based_doublet_detection $PAR_ensemblex_preliminary_ensemble_independent_doublet $PAR_ensemblex_ensemble_independent_doublet $PAR_ensemblex_doublet_Demuxalot_threshold $PAR_ensemblex_doublet_Demuxalot_no_threshold $PAR_ensemblex_doublet_Demuxlet_threshold $PAR_ensemblex_doublet_Demuxlet_no_threshold $PAR_ensemblex_doublet_Souporcell_threshold $PAR_ensemblex_doublet_Souporcell_no_threshold $PAR_ensemblex_doublet_Vireo_threshold $PAR_ensemblex_doublet_Vireo_no_threshold $PAR_ensemblex_compute_singlet_confidence $PAR_ensemblex_nCD $PAR_ensemblex_pT
echo "End of ensemblexing"
exit 0
