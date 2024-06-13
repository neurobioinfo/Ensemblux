#!/bin/bash

source $OUTPUT_DIR/job_info/configs/ensemblex_config.ini
source $OUTPUT_DIR/job_info/.tmp/temp_config.ini
export QUEUE=bash

MODE0=$MODE

# ===============================================
# STEP Vireo: 
# ===============================================
#
STEP=vireo
if [[  ${MODE0[@]}  =~  vireo ]] || [[  ${MODE0[@]}  =~  all ]]  ; then
  echo -e "\n\n-----------------------------------------------------------" >> $EXPECTED_DONE_FILES
  echo -e  "--------Job submitted using pipeline version $VERSION--------"  >> $EXPECTED_DONE_FILES
  echo "* It submitted at `date +%FT%H.%M.%S`" >> $EXPECTED_DONE_FILES
  echo "-----------------------------------------------------------"  >> $EXPECTED_DONE_FILES
  echo "STEP  submitted $STEP GT "  >> $EXPECTED_DONE_FILES
  $QUEUE $PIPELINE_HOME/gt/scripts/vireo/pipeline_vireo.sh
  sleep 10
  $QUEUE $JOB_OUTPUT_DIR/.tmp/pipeline_vireo.sh &> $JOB_OUTPUT_DIR/logs/${STEP}.$(date +%FT%H.%M.%S) 
  echo "-------------------------------------------" >> $EXPECTED_DONE_FILES
  echo "The output is under ${OUTPUT_DIR}/vireo" >> $EXPECTED_DONE_FILES  
fi 

# ===============================================
# STEP souporcell: 
# ===============================================
#
STEP=souporcell
if [[  ${MODE0[@]}  =~  souporcell ]] || [[  ${MODE0[@]}  =~  all ]]  ; then
  echo -e "\n\n-----------------------------------------------------------" >> $EXPECTED_DONE_FILES
  echo -e  "--------Job submitted using pipeline version $VERSION--------"  >> $EXPECTED_DONE_FILES
  echo "* It submitted at `date +%FT%H.%M.%S`" >> $EXPECTED_DONE_FILES
  echo "-----------------------------------------------------------"  >> $EXPECTED_DONE_FILES
  echo "STEP $STEP  GT submitted"  >> $EXPECTED_DONE_FILES
  $QUEUE $PIPELINE_HOME/gt/scripts/souporcell/pipeline_souporcell_generate.sh
  sleep 10
  $QUEUE $JOB_OUTPUT_DIR/.tmp/pipeline_souporcell.sh &> $JOB_OUTPUT_DIR/logs/${STEP}.$(date +%FT%H.%M.%S) 
  echo "-------------------------------------------" >> $EXPECTED_DONE_FILES
  echo "The output is under ${OUTPUT_DIR}/souporcell" >> $EXPECTED_DONE_FILES  
fi 


# ===============================================
# STEP Demuxlet: 
# ===============================================
#
STEP=demuxlet
if [[  ${MODE0[@]}  =~  demuxlet ]] || [[  ${MODE0[@]}  =~  all ]] ; then
  echo -e "\n\n-----------------------------------------------------------" >> $EXPECTED_DONE_FILES
  echo -e  "--------Job submitted using pipeline version $VERSION--------"  >> $EXPECTED_DONE_FILES
  echo "* It submitted at `date +%FT%H.%M.%S`" >> $EXPECTED_DONE_FILES
  echo "-----------------------------------------------------------"  >> $EXPECTED_DONE_FILES
  echo "STEP $STEP GT submitted"  >> $EXPECTED_DONE_FILES

  $QUEUE $PIPELINE_HOME/gt/scripts/demuxlet/pipeline_demuxlet.sh &> $JOB_OUTPUT_DIR/logs/${STEP}.$(date +%FT%H.%M.%S) 
  echo "-------------------------------------------" >> $EXPECTED_DONE_FILES
  echo "The output is under ${OUTPUT_DIR}/demuxlet" >> $EXPECTED_DONE_FILES  
fi 

# ===============================================
# STEP Demuxalot: 
# ===============================================
#
STEP=demuxalot
if [[  ${MODE0[@]}  =~  demuxalot ]] || [[  ${MODE0[@]}  =~  all ]] ; then
  echo -e "\n\n-----------------------------------------------------------" >> $EXPECTED_DONE_FILES
  echo -e  "--------Job submitted using pipeline version $VERSION--------"  >> $EXPECTED_DONE_FILES
  echo "* It submitted at `date +%FT%H.%M.%S`" >> $EXPECTED_DONE_FILES
  echo "-----------------------------------------------------------"  >> $EXPECTED_DONE_FILES
  echo "STEP $STEP GT submitted"  >> $EXPECTED_DONE_FILES

  $QUEUE $PIPELINE_HOME/gt/scripts/demuxalot/pipeline_demuxalot.sh &> $JOB_OUTPUT_DIR/logs/${STEP}.$(date +%FT%H.%M.%S) 
  echo "-------------------------------------------" >> $EXPECTED_DONE_FILES
  echo "The output is under ${OUTPUT_DIR}/demuxalot" >> $EXPECTED_DONE_FILES  
fi 


# ===============================================
# STEP ensemblexing: 
# ===============================================
#

STEP=ensemblexing
if [[  ${MODE0[@]}  =~  ensemblexing ]] || [[  ${MODE0[@]}  =~  all ]] ; then

  echo -e "\n\n-----------------------------------------------------------" >> $EXPECTED_DONE_FILES
  echo -e  "--------Job submitted using pipeline version $VERSION--------"  >> $EXPECTED_DONE_FILES
  echo "* It submitted at `date +%FT%H.%M.%S`" >> $EXPECTED_DONE_FILES
  echo "-----------------------------------------------------------"  >> $EXPECTED_DONE_FILES
  echo "STEP $STEP GT submitted"  >> $EXPECTED_DONE_FILES

  $QUEUE $PIPELINE_HOME/gt/scripts/ensemblexing/pipeline_ensemblexing.sh &> $JOB_OUTPUT_DIR/logs/${STEP}.$(date +%FT%H.%M.%S) 
  echo "-------------------------------------------" >> $EXPECTED_DONE_FILES
  echo "The output is under ${OUTPUT_DIR}/ensemblexing" >> $EXPECTED_DONE_FILES  
fi 

exit 0



