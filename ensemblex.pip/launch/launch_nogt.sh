#!/bin/bash

source $OUTPUT_DIR/job_info/configs/ensemblex_config.ini
source $OUTPUT_DIR/job_info/.tmp/temp_config.ini

export QUEUE=bash

MODE0=$MODE



# ===============================================
# STEP Demuxlet: 
# ===============================================
#
STEP=freemuxlet
if [[  ${MODE0[@]}  =~  freemuxlet ]] || [[  ${MODE0[@]}  =~  all ]] ; then
  echo -e "\n\n-----------------------------------------------------------" >> $EXPECTED_DONE_FILES
  echo -e  "--------Job submitted using pipeline version $VERSION--------"  >> $EXPECTED_DONE_FILES
  echo "* It submitted at `date +%FT%H.%M.%S`" >> $EXPECTED_DONE_FILES
  echo "-----------------------------------------------------------"  >> $EXPECTED_DONE_FILES
  echo "STEP $STEP freemuxlet submitted"  >> $EXPECTED_DONE_FILES

  $QUEUE $PIPELINE_HOME/nogt/scripts/freemuxlet/pipeline_freemuxlet.sh &> $JOB_OUTPUT_DIR/logs/${STEP}.$(date +%FT%H.%M.%S) 
  echo "-------------------------------------------" >> $EXPECTED_DONE_FILES
  echo "The output is under ${OUTPUT_DIR}/freemuxlet" >> $EXPECTED_DONE_FILES  
fi

# ===============================================
# STEP Vireo: 
# ===============================================
#
STEP=vireo
if [[  ${MODE0[@]}  =~  vireo ]] || [[  ${MODE0[@]}  =~  all ]] ; then
  echo -e "\n\n-----------------------------------------------------------" >> $EXPECTED_DONE_FILES
  echo -e  "--------Job submitted using pipeline version $VERSION--------"  >> $EXPECTED_DONE_FILES
  echo "* It submitted at `date +%FT%H.%M.%S`" >> $EXPECTED_DONE_FILES
  echo "-----------------------------------------------------------"  >> $EXPECTED_DONE_FILES
  echo "STEP  submitted $STEP No-GT "  >> $EXPECTED_DONE_FILES
  $QUEUE $PIPELINE_HOME/nogt/scripts/vireo/pipeline_vireo.sh
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
if [[  ${MODE0[@]}  =~  souporcell ]] || [[  ${MODE0[@]}  =~  all ]] ; then
  echo -e "\n\n-----------------------------------------------------------" >> $EXPECTED_DONE_FILES
  echo -e  "--------Job submitted using pipeline version $VERSION--------"  >> $EXPECTED_DONE_FILES
  echo "* It submitted at `date +%FT%H.%M.%S`" >> $EXPECTED_DONE_FILES
  echo "-----------------------------------------------------------"  >> $EXPECTED_DONE_FILES
  echo "STEP $STEP No-GT submitted"  >> $EXPECTED_DONE_FILES
  $QUEUE $PIPELINE_HOME/nogt/scripts/souporcell/pipeline_souporcell_generate.sh
  sleep 10
  $QUEUE $JOB_OUTPUT_DIR/.tmp/pipeline_souporcell.sh &> $JOB_OUTPUT_DIR/logs/${STEP}.$(date +%FT%H.%M.%S) 
  echo "-------------------------------------------" >> $EXPECTED_DONE_FILES
  echo "The output is under ${OUTPUT_DIR}/souporcell" >> $EXPECTED_DONE_FILES  
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
  echo "STEP $STEP submitted"  >> $EXPECTED_DONE_FILES

  $QUEUE $PIPELINE_HOME/nogt/scripts/demuxalot/pipeline_demuxalot.sh &> $JOB_OUTPUT_DIR/logs/${STEP}.$(date +%FT%H.%M.%S) 
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
  echo "STEP $STEP submitted"  >> $EXPECTED_DONE_FILES

  $QUEUE $PIPELINE_HOME/nogt/scripts/ensemblexing/pipeline_ensemblexing.sh &> $JOB_OUTPUT_DIR/logs/${STEP}.$(date +%FT%H.%M.%S) 
  echo "-------------------------------------------" >> $EXPECTED_DONE_FILES
  echo "The output is under ${OUTPUT_DIR}/ensemblexing" >> $EXPECTED_DONE_FILES  
fi 

exit 0
