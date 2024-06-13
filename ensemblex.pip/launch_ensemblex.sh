#!/bin/bash

# The Ensemblex pipeline was produced for projects funded by the Canadian Institute of Health Research and Michael J. Fox Foundation Parkinson's Progression Markers Initiative (MJFF PPMI) in collaboration with The Neuro's Early Drug Discovery Unit (EDDU), McGill University.
# Copyright belong MNI BIOINFO CORE (https://github.com/neurobioinfo)
# The pipeline is scripted by Saeid Amiri (saeid.amiri@mcgill.ca) and  Michael Fiorini( michael.fiorini@mail.mcgill.ca )

VERSION=0.0.02
DATE0=2024-06-13
echo -e "ensemblex pipeline version $VERSION"

# ===============================================
# default variables values
# ===============================================
unset OUTPUT_DIR PIPELINE_HOME

PIPELINE_HOME0=`realpath ${BASH_SOURCE[0]}`
export PIPELINE_HOME=$(cd $(dirname $PIPELINE_HOME0) && pwd -P)

TIMESTAMP=`date +%FT%H.%M.%S`


# create function to handle error messages
# ===============================================
Usage() {
	echo
  echo "------------------- " 
	echo -e "Usage:\t$0 [arguments]"
	echo -e "\tmandatory arguments:\n" \
          "\t\t-d  (--dir)  = Working directory (where all the outputs will be printed) (give full path)\n" \
          "\t\t--steps  =  Specify the steps to execute. Begin by selecting either init-GT or init-noGT to establish the working directory. \n" \
          "\t\t       For GT: vireo, demuxalot, demuxlet, souporcell, ensemblexing \n" \
          "\t\t       For noGT: vireo, demuxalot, freemuxlet, souporcell, ensemblexing \n" 
	echo -e "\toptional arguments:\n " \
          "\t\t-h  (--help)  = See helps regarding the pipeline arguments \n" \
          "\t\t--vcf  = The path of vcf file \n" \
          "\t\t--bam  = The path of bam file \n" \
          "\t\t--sortout  = The path snd nsme of vcf generated using sort  \n" \
          "------------------- \n" \
          "For a comprehensive help, visit  https://neurobioinfo.github.io/ensemblex/site/ for documentation. "

echo 
}


# ===============================================
# PARSING ARGUMENTS
# ===============================================
if ! options=$(getopt --name pipeline --alternative --unquoted --options hs:d:t:m:vw:f:S:c:a:x: --longoptions dir:,steps:,method:,container:,vcf:,bam:,sortout: -- "$@")
then
    # something went wrong, getopt will put out an error message for us
    echo "Error processing options."
    exit 42
fi

# ===============================================
# LOAD & OVERRIDE EXTRA CONFIG FILE FOR PROJECT
# ===============================================
set -- $options

while [ $# -gt 0 ]

do
    case $1 in
    -x| --extra) 
      EXTRA_CONF="$2" ;
      if [ -f $EXTRA_CONF ]; then
        echo "* LOADING EXTRA CONFIG FILE $EXTRA_CONF";
        . $EXTRA_CONF
      else
        echo "ERROR: invalid EXTRA CONFIG file: $EXTRA_CONF";
        echo "Please check options and try again"; exit 42;
      fi
    esac
    shift
done

# ===============================================
# LOAD ALL OTHER OPTIONS
# ===============================================
set -- $options

while [ $# -gt 0 ]
do
    case $1 in
    -h| --help) Usage; exit 0;;
    -d| --dir) OUTPUT_DIR="$2" ; shift ;;
    -v| --verbose) VERBOSE=1 ;; 
    --steps) MODE="$2"; shift ;; 
    --vcf) VCF="$2"; shift ;;
    --bam) BAM="$2"; shift ;; 
    --sortout) SORTOUT="$2"; shift ;; 
    --method) METHOD0="$2"; shift ;; 
    --container) CONTAINER0="$2"; shift ;;
    (--) shift; break;;
    (-*) echo "$0: error - unrecognized option $1" 1>&2; exit 42;;
    (*) break;;
    esac
    shift
done


MODE0=$MODE
FOUND_ERROR=0

# ===============================================
# CHECKING VARIABLES
# ===============================================
#check to ensure all mandatory arguments have been entered

if [ -z $OUTPUT_DIR ]; then echo "ERROR: missing mandatory option: -d (--dir) must be specified"; FOUND_ERROR=1; fi
if (( $FOUND_ERROR )); then echo "Please check options and try again"; exit 42; fi


# STEP 0: RUN setting 
# ===============================================
  JOB_OUTPUT_DIR=$OUTPUT_DIR/job_info


if [[ ${MODE0[@]} =~ init-GT ]]; then
  if [[ -s $OUTPUT_DIR/job_info ]]; then
    read -p "folder\files already exist in $OUTPUT_DIR. Overwrite? y|n: "$'\n'  answer
    if [[ $answer =~ y ]]; then
      echo "NOTE: the folder\files are Overwritten."
            rm -rf $OUTPUT_DIR/job_info; mkdir -p $OUTPUT_DIR/job_info
            rm -rf $OUTPUT_DIR/job_info/configs; mkdir -p $OUTPUT_DIR/job_info/configs
            cp -r $PIPELINE_HOME/gt/configs $OUTPUT_DIR/job_info/
            rm -rf $OUTPUT_DIR/job_info/logs; mkdir -p $OUTPUT_DIR/job_info/logs
            rm -rf $OUTPUT_DIR/input_files; mkdir -p $OUTPUT_DIR/input_files
            rm -rf $OUTPUT_DIR/demuxalot; mkdir -p $OUTPUT_DIR/demuxalot
            rm -rf $OUTPUT_DIR/demuxlet; mkdir -p $OUTPUT_DIR/demuxlet
            rm -rf $OUTPUT_DIR/ensemblex_gt; mkdir -p $OUTPUT_DIR/ensemblex_gt
            rm -rf $OUTPUT_DIR/souporcell; mkdir -p $OUTPUT_DIR/souporcell
            rm -rf $OUTPUT_DIR/vireo_gt ; mkdir -p $OUTPUT_DIR/vireo_gt
            rm -rf $OUTPUT_DIR/job_info/.tmp ; mkdir -p $OUTPUT_DIR/job_info/.tmp
    else
     echo "NOTE: the pipeline is using the existing config file and parameters."
    fi
  else
    echo "The configuration files do not exist, the pipeline will create them during execution."
        mkdir -p $OUTPUT_DIR/job_info
        mkdir -p $OUTPUT_DIR/input_files
        mkdir -p $OUTPUT_DIR/job_info/configs
        cp -r $PIPELINE_HOME/gt/configs $OUTPUT_DIR/job_info/
        mkdir -p $OUTPUT_DIR/job_info/logs
        mkdir -p $OUTPUT_DIR/demuxalot
        mkdir -p $OUTPUT_DIR/demuxlet
        mkdir -p $OUTPUT_DIR/ensemblex_gt
        mkdir -p $OUTPUT_DIR/souporcell
        mkdir -p $OUTPUT_DIR/vireo_gt
        rm -rf $OUTPUT_DIR/job_info/.tmp ; mkdir -p $OUTPUT_DIR/job_info/.tmp
  fi
    touch $JOB_OUTPUT_DIR/summary_report.txt
fi 

if [[ ${MODE0[@]} =~ init-noGT ]]; then
  if [[ -s $OUTPUT_DIR/job_info ]]; then
    read -p "folder\files already exist in $OUTPUT_DIR. Overwrite? y|n: "$'\n'  answer
    if [[ $answer =~ y ]]; then 
      echo "NOTE: the folder\files are Overwritten."
      rm -rf $OUTPUT_DIR/job_info; mkdir -p $OUTPUT_DIR/job_info
      rm -rf $OUTPUT_DIR/job_info/configs; mkdir -p $OUTPUT_DIR/job_info/configs
      cp -r $PIPELINE_HOME/nogt/configs $OUTPUT_DIR/job_info/
      rm -rf $OUTPUT_DIR/job_info/logs; mkdir -p $OUTPUT_DIR/job_info/logs
      rm -rf  $OUTPUT_DIR/input_files; mkdir -p $OUTPUT_DIR/input_files
      rm -rf $OUTPUT_DIR/demuxalot; mkdir -p $OUTPUT_DIR/demuxalot
      rm -rf $OUTPUT_DIR/freemuxlet; mkdir -p $OUTPUT_DIR/freemuxlet
      rm -rf $OUTPUT_DIR/souporcell; mkdir -p $OUTPUT_DIR/souporcell
      rm -rf $OUTPUT_DIR/vireo; mkdir -p $OUTPUT_DIR/vireo
      rm -rf $OUTPUT_DIR/job_info/.tmp ; mkdir -p $OUTPUT_DIR/job_info/.tmp
    else
     echo "NOTE: the pipeline is using the existing config file and parameters."
    fi
  else 
    echo "The configuration files do not exist, the pipeline will create them during execution."
    mkdir -p $OUTPUT_DIR/job_info
    mkdir -p $OUTPUT_DIR/input_files
    mkdir -p $OUTPUT_DIR/job_info/configs
    cp -r $PIPELINE_HOME/nogt/configs $OUTPUT_DIR/job_info/
    mkdir -p $OUTPUT_DIR/job_info/logs
    mkdir -p $OUTPUT_DIR/demuxalot
    mkdir -p $OUTPUT_DIR/freemuxlet
    mkdir -p $OUTPUT_DIR/souporcell
    mkdir -p $OUTPUT_DIR/vireo
    rm -rf $OUTPUT_DIR/job_info/.tmp ; mkdir -p $OUTPUT_DIR/job_info/.tmp
  fi
  touch $JOB_OUTPUT_DIR/summary_report.txt
fi

export JOB_OUTPUT_DIR=$OUTPUT_DIR/job_info
export EXPECTED_DONE_FILES=$JOB_OUTPUT_DIR/summary_report.txt
chmod 775 $EXPECTED_DONE_FILES
export OUTPUT_DIR=$OUTPUT_DIR

echo -e "NOTE: the pipeline is running  << $MODE  >>"

TEMPCONFIG=$OUTPUT_DIR/job_info/.tmp/temp_config.ini
touch $TEMPCONFIG
if [[ -f "TEMPCONFIG" ]]; then 
  rm $TEMPCONFIG
fi

echo " # IT IS A temp FILE. DO NOT EDIT THIS FILE DIRECTLY."  > $TEMPCONFIG

if [[ ${MODE0[@]} =~ init-GT ]]; then
  echo METHODGT=GT >> $TEMPCONFIG
  echo -e  "-------------------------------------------" >> $EXPECTED_DONE_FILES  
  echo -e  "--------Pipeline is set up for GT-----------------" $VERSION >> $EXPECTED_DONE_FILES
  echo "-----------------------------------------------------------"  >> $EXPECTED_DONE_FILES
  echo "The Output is under ${OUTPUT_DIR}/job_info/" >> $EXPECTED_DONE_FILES
fi 
if [[ ${MODE0[@]} =~ init-noGT ]]; then
  echo METHODGT=noGT >> $TEMPCONFIG
  echo -e  "-------------------------------------------" >> $EXPECTED_DONE_FILES  
  echo -e  "--------Pipeline is set up for noGT-----------------" $VERSION >> $EXPECTED_DONE_FILES
  echo "-----------------------------------------------------------"  >> $EXPECTED_DONE_FILES
  echo "The Output is under ${OUTPUT_DIR}/job_info/" >> $EXPECTED_DONE_FILES
fi 

echo OUTPUT_DIR=$OUTPUT_DIR >> $TEMPCONFIG
echo JOB_OUTPUT_DIR=$JOB_OUTPUT_DIR  >> $TEMPCONFIG
echo EXPECTED_DONE_FILES=$EXPECTED_DONE_FILES   >> $TEMPCONFIG
echo MODE=$MODE >> $TEMPCONFIG
echo VERSION=$VERSION >> $TEMPCONFIG

if [[ $CONTAINER_CMD ]]; then 
    echo CONTAINER_CMD=$CONTAINER_CMD >> $TEMPCONFIG
fi
if [[ $CONTAINER_MODULE ]]; then 
    echo CONTAINER_MODULE=`echo \'$CONTAINER_MODULE\'` >> $TEMPCONFIG
fi
if [[ $CONTAINER ]]; then 
    echo CONTAINER=$CONTAINER >> $TEMPCONFIG
fi

if [[ $R_MODULE ]]; then 
    echo R_MODULE=$R_MODULE >> $TEMPCONFIG
fi

if [[ $R_CMD ]]; then 
    echo R_CMD=$R_CMD >> $TEMPCONFIG
fi

R_LIB_PATH_CONT=/opt/R/4.3.1/lib/R/library
echo R_LIB_PATH_CONT=$R_LIB_PATH_CONT >> $TEMPCONFIG
PIPELINE_HOME_CONT=/opt/ensemblex.pip
echo PIPELINE_HOME_CONT=$PIPELINE_HOME_CONT >> $TEMPCONFIG

source $OUTPUT_DIR/job_info/configs/ensemblex_config.ini
source $OUTPUT_DIR/job_info/.tmp/temp_config.ini



if  [[  ${METHOD} == GT  ]] ; then
   bash ${PIPELINE_HOME}/launch/launch_gt.sh
elif  [[  ${METHOD} == noGT  ]] ; then
   bash ${PIPELINE_HOME}/launch/launch_nogt.sh
fi

STEP=sort
if [[  ${MODE0[@]}  =~  sort ]]  ; then
  echo -e "\n\n-----------------------------------------------------------" >> $EXPECTED_DONE_FILES
  echo -e  "--------Job submitted using pipeline version $VERSION--------"  >> $EXPECTED_DONE_FILES
  echo "-----------------------------------------------------------"  >> $EXPECTED_DONE_FILES
  echo "STEP $STEP: sort vcf same as bam"  >> $EXPECTED_DONE_FILES
  CONTAINER1=$PIPELINE_HOME/soft/ensemblex.sif
  $CONTAINER_CMD exec --bind  ${PIPELINE_HOME},$OUTPUT_DIR,$VCF,$BAM ${CONTAINER1}  ${PIPELINE_HOME}/tools/sort_vcf_same_as_bam.sh $BAM $VCF > $SORTOUT  2> $JOB_OUTPUT_DIR/logs/${STEP}.$(date +%FT%H.%M.%S) 
  echo "This function prepared by Gert Hulselmans" >> $EXPECTED_DONE_FILES
  echo "-------------------------------------------" >> $EXPECTED_DONE_FILES
  echo "The output is under ${OUTPUT_DIR}/souporcell" >> $EXPECTED_DONE_FILES
  echo -e " \n"
  exit 0
fi 


echo -e " \n"
exit 0

