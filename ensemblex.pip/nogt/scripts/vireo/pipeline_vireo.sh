umask 002
source $OUTPUT_DIR/job_info/configs/ensemblex_config.ini
source $OUTPUT_DIR/job_info/.tmp/temp_config.ini

TEMPFILERUN=$OUTPUT_DIR/job_info/.tmp/pipeline_vireo.sh

cat <<EOF > $TEMPFILERUN
  #!/bin/bash

  umask 002
  source \$OUTPUT_DIR/job_info/configs/ensemblex_config.ini
  source \$OUTPUT_DIR/job_info/.tmp/temp_config.ini
  if [ -d $OUTPUT_DIR/ensemblex ]; then
      rm -rf  $OUTPUT_DIR/vireo
      mkdir -p $OUTPUT_DIR/vireo
  else 
      mkdir -p $OUTPUT_DIR/vireo
  fi
  echo "-------------------------------------------"
  echo "* step vireo No-GT submitted at \`date +%FT%H.%M.%S\`"
  echo "-------------------------------------------"
  echo "* PIPELINE_HOME:           $PIPELINE_HOME"
  echo "* OUTPUT_DIR:              $OUTPUT_DIR"
  echo "-------------------------------------------"
  echo "------Parameters used in this step---------"
  echo "* PAR_vireo_N:                       $PAR_vireo_N"
  echo "* PAR_vireo_h:                       $PAR_vireo_h"
  echo "* PAR_vireo_processes:               $PAR_vireo_processes"
  echo "* PAR_vireo_minMAF:                  $PAR_vireo_minMAF"
  echo "* PAR_vireo_minCOUNT:                $PAR_vireo_minCOUNT"
  echo "* PAR_vireo_forcelearnGT:            $PAR_vireo_forcelearnGT"
  echo "-------------------------------------------"
  echo -e "------Output of Run------------------------\n\n"
  CONTAINER1=$PIPELINE_HOME/soft/ensemblex.sif
  echo "Start of cellSNP step"
  $CONTAINER_CMD exec  --bind  $OUTPUT_DIR \${CONTAINER1} cellSNP -s $OUTPUT_DIR/input_files/pooled_bam.bam \\
      -b $OUTPUT_DIR/input_files/pooled_barcodes.tsv \\
      -O $OUTPUT_DIR/vireo \\
      -R $OUTPUT_DIR/input_files/reference.vcf \\
      -p $PAR_vireo_processes \\
      --minMAF $PAR_vireo_minMAF \\
      --minCOUNT $PAR_vireo_minCOUNT
EOF
  echo "" >> $TEMPFILERUN
  echo "End of cellSNP step" >> $TEMPFILERUN
  echo "" >> $TEMPFILERUN
  echo "Start of vireo step" >> $TEMPFILERUN
  echo "$CONTAINER_CMD exec  --bind  $OUTPUT_DIR \${CONTAINER1} vireo -c $OUTPUT_DIR/vireo \\" >> $TEMPFILERUN
  echo "        --forceLearnGT \\" >>  $TEMPFILERUN
  echo "        -t GT \\" >> $TEMPFILERUN
  echo "        -o $OUTPUT_DIR/vireo \\" >> $TEMPFILERUN
  if [[ $PAR_vireo_N ]]  ;then
     echo "        -N $PAR_vireo_N \\" >> $TEMPFILERUN
  fi
  if [[ $PAR_vireo_h ]]  ;then
     echo "        -h $PAR_vireo_h \\" >> $TEMPFILERUN
  fi
  echo "" >> $TEMPFILERUN
  echo "End of vireo step" >> $TEMPFILERUN
  echo "" >> $TEMPFILERUN
  