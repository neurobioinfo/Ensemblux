umask 002
source $OUTPUT_DIR/job_info/configs/ensemblex_config.ini
source $OUTPUT_DIR/job_info/.tmp/temp_config.ini


TEMPFILERUN=$OUTPUT_DIR/job_info/.tmp/pipeline_vireo.sh

cat <<EOF > $TEMPFILERUN
  #!/bin/bash

  umask 002
  source \$OUTPUT_DIR/job_info/configs/ensemblex_config.ini
  source \$OUTPUT_DIR/job_info/.tmp/temp_config.ini

  echo "-------------------------------------------"
  echo "* step vireo submitted at \`date +%FT%H.%M.%S\`"
  echo "-------------------------------------------"
  echo "* PIPELINE_HOME:           $PIPELINE_HOME"
  echo "* OUTPUT_DIR:              $OUTPUT_DIR"
  echo "-------------------------------------------"
  echo "------Parameters used in this step---------"
  echo "* N:                       $PAR_vireo_N"
  echo "* type:                    $PAR_vireo_type"
  echo "* h:                       $PAR_vireo_h"
  echo "* processes:               $PAR_vireo_processes"
  echo "* minMAF:                  $PAR_vireo_minMAF"
  echo "* minCOUNT:                $PAR_vireo_minCOUNT"
  echo "* forcelearnGT:            $PAR_vireo_forcelearnGT"
  echo "-------------------------------------------"
  echo -e "------Output of Run------------------------\n\n"
  CONTAINER1=$PIPELINE_HOME/soft/ensemblex.sif
  $CONTAINER_CMD exec  --bind  $OUTPUT_DIR \${CONTAINER1} cellSNP -s $OUTPUT_DIR/input_files/pooled_bam.bam \\
      -b $OUTPUT_DIR/input_files/pooled_barcodes.tsv \\
      -O $OUTPUT_DIR/vireo_gt \\
      -R $OUTPUT_DIR/input_files/reference.vcf \\
      -p $PAR_vireo_processes \\
      --minMAF $PAR_vireo_minMAF \\
      --minCOUNT $PAR_vireo_minCOUNT
EOF

  echo "" >> $TEMPFILERUN
  echo "" >> $TEMPFILERUN
  echo "$CONTAINER_CMD exec  --bind  $OUTPUT_DIR \${CONTAINER1} vireo -c $OUTPUT_DIR/vireo_gt \\" >> $TEMPFILERUN
  echo "        -d $OUTPUT_DIR/input_files/pooled_samples.vcf \\" >> $TEMPFILERUN
  echo "        --forceLearnGT \\" >> $TEMPFILERUN
  echo "        -t GT \\" >> $TEMPFILERUN
  echo "        -o $OUTPUT_DIR/vireo_gt \\" >> $TEMPFILERUN
  if [[ $PAR_vireo_N ]]  ;then
     echo "        -N $PAR_vireo_N \\" >> $TEMPFILERUN
  fi
  if [[ $PAR_vireo_h ]]  ;then
     echo "        -h $PAR_vireo_h \\" >> $TEMPFILERUN
  fi
  echo "" >> $TEMPFILERUN
  echo "" >> $TEMPFILERUN
  