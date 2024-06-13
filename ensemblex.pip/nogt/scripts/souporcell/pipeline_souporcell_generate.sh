umask 002
source $OUTPUT_DIR/job_info/configs/ensemblex_config.ini
source $OUTPUT_DIR/job_info/.tmp/temp_config.ini


TEMPFILERUN=$JOB_OUTPUT_DIR/.tmp/pipeline_souporcell.sh
cat <<EOF > $TEMPFILERUN
#!/bin/bash

umask 002
source \$OUTPUT_DIR/job_info/configs/ensemblex_config.ini
source \$OUTPUT_DIR/job_info/.tmp/temp_config.ini

if [ -d $OUTPUT_DIR/souporcell ]; then
   rm -rf  $OUTPUT_DIR/souporcell
   mkdir -p $OUTPUT_DIR/souporcell
else 
   mkdir -p $OUTPUT_DIR/souporcell
fi

cd $OUTPUT_DIR/souporcell

echo "-------------------------------------------"
echo "* step souporcell No-GT submitted at \`date +%FT%H.%M.%S\`"
echo "-------------------------------------------"
echo "* PIPELINE_HOME:           \$PIPELINE_HOME"
echo "* OUTPUT_DIR:              \$OUTPUT_DIR"
echo "-------------------------------------------"
echo "------Parameters used in this step---------"
echo "* souporcell_N:                    \$PAR_souporcell_N"
echo "* souporcell_h:                    \$PAR_souporcell_h"
echo "* souporcell_threads:              \$PAR_souporcell_threads"
echo "-------------------------------------------"
echo -e "------Output of Run------------------------\n\n"
CONTAINER1=\$PIPELINE_HOME/soft/ensemblex.sif
SOFT_SOUP=/opt/souporcell
#----------------------------------------------------------------#
# START PIPELINE                                      #
#----------------------------------------------------------------#
echo "#----------------------------------------------------------------#"

echo "Start of Renamer step"
$CONTAINER_CMD exec  --bind  $OUTPUT_DIR \${CONTAINER1} \$SOFT_SOUP/renamer.py \\
  --bam $OUTPUT_DIR/input_files/pooled_bam.bam \\
  --barcodes $OUTPUT_DIR/input_files/pooled_barcodes.tsv \\
  --out $OUTPUT_DIR/souporcell/fq.fq

echo "End of Renamer Step"


echo "#----------------------------------------------------------------#"
echo "Start of Re-align step"
$CONTAINER_CMD exec  --bind  $OUTPUT_DIR \${CONTAINER1} minimap2 $PAR_minimap2 \\
$OUTPUT_DIR/input_files/reference.fa $OUTPUT_DIR/souporcell/fq.fq > $OUTPUT_DIR/souporcell/minimap.sam

echo "End of Step Re-align "

echo "#----------------------------------------------------------------#"
echo "Start of Retag step"
$CONTAINER_CMD exec  --bind  $OUTPUT_DIR \${CONTAINER1} \$SOFT_SOUP/retag.py --sam $OUTPUT_DIR/souporcell/minimap.sam --out $OUTPUT_DIR/souporcell/minitagged.bam
echo "End of Retag step"


echo "#----------------------------------------------------------------#"
echo "Start of Sorting step"
$CONTAINER_CMD exec  --bind  $OUTPUT_DIR \${CONTAINER1} samtools sort $OUTPUT_DIR/souporcell/minitagged.bam -o $OUTPUT_DIR/souporcell/minitagged_sorted.bam
$CONTAINER_CMD exec  --bind  $OUTPUT_DIR \${CONTAINER1} samtools index $OUTPUT_DIR/souporcell/minitagged_sorted.bam

echo "End of Sorting step"

echo "#----------------------------------------------------------------#"
echo "Start of Call variants step"
$CONTAINER_CMD exec  --bind  $OUTPUT_DIR \${CONTAINER1} freebayes -f $OUTPUT_DIR/input_files/reference.fa \\
$PAR_freebayes $OUTPUT_DIR/souporcell/minitagged_sorted.bam > $OUTPUT_DIR/souporcell/Pool.vcf
echo "End of Call variants step"

echo "#----------------------------------------------------------------#"
echo "Start of Vartrix step"
$CONTAINER_CMD exec  --bind  $OUTPUT_DIR \${CONTAINER1} vartrix \\
  -v $OUTPUT_DIR/souporcell/Pool.vcf \\
  -b $OUTPUT_DIR/input_files/pooled_bam.bam \\
  -f $OUTPUT_DIR/input_files/reference.fa \\
  -c $OUTPUT_DIR/input_files/pooled_barcodes.tsv \\
  --ref-matrix $OUTPUT_DIR/souporcell/ref.mtx \\
  --out-matrix $OUTPUT_DIR/souporcell/alt.mtx \\
  --scoring-method coverage \\
EOF

if [[ $PAR_vartrix_umi ]]; then
  echo "  --umi \\" >> $TEMPFILERUN
fi

if [[ $PAR_vartrix_mapq ]]  ; then
  echo "  --mapq $PAR_vartrix_mapq \\" >> $TEMPFILERUN
fi
if [[ $PAR_vartrix_threads ]]  ; then
  echo "  --threads $PAR_vartrix_threads \\" >> $TEMPFILERUN
fi
echo " " >> $TEMPFILERUN
echo " echo \" End of Vartrix \" " >> $TEMPFILERUN
echo " " >> $TEMPFILERUN

echo "#----------------------------------------------------------------#" >> $TEMPFILERUN
echo " echo \"Start of Clustering cells by genotype step\" "  >> $TEMPFILERUN
echo " $CONTAINER_CMD exec  --bind  $OUTPUT_DIR \${CONTAINER1} \$SOFT_SOUP/souporcell/target/release/souporcell \\" >> $TEMPFILERUN
echo " -a $OUTPUT_DIR/souporcell/alt.mtx \\" >> $TEMPFILERUN
echo " -r $OUTPUT_DIR/souporcell/ref.mtx \\" >> $TEMPFILERUN
echo " -b $OUTPUT_DIR/input_files/pooled_barcodes.tsv \\" >> $TEMPFILERUN
echo " -k $PAR_souporcell_k \\" >> $TEMPFILERUN
echo " -t $PAR_souporcell_t > $OUTPUT_DIR/souporcell/clusters_tmp.tsv " >> $TEMPFILERUN
echo " " >> $TEMPFILERUN
echo " echo \"End Clustering cells by genotype step\" " >> $TEMPFILERUN

echo " " >> $TEMPFILERUN
echo "#----------------------------------------------------------------#" >> $TEMPFILERUN
echo " echo \"Start of Step Calling doublets step\" " >> $TEMPFILERUN
echo " $CONTAINER_CMD exec  --bind  $OUTPUT_DIR \${CONTAINER1} \$SOFT_SOUP/troublet/target/release/troublet \\" >> $TEMPFILERUN
echo " -a $OUTPUT_DIR/souporcell/alt.mtx \\" >> $TEMPFILERUN
echo " -r $OUTPUT_DIR/souporcell/ref.mtx \\" >> $TEMPFILERUN
echo " --clusters $OUTPUT_DIR/souporcell/clusters_tmp.tsv > $OUTPUT_DIR/souporcell/clusters.tsv " >> $TEMPFILERUN
echo " " >> $TEMPFILERUN
echo " echo \"End of Calling doublets step\" " >> $TEMPFILERUN
echo " " >> $TEMPFILERUN

echo "#----------------------------------------------------------------#" >> $TEMPFILERUN
echo " echo \" Start of Genotype and ambient RNA coinference \" " >> $TEMPFILERUN
echo "$CONTAINER_CMD exec  --bind  $OUTPUT_DIR \${CONTAINER1} \$SOFT_SOUP/consensus.py \\" >> $TEMPFILERUN
echo " -c $OUTPUT_DIR/souporcell/clusters.tsv \\" >> $TEMPFILERUN
echo " -a $OUTPUT_DIR/souporcell/alt.mtx \\" >> $TEMPFILERUN
echo " -r $OUTPUT_DIR/souporcell/ref.mtx \\" >> $TEMPFILERUN
echo " --soup_out $OUTPUT_DIR/souporcell/soup.txt \\" >> $TEMPFILERUN
echo " -v $OUTPUT_DIR/souporcell/Pool.vcf \\" >> $TEMPFILERUN
echo " --vcf_out $OUTPUT_DIR/souporcell/cluster_genotypes.vcf \\" >> $TEMPFILERUN
echo " --output_dir $OUTPUT_DIR/souporcell " >> $TEMPFILERUN
echo " " >> $TEMPFILERUN

echo "echo \"End of Genotype and ambient RNA coinference step\" " >> $TEMPFILERUN
echo " " >> $TEMPFILERUN
echo "exit 0 " >> $TEMPFILERUN
echo " " >> $TEMPFILERUN
