import argparse
import sys
import os


parser = argparse.ArgumentParser()

parser.add_argument('-fl', '--folder',action='store', type=str,required=False, help='outputfolder') 
parser.add_argument('-p1', '--pdgn',action='store',required=False, help='PAR_demuxalot_genotype_names') 
parser.add_argument('-p2', '--pdps',action='store', type=int,required=False, help='PAR_demuxalot_prior_strength') 
parser.add_argument('-p3', '--pdmc',action='store', type=int,required=False, help='PAR_demuxalot_minimum_coverage') 
parser.add_argument('-p4', '--pdmac',action='store', type=int,required=False, help='PAR_demuxalot_minimum_alternative_coverage') 
parser.add_argument('-p5', '--pdnbspd',action='store', type=int, required=False, help='PAR_demuxalot_n_best_snps_per_donor')
parser.add_argument('-p6', '--pdgps',action='store',type=int,required=False, help='PAR_demuxalot_genotypes_prior_strength') 
parser.add_argument('-p7', '--pddp',action='store', type=float,required=False, help='PAR_demuxalot_doublet_prior') 
args = parser.parse_args()
folder=args.folder
PAR_demuxalot_genotype_names=args.pdgn.split(",")

PAR_demuxalot_prior_strength=args.pdps
PAR_demuxalot_minimum_coverage=args.pdmc
PAR_demuxalot_minimum_alternative_coverage=args.pdmac
PAR_demuxalot_n_best_snps_per_donor=args.pdnbspd
PAR_demuxalot_genotypes_prior_strength=args.pdgps
PAR_demuxalot_doublet_prior=args.pddp


from pathlib import Path
import pandas as pd
import pysam

from demuxalot.utils import download_file
from demuxalot import BarcodeHandler, ProbabilisticGenotypes, Demultiplexer, count_snps, detect_snps_positions
from pysam import VariantFile
import pandas as pd
import io
from demuxalot import utils

print('Part I')
handler = BarcodeHandler.from_file(os.path.join(folder,'input_files/pooled_barcodes.tsv')) ### modify
print(handler)
genotype_names = PAR_demuxalot_genotype_names
genotype_names.sort()
genotypes = ProbabilisticGenotypes(genotype_names=genotype_names)     

print('Part II')

outfolder=os.path.join(folder,'freemuxlet')
cmd_unzip=f'gunzip  {outfolder}/outs.clust1.vcf.gz'
os.system(cmd_unzip)
genotypes.add_vcf(os.path.join(folder,'freemuxlet/outs.clust1.vcf'), prior_strength=PAR_demuxalot_prior_strength) ### modify
pysam.index(os.path.join(folder,'input_files/pooled_bam.bam')) ### modify
print('Part III')
counts = count_snps(
    bamfile_location=os.path.join(folder,'input_files/pooled_bam.bam'), ### modify
    chromosome2positions=genotypes.get_chromosome2positions(),
    barcode_handler=handler,
)

utils.summarize_counted_SNPs(counts)
print('Part IV')

new_snps_filename = 'new_snps_single_file.betas'
_ = detect_snps_positions(
    bamfile_location=str(os.path.join(folder,'input_files/pooled_bam.bam')), ### modify
    genotypes=genotypes,
    barcode_handler=handler,
    minimum_coverage=PAR_demuxalot_minimum_coverage,
    minimum_alternative_coverage=PAR_demuxalot_minimum_alternative_coverage,
    result_beta_prior_filename=str(os.path.join(folder,'demuxalot/new_snps_single_file.betas')),
    n_best_snps_per_donor=PAR_demuxalot_n_best_snps_per_donor,
)

print('Part V')

genotypes_with_new_snps = genotypes.clone()
genotypes_with_new_snps.add_prior_betas(str(os.path.join(folder,'demuxalot/new_snps_single_file.betas')), prior_strength=PAR_demuxalot_genotypes_prior_strength)
counts_enriched = count_snps(
    bamfile_location=str(os.path.join(folder,'input_files/pooled_bam.bam')), ### modify
    chromosome2positions=genotypes_with_new_snps.get_chromosome2positions(),
    barcode_handler=handler,
)
print('Part VI')
learnt_enriched_genotypes, probs_learning_new_snps = Demultiplexer.learn_genotypes(
    counts_enriched,
    genotypes_with_new_snps,
    barcode_handler=handler,
    doublet_prior=PAR_demuxalot_doublet_prior,
)

print('Part V')

probs_learning_new_snps.head()
# print('save to the following')
# print(os.path.join(folder,'demuxalot/Demuxalot_result.csv'))
probs_learning_new_snps.to_csv(os.path.join(folder,'demuxalot/Demuxalot_result.csv')) ### modify


cmd_zip=f'gunzip  {outfolder}/outs.clust1.vcf'
os.system(cmd_zip)

