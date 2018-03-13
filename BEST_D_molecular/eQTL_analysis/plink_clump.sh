#!/usr/bin/env bash

###########################
# Runs basic LD based, p-value index SNP clumping
# You need plink 1.90 installed to run as plink
# Takes a list of SNPs with p-values from e.g. an eQTL study
# and returns plnk's clump results
# run as e.g.
# bash plink_clump.sh bfile=hapmap_ref_file \
#                     clump_file=my_file \   
#                     clump_field=p-value \
#                     clump_snp_field=SNP \
#                     clump_p1=0.0001 \ # miminum p-value for index SNP
#                     clump_p2=1.0 \ # minimum p-value for SNPs being added to index clump
#                     clump_r2=0.1 \
#                     clump_kb=10000 \
#                     out_file=my_results
###########################

###########################
# Some references to check:
# https://www.cog-genomics.org/plink/1.9/postproc#clump
# Report postprocessing - PLINK 1.9

# For bash:
# https://kvz.io/blog/2013/11/21/bash-best-practices/
# http://jvns.ca/blog/2017/03/26/bash-quirks/
# Bash traps:
# http://aplawrence.com/Basics/trapping_errors.html
# https://stelfox.net/blog/2013/11/fail-fast-in-bash-scripts/
###########################

###########################
# Set bash script options

# exit when a command fails
set -o errexit

# exit if any pipe commands fail
set -o pipefail

# exit when your script tries to use undeclared variables
set -o nounset

# trace what gets executed
set -o xtrace

set -o errtrace
###########################

###########################
# Variables to substitute:
bfile=${1:-hapmap3_r2_b36_fwd_CEU_qc_poly}
#${bfile:-hapmap3_r2_b36_fwd_CEU_qc_poly}
#hapmap3_r2_b36_fwd_CEU_qc_poly
clump_file=${2:-my_file}
#${clump_file:-my_file}
#my file with SNPs and p-values
clump_field=${3:=P}
#${clump_field:=P}
#p-value
clump_snp_field=${4:=SNP}
#${clump_snp_field:=SNP}
#SNP
clump_p1=${5:=0.0001}
#${clump_p1:=0.0001}
#0.0001
clump_p2=${6:=0.01}
#${clump_p2:=0.01}
#1.0
clump_r2=${7:=0.5}
#${clump_r2:=0.5}
#0.1
clump_kb=${8:=250}
#${clump_kb:=250}
#10000
out_file=${9:=my_results}
#${out_file:=my_results}
#prefix to save output
###########################

###########################
# Run command:
plink \
--bfile ${bfile} \
--clump ${clump_file} \
--clump-field ${clump_field}  \
--clump-snp-field ${clump_snp_field}  \
--clump-p1 ${clump_p1} \
--clump-p2 ${clump_p2} \
--clump-r2 ${clump_r2} \
--clump-kb ${clump_kb} \
--out ${out_file} \
--clump-verbose \
--clump-best

echo 'Done clumpling SNPs, the command included --clump-best and --clump-verbose'
###########################
