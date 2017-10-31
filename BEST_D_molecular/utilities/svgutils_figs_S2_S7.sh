#!/bin/bash

set -e

python svgutils_BESTD.py \
                     --plotA=S1_missing_variants_lmiss.svg \
                     --plotB=S1_missing_genotypes_vs_het_rate.svg \
                     --plotC=S1_IBD-plot.svg \
                     --plotD=S1_pca-ancestry-plot.svg \
                     -O S1_genotype_QC

python svgutils_BESTD.py \
                     --plotA=S2A_plot_PCA_raw.svg \
                     --plotB=S2B_density_plot_normalised_VSN.svg \
                     --plotC=S2C_SD_vs_means_normalised_GEx.svg \
                     -O S2_GEx_QC

python svgutils_BESTD.py \
                     --plotA=S7A_pval_diff_in_diff_GEx_VSN.svg \
                     --plotB=S7B_pval_diff_in_diff_GEx_quantile.svg \
                     --plotC=S7C_pval_joint_paired_GEx_VSN.svg \
                     --plotD=S7D_pval_joint_paired_GEx_quantile.svg \
                     -O S7_GEx_pvals


