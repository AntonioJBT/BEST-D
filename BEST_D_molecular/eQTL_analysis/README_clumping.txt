Obtained Hapmap data from:
https://neurogenetics.qimrberghofer.edu.au/iSECA/LD_clumping_tutorial.html

as well as clumping example.

See also:

https://neurogenetics.qimrberghofer.edu.au/iSECA/LD_clumping_tutorial.html
LD clumping tutorial
https://github.com/molgenis/systemsgenetics/wiki/eQTL-mapping-analysis-cookbook-for-RNA-seq-data
eQTL mapping analysis cookbook for RNA seq data · molgenis/systemsgenetics Wiki
https://privefl.github.io/bigsnpr/articles/pruning-vs-clumping.html
Why clumping should be preferred over pruning • bigsnpr
https://privefl.github.io/bigsnpr/index.html
Analysis of Massive SNP Arrays • bigsnpr
http://zhengxwen.github.io/SNPRelate/release/help/snpgdsLDpruning.html
R: Linkage Disequilibrium (LD) based SNP pruning
https://www.cog-genomics.org/plink/1.9/ld
Linkage disequilibrium - PLINK 1.9
https://www.cog-genomics.org/plink/1.9/postproc#clump
Report postprocessing - PLINK 1.9

Plink clumping commands are to obtain a simple number of LD independent SNPs
from eQTL results without reference to MAF< p-value, effect size, etc.

Used only FDR under 5%:
awk '{if ($6 < 0.05) print $0}' 2000+4000-12months-1.eQTL_trans > 2000+4000-12months-1.eQTL_trans_FDR5

then inserted header with e.g.:
gsed -i '1i\SNP\tgene\tbeta\tt-stat\tp-value\tFDR' 2000+4000-12months-1.eQTL_trans_FDR5

Ran plink clump as bash script with:
bash plink_clump.sh hapmap3_r2_b36_fwd_CEU_qc_poly 2000+4000-baseline-1.eQTL_cis_FDR5  p-value  SNP  0.0001  1.0  0.1  10  2000+4000-baseline-1.eQTL_cis_FDR5

