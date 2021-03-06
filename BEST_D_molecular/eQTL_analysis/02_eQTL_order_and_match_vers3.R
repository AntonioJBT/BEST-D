#############################################
# eQTL file processing script
# Antonio J Berlanga-Taylor
# 12 May 2016


# Script takes three files, orders and matches them to each other for MatrixEQTL input
# Use after 00 and 01 scripts
########################


#############################################
##Set working directory and file locations and names of required inputs:

# Working directory:
# setwd('/ifs/projects/proj043/analysis.dir/eqtl_analysis_5.dir/')
# setwd('/Users/antoniob/Desktop/BEST_D.DIR/mac_runs_to_upload/gex_FC_tests/')

#Direct output to file as well as printing to screen (plots aren't redirected though, each done separately). 
#Input is not echoed to the output file either.
options(echo = TRUE)

output_file <- file(paste("R_session_output_order_match_eQTL",Sys.Date(),".txt", sep=""), open = 'a')
output_file
sink(output_file, append=TRUE, split=TRUE, type = c("output", "message"))

#If the script can run from end to end, use source() with echo to execute and save all input 
#to the output file (and not truncate 150+ character lines):
#source(script_file, echo=TRUE, max.deparse.length=10000)

#Record start of session and locations:
Sys.time()
print(paste('Working directory :', getwd()))
getwd()

##TO DO extract parameters:

# Re-load a previous R session, data and objects:
# load('R_session_saved_image_processed_files.RData', verbose=T)

# Filename to save current R session, data and objects at the end:
# R_session_saved_image <- paste('R_session_saved_image_order_and_match_2','.RData', sep='')
R_session_saved_image <- paste('R_session_saved_image_order_and_match_eQTL','.RData', sep='')
R_session_saved_image
#############################################


########################
## Load packages:
library(data.table)

# Get script with functions needed:
source('/ifs/devel/antoniob/projects/BEST-D/moveme.R')
# source('/Users/antoniob/Documents/github.dir/cgat_projects/BEST-D/moveme.R')
source('/ifs/devel/antoniob/projects/BEST-D/BEST-D/eQTL_analysis/functions_for_MatrixeQTL.R')
# source('/Users/antoniob/Documents/github.dir/cgat_projects/BEST-D/eQTL_analysis/functions_for_MatrixeQTL.R')
########################


########################
# # TO DO for other files: Master file ordered:
# phenotype_file <- 'BEST-D_phenotype_file_final.tsv'
# phenotype_data <- fread(phenotype_file, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
# phenotype_data[1:5, 1:5, with = F]
# master_IDs <- phenotype_data[, c('pt_id', 'kit_id_randomisation', 'kit_id_finalVisit'), with = F]
# master_IDs <- master_IDs[order(kit_id_randomisation)] # Watch the lack of comma as DT not DF
# # master_IDs <- transpose(master_IDs) # Looses headers
# # master_IDs <- as.list(master_IDs)
# master_IDs <- as.data.frame(master_IDs)
# master_IDs <- as.data.frame(t(master_IDs))
# class(master_IDs)
# head(master_IDs)
# row.names(master_IDs)[2] <- 'FID'
# row.names(master_IDs)

#########################


#############################################
# Run with command line arguments:
options(echo=TRUE) # to see commands in output file. TO DO: check how it works with sink() above.
args <- commandArgs(trailingOnly = TRUE)

# TO DO: pass to configuration file

geno_file <- as.character(args[1])
# geno_file <- 'genotype_data_all_treated_baseline.tsv'
# geno_file <- 'genotype_data_all_treated_final.tsv'


expr_file <- as.character(args[2])
# expr_file <- 'GEx_baseline_4000_and_2000.tsv'
# expr_file <- 'GEx_treated_4000_and_2000.tsv'
# expr_file <- 'Airwave_CPMG_Plasma.txt'
# expr_file <- 'gex_FC_GEx_treated_4000_and_2000.tsv_over_GEx_baseline_4000_and_2000.tsv'

covar_PCs_file <- as.character(args[3])
# covar_PCs_file <- 'principal_components_normalised_filtered_PC20.tsv'
# covar_PCs_file <- ''

probe_pos_file <- as.character(args[4])
# probe_pos_file <- 'biomart_QCd_probes_genomic_locations_annot_MatrixeQTL.txt'

snp_pos_file <- as.character(args[5])
# snp_pos_file <- 'snp146Common_MatrixEQTL_snp_pos.txt'

# See qsub 02 xx file for files run.
print(args)
########################

########################
# Read geno baseline and order:
geno_data <- fread(geno_file, sep = '\t', header = TRUE, stringsAsFactors = FALSE, drop = 1)
# setkey(geno_data, )
geno_data[1:5, 1:5, with = F]
# TO DO for data.table ordering if large files:
# col_order <- order(colnames(geno_data))
# col_order
# setcolorder(geno_data, c("SNP", master_IDs_read))
geno_data <- as.data.frame(geno_data[, order(colnames(geno_data)), with = F])
class(geno_data)
# Change 'SNP' column to FID after reordering:
colnames(geno_data)[ncol(geno_data)] <- 'FID'
names(geno_data)
geno_data <- geno_data[, moveme(names(geno_data), 'FID first')]
head(geno_data)
dim(geno_data)
geno_data[1:5, 1:5]
########################

########################
# Read expr baseline and order:
# Should samples as columns and probe IDs as first column
expr_data <- fread(expr_file, sep = '\t', header = TRUE, stringsAsFactors = FALSE)#, drop = c(1))
expr_data[1:5, 1:5, with = F]
expr_data <- as.data.frame(expr_data[, order(colnames(expr_data)), with = F])
class(expr_data)
# Change 'SNP' or 'rownames' column to FID after reordering:
colnames(expr_data)[ncol(expr_data)] <- 'FID'
names(expr_data)
expr_data <- expr_data[, moveme(names(expr_data), 'FID first')]
head(expr_data)
dim(expr_data)
expr_data[1:5, 1:5]
########################

########################
# Read PCs from expression data:
covar_PCs <- fread(covar_PCs_file, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
covar_PCs <- data.table::transpose(covar_PCs) # Tranpose looses headers but these were in order already
# TO DO: transpose() gives odd error (not in namespace...), runs OK outside of Rscript but don't know where the error is.
# See:
# http://stackoverflow.com/questions/23252231/r-data-table-breaks-in-exported-functions
# http://stackoverflow.com/questions/29549690/function-on-data-table-environment-errors
# http://stackoverflow.com/questions/15223367/wrapping-data-table-using-an-evaluated-call-in-a-package
covar_PCs <- as.data.frame(covar_PCs)
names(covar_PCs) <- covar_PCs[1,]
covar_PCs <- covar_PCs[-1,]
covar_PCs[1:5, 1:5]
class(covar_PCs)
covar_PCs <- covar_PCs[, order(colnames(covar_PCs))]
covar_PCs[1:5, 1:5]
dim(covar_PCs)
########################

########################
# Match each
# Geno and expr:
length(which((colnames(expr_data) %in% colnames(geno_data))))
expr_data <- expr_data[, which(colnames(expr_data) %in% colnames(geno_data))]
geno_data <- geno_data[, which(colnames(geno_data) %in% colnames(expr_data))]
identical(colnames(expr_data), colnames(geno_data))
expr_data[1:5, 1:5]
geno_data[1:5, 1:5]

# Covar PCs to expr:
length(which((colnames(covar_PCs) %in% colnames(expr_data)[-1])))
covar_PCs <- covar_PCs[, which(colnames(covar_PCs) %in% colnames(expr_data)[-1])]

# TO DO: stop if not 'TRUE':
identical(colnames(covar_PCs), colnames(expr_data)[-1])
expr_data[1:5, 1:5]
covar_PCs[1:5, 1:5]
########################

########################
# # TO DO:
# # Match gene expression probes measured to location file:
# snp_pos_data <- fread(snp_pos_file, sep = ' ', header = TRUE, stringsAsFactors = FALSE)
# head(snp_pos_data)
# dim(snp_pos_data)
# probe_pos_data <- fread(probe_pos_file, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
# head(probe_pos_data)
# dim(probe_pos_data)
# 
# # Switch classes as these are read by data.table and cause problems later:
# snp_pos_data <- as.data.frame(snp_pos_data)
# probePos <- as.data.frame(probe_pos_data)
# 
# length(which(as.character(probe_pos_data[, 1]) %in% as.character(expr_data[, 1])))
# probe_pos_data <- probe_pos_data[which(as.character(probe_pos_data[, 1]) %in% as.character(expr_data[, 1])), ]
# dim(probe_pos_data)
# head(probe_pos_data)
# 
# length(which(as.character(snp_pos_data[, 1]) %in% as.character(geno_data[, 1])))
# snp_pos_data <- snp_pos_data[which(as.character(snp_pos_data[, 1]) %in% as.character(geno_data[, 1])), ]
# dim(snp_pos_data)
# head(snp_pos_data)
########################

########################
# Write to file each, this saves with an empty first header and row names as first column, cut when read next:
# col.names = NA makes headers match but can't be used with row.names = F
write.table(expr_data, paste(expr_file, '_matched.tsv', sep = ''), sep='\t', quote = FALSE, col.names = NA)
write.table(geno_data, paste(geno_file, '_matched.tsv', sep = ''), sep='\t', quote = FALSE, col.names = NA)

covar_PCs_file_substr <- substring(as.character(covar_PCs_file), 1, 21)
expr_file_substr <- substring(as.character(expr_file), 1, 12)
write.table(covar_PCs, paste(covar_PCs_file_substr, expr_file_substr, '_matched.tsv', sep = ''), 
            sep='\t', quote = FALSE, col.names = NA)

# Cut first column:
expr_written <- paste(expr_file, '_matched.tsv', sep = '')
cmd_cut <- sprintf('cat %s | cut -f2- > cut_%s', expr_written, expr_written)
system(as.character(cmd_cut))
geno_written <- paste(geno_file, '_matched.tsv', sep = '')
cmd_cut <- sprintf('cat %s | cut -f2- > cut_%s', geno_written, geno_written)
system(as.character(cmd_cut))
covar_written <- paste(covar_PCs_file_substr, expr_file_substr, '_matched.tsv', sep = '')
cmd_cut <- sprintf('cat %s | cut -f2- > cut_%s', covar_written, covar_written)
system(as.character(cmd_cut))

#######
# More processing:
# Cut first column:
# system('cat genotype_data_all_treate | cut -f2- > cut_genotype_data_all_treated_baseline_matched.tsv')
# # Produce test file:
# system('head -n 100 genotype_data_all_treated_baseline_matched.tsv > test_genotype_data_all_treated_baseline_matched.tsv')
# 
# # # Cut first column:
# system('cat GEx_baseline_4000_and_2000_matched.tsv | cut -f2- > cut_GEx_baseline_4000_and_2000_matched.tsv')
# # Produce test file:
# system('head -n 100 GEx_baseline_4000_and_2000_matched.tsv > test_GEx_baseline_4000_and_2000_matched.tsv')

########################

########################
# Print session info:
sessionInfo()

# To save R workspace with all objects to use at a later time:
save.image(file=R_session_saved_image, compress='gzip')

q()
# Next: run the script 04_eQTL_xxx or runMepipe xxx.
#############################################