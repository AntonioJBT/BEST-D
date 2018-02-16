#############################
# To be run after 02 normalisation of array data, pheno file processing
# Antonio J Berlanga-Taylor
# 12 Feb 2018
# BEST-D project differential expression additional analyses - e.g. gender interaction
#############################


#############################
##Set working directory and file locations and names of required inputs:

# Working directory:
# setwd('/Users/antoniob/Documents/quickstart_projects/BEST_D_molecular.p_q/results/reviewers/')

#Direct output to file as well as printing to screen (plots aren't redirected though, each done separately). 
#Input is not echoed to the output file either.
output_file <- file(paste("R_session_output_diff_expression_gender",".txt", sep=""), open='a')
output_file
sink(output_file, append=TRUE, split=TRUE, type = c("output", "message"))

#If the script can run from end to end, use source() with echo to execute and save all input 
#to the output file (and not truncate 150+ character lines):
#source(script_file, echo=TRUE, max.deparse.length=10000)

#Record start of session and locations:
Sys.time()
print(paste('Working directory :', getwd()))

##TO DO extract parameters. This file is project specific. Check ways of making count comparisons.

# Load results from 02_microarrayxxx file, saved as RData object:
# Re-load a previous R session, data and objects:
#load('R_session_saved_image_probe_filtering.RData', verbose=T)
load('../../data/raw/R_session_saved_image_pheno_file_check.RData')
#load('R_session_saved_image_diff_expression_3.RData', verbose=T)
# load('../../data/raw/R_session_saved_image_diff_expression.RData', verbose=T)

# If plotting quantile normalisation data (neqc), load:
# load('R_session_saved_image_normalisation_NEQC.RData', verbose=T)
# where the object 'normalised' is run_neqc with probe and pheno file processing up to this file.
# This is for the pval_plots of SF7

# For test without SNP filter:
# load('R_session_saved_image_pheno_file_check_noSNP_filter.RData', verbose=T)

# Filename to save current R session, data and objects at the end:
R_session_saved_image <- paste('R_session_saved_image_diff_expression', '.RData', sep='')
# For test without SNP filter:
# R_session_saved_image <- paste('R_session_saved_image_diff_expression_noSNP_filter', '.RData', sep='')

#############################


#############################
## Update packages if necessary and load them:
#install.packages('ellipse')
#install.packages("statmod")

library(limma)
library(plyr)
# library(ggplot2)
# library(gridExtra)
# library(illuminaHumanv4.db)
# library(ellipse)
# library(Hmisc)
# library(splines)
# library(statmod)
# library(cowplot)

# Get functions from other scripts (eg to add annotations to topTable results):
# source('/ifs/devel/antoniob/projects/BEST-D/gene_expression_functions.R')
# source('/ifs/devel/antoniob/projects/BEST-D/moveme.R')
source('../../code/BEST_D_molecular/utilities/moveme.R')
source('../../code/BEST_D_molecular/utilities/gene_expression_functions.R')
source('../../code/BEST_D_molecular/utilities/ggtheme.R')
#############################

#############################
# Read in the file with the sample membership information (group membership, replicate number, 
# treatment status, etc.) to be able to create a design and contrast (done above for PCA plots).

# # For test without SNP filter, drop 'Status' column:
# row.names(membership_file_cleaned)
# colnames(normalised_filtered)
# normalised_filtered$Status <- NULL

# Sanity check:
# TO DO: Raise error and stop if false:
identical(row.names(membership_file_cleaned), colnames(normalised_filtered))
length(which(row.names(membership_file_cleaned) %in% colnames(normalised_filtered)))

head(membership_file_cleaned)
length(membership_file_cleaned$pt_id)
tail(membership_file_cleaned)
str(membership_file_cleaned)
#############################

#############################
# Read in phenotype data to add gender and other variables of possible interest
# Load full phenotype data for covariates adjustment:
phenotype_data <- read.csv('../../data/raw/BEST-D_phenotype_file_final.tsv', sep = '\t', 
                           header = TRUE, na.string = c('-99999', "", " ", "NA"))
dim(phenotype_data)
length(which(complete.cases(phenotype_data)))
#View(phenotype_data)
head(phenotype_data)
tail(phenotype_data)
summary(phenotype_data)
str(phenotype_data)
class(phenotype_data)
names(phenotype_data)
summary(phenotype_data[, c('kit_id_randomisation', 'vitd0', 'vitd12')])
count(phenotype_data$male)
#############################

#############################
# Get covariates and IDs from phenotype data:
head(membership_file_cleaned)
head(phenotype_data)
names(phenotype_data)
covars <- c('incident_fracture',
            'incident_resp_infection',
            'diabetes',
            'heart_disease',
            'copd_asthma',
            'basemed_vitamind',
            'currsmoker',
            'bmi0',
            'calendar_age_ra',
            'season_randomisation_2',
            'male')
IDs <- c('pt_id',
         'kit_id_randomisation',
         'kit_id_finalVisit'
)

# Merge with membership file for processing with GEx data:
membership_file_cleaned <- merge(membership_file_cleaned, phenotype_data[, c(IDs, covars)])
head(membership_file_cleaned)
# View(membership_file_cleaned)
dim(membership_file_cleaned)

# Sanity
# Number of unique individuals:
head(membership_file_cleaned$pt_id)
# View(membership_file_cleaned)
length(unique(membership_file_cleaned$pt_id)) # unique removes duplicated elements, total individuals
length(which(!duplicated(membership_file_cleaned$pt_id))) # Number of individuals with at least 1 sample
length(which(duplicated(membership_file_cleaned$pt_id))) # Number of individuals with 2 samples
length(which(grepl('100173', membership_file_cleaned[, 'pt_id'])))
# Sanity:
duplicated_df <- as.data.frame(count(membership_file_cleaned$pt_id))
head(duplicated_df)
length(which(duplicated_df$freq == 2))
length(which(duplicated_df$freq == 1))
count(duplicated_df$freq)
dim(duplicated_df)

# Counts per group:
dim(membership_file_cleaned)
head(membership_file_cleaned)
count(membership_file_cleaned$group_membership)
count(membership_file_cleaned$arm)
count(membership_file_cleaned$visit_type)
count(phenotype_data$male) # male is 1
count(membership_file_cleaned$male)
head(phenotype_data[, c('pt_id', 'male')])
head(membership_file_cleaned[, c('pt_id', 'male')])
head(unique(membership_file_cleaned[, c('pt_id', 'male')]))

tail(phenotype_data[, c('pt_id', 'male')])
tail(unique(membership_file_cleaned[, c('pt_id', 'male')]))

dim(phenotype_data[, c('pt_id', 'male')])
dim(unique(membership_file_cleaned[, c('pt_id', 'male')]))

dim(membership_file_cleaned)
#############################

#############################
# Check gene expression data:
dim(normalised_filtered)
normalised_filtered[1:5, 1:5]
length(which(complete.cases(normalised_filtered)))
#############################


#############################
# Subset by gender
dim(membership_file_cleaned)
head(membership_file_cleaned)
count(membership_file_cleaned$group_membership)
count(membership_file_cleaned$male)
summary(membership_file_cleaned$male)
which(is.na(membership_file_cleaned$male) == TRUE) # should be zero
# Final visit with array data counts:
gender_code <- 1 # 0 = female (should be 144)
count_gender <- membership_file_cleaned[which(membership_file_cleaned$male == gender_code &
                                !is.na(membership_file_cleaned$kit_id_finalVisit)), c('male',
                                                                                      'kit_id_finalVisit')
                                ]
# Remove duplicates as membership file has pt_id twice:
length(unique(count_gender$kit_id_finalVisit)) # female should be ~144, male ~

# Get women only:
summary(phenotype_data$male == gender_code)
membership_file_cleaned_women <- membership_file_cleaned[which(membership_file_cleaned$male == gender_code), ]
dim(membership_file_cleaned_women)
head(membership_file_cleaned_women)
names(membership_file_cleaned_women)
head(membership_file_cleaned_women$kit_id_randomisation)
head(membership_file_cleaned_women$kit_id_finalVisit)
is.na(membership_file_cleaned_women$kit_id_randomisation)
is.na(membership_file_cleaned_women$kit_id_finalVisit)
count(membership_file_cleaned_women$group_membership)
count(membership_file_cleaned_women$male)

# Subset array file:
kit_IDs <- cbind(membership_file_cleaned_women$kit_id_randomisation,
                 membership_file_cleaned_women$kit_id_finalVisit)
kit_IDs <- as.vector(as.matrix(membership_file_cleaned_women[,c('kit_id_randomisation',
                                                     'kit_id_finalVisit')]))
kit_IDs
dim(normalised_filtered)
colnames(normalised_filtered) # kit IDs at baseline and final visit
normalised_filtered_women <- normalised_filtered[, which(colnames(normalised_filtered) %in% kit_IDs)]
dim(normalised_filtered_women)
head(normalised_filtered_women)
#############################


#############################
# Compare joint 2000+4000 vs baseline and placebo with pairing
# Pairing and treatment:
pairing_joint <- factor(membership_file_cleaned_women$pt_id)
head(pairing_joint)
length(pairing_joint)

# Get treatment groups:
treatment_joint <- factor(membership_file_cleaned_women$all_treated, levels = c('untreated', 'treated_2000+4000', 'treated_placebo'))
dim(membership_file_cleaned_women)
length(treatment_joint)
head(treatment_joint)
summary(treatment_joint)

#Define design:
design_all_treated_pairs <- model.matrix(~pairing_joint+treatment_joint)
head(design_all_treated_pairs)[1:5, 1:5]
tail(design_all_treated_pairs)[, -1]
dim(design_all_treated_pairs)
colnames(design_all_treated_pairs)

#Run linear model and obtain differentially expressed genes based on all pairs:
fit_all_treated <- lmFit(normalised_filtered_women, design_all_treated_pairs)
fit_all_treated_2 <- eBayes(fit_all_treated)
dim(fit_all_treated_2)
str(fit_all_treated_2)
summary(fit_all_treated_2)
names(fit_all_treated_2)
head(fit_all_treated_2$coefficients)#[1:5, 1:5]
head(fit_all_treated_2$lods)#[1:5, 1:5]
head(fit_all_treated_2$cov.coefficients)#[1:5, 1:5]
fit_all_treated_2$coefficients
colnames(fit_all_treated_2)

# Get results:
topTable(fit_all_treated_2, adjust='BH')
topTable(fit_all_treated_2, coef="treatment_jointtreated_2000+4000", adjust='BH')
topTable(fit_all_treated_2, coef="treatment_jointtreated_placebo", adjust='BH')

topTable_pairing_joint_treated <- topTable(fit_all_treated_2, coef="treatment_jointtreated_2000+4000", adjust='BH', n=Inf)
topTable_pairing_joint_placebo <- topTable(fit_all_treated_2, coef="treatment_jointtreated_placebo", adjust='BH', n=Inf)
class(topTable_pairing_joint_placebo)
head(topTable_pairing_joint_treated)
dim(topTable_pairing_joint_treated)
head(topTable_pairing_joint_placebo)
dim(topTable_pairing_joint_placebo)

# Basic counts
count(topTable_pairing_joint_treated$adj.P.Val < 0.05)
count(topTable_pairing_joint_placebo$adj.P.Val < 0.05)

count(topTable_pairing_joint_treated$adj.P.Val < 0.10)
count(topTable_pairing_joint_placebo$adj.P.Val < 0.10)

count(topTable_pairing_joint_treated$adj.P.Val < 0.05 & 2^(topTable_pairing_joint_treated$logFC) > 1.1)
count(topTable_pairing_joint_treated$adj.P.Val < 0.05 & 2^(topTable_pairing_joint_treated$logFC) < 0.9)

# Interpretation: There are no significant differences for paired tests (before vs after) when only looking at women
# when joining treatment groups (placebo and 2000 + 4000).

# Write results to disk:
write.table(x=topTable_pairing_joint_treated, sep='\t', quote = FALSE, col.names = NA, row.names = TRUE,
            file='full_topTable_pairing_all_treated_women.txt')

write.table(x=topTable_pairing_joint_placebo, sep='\t', quote = FALSE, col.names = NA, row.names = TRUE,
            file='full_topTable_pairing_all_placebo_women.txt')
#############################



#############################
# Compare individual arms vs baseline and placebo
# Pairing and treatment:
# pairing_arms_women <- factor(membership_file_cleaned_women$pt_id)
# head(pairing_arms_women)
# length(pairing_arms_women)

# Get treatment groups, exclude placebo as these have no delta with > median:
count(membership_file_cleaned_women$group_membership)
treatment_arms_women <- factor(membership_file_cleaned_women$group_membership, 
                                  levels = c('baseline_2000',
                                             'baseline_4000',
                                             'final_2000',
                                             'final_4000',
                                             'baseline_placebo',
                                             'final_placebo'))
dim(membership_file_cleaned_women)
length(treatment_arms_women)
head(treatment_arms_women)
summary(treatment_arms_women)

#Define design:
design_arms_pairs_women <- model.matrix(~0+treatment_arms_women)
head(design_arms_pairs_women)
tail(design_arms_pairs_women)
dim(design_arms_pairs_women)
colnames(design_arms_pairs_women)
colSums(design_arms_pairs_women)

#Run linear model and obtain differentially expressed genes based on all pairs:
fit_arms_women <- lmFit(normalised_filtered_women, design_arms_pairs_women)

# Set contrasts:
fit_arms_women
cont.matrix <- makeContrasts(f4000vsb4000=treatment_arms_womenfinal_4000-treatment_arms_womenbaseline_4000, 
                             f2000vsb2000=treatment_arms_womenfinal_2000-treatment_arms_womenbaseline_2000,
                             fplacebovsbplacebo=treatment_arms_womenfinal_placebo-treatment_arms_womenbaseline_placebo,
                             levels=design_arms_pairs_women)
cont.matrix

# Obtain differentially expressed genes based on contrasted factors:
fit_arms_2_women <- contrasts.fit(fit_arms_women, cont.matrix)
fit_arms_2_women <- eBayes(fit_arms_2_women)
names(fit_arms_2_women)
colnames(fit_arms_2_women)

# Get results:
topTable(fit_arms_2_women, adjust='BH')
topTable(fit_arms_2_women, coef="f4000vsb4000", adjust='BH')
topTable(fit_arms_2_women, coef="f2000vsb2000", adjust='BH')
topTable(fit_arms_2_women, coef="fplacebovsbplacebo", adjust='BH')

# Interpretation: There are no significant differences for tests (before vs after) when only looking at
# women within each treatment arm.

# Write results to disk:
# write_topTable("treatment_arms_womentreated_2000+4000", fit_arms_2_women)
#############################

#############################
#The end:

# To save R workspace with all objects to use at a later time:
save.image(file=R_session_saved_image, compress='gzip')
sessionInfo()

q()
# Next: run script for higher level analyses.
#############################