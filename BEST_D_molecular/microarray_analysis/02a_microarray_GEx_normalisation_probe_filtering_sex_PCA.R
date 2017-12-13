#############################
# To be run after 02 normalisation and filtering script
# Antonio J Berlanga-Taylor
# Nov 2016
# Runs PCA on gene expression data, code for BESTD project
#############################

# TO DO: most of this is project specific, clean up code...

#############################
##Set working directory and file locations and names of required inputs:

# Working directory:
#setwd("/ifs/projects/proj043/analysis.dir/gene_expression_2.dir")
# setwd('/Users/antoniob/Documents/quickstart_projects/BEST_D_molecular.p_q/results/repro_re_runs/')


#Direct output to file as well as printing to screen (plots aren't redirected though, each done separately). 
#Input is not echoed to the output file either.
output_file <- file(paste("R_session_output_PCA",".txt", sep=""), open='a')
output_file
sink(output_file, append=TRUE, split=TRUE, type = c("output", "message"))

#If the script can run from end to end, use source() with echo to execute and save all input 
#to the output file (and not truncate 150+ character lines):
#source(script_file, echo=TRUE, max.deparse.length=10000)

#Record start of session and locations:
Sys.time()
print(paste('Working directory :', getwd()))

##TO DO extract parameters:

# Read results from 01_microarrayxxx file (saved as RData object):
# Load a previous R session, data and objects:
# load('../../data/raw/R_session_saved_image_probe_filtering_full_noSNP_filter.RData')
load('../../data/raw/R_session_saved_image_read_and_QC_full.RData')
load('../../data/raw/R_session_saved_image_probe_filtering.RData')


# Filename to save current R session, data and objects at the end:
R_session_saved_image <- paste('R_session_saved_image_PCA', '.RData', sep='')

#############################


#############################
## Update packages if necessary and load them:

library(limma)
library(plyr)
library(ggplot2)
# library(gridExtra)
# library(grid)
library(cowplot)

# source('/ifs/devel/antoniob/projects/BEST-D/moveme.R')
source('../../code/BEST_D_molecular/utilities/moveme.R')
# source('/ifs/devel/antoniob/projects/BEST-D/microarray_analysis/gene_expression_functions.R')
source('../../code/BEST_D_molecular/utilities/gene_expression_functions.R')
source('../../code/BEST_D_molecular/utilities/ggtheme.R')
#############################


#############################
# Check gene expression values:
dim(normalised_filtered)
head(normalised_filtered)
dim(normalised_filtered)
class(normalised_filtered)
summary(normalised_filtered)
#############################

#############################
## Visual inspection of normalised and filtered dataset: 
# TO DO: I changed the normalised_filtered object from an Elist to a dataframe, downstream will now error.
# TO DO: move PCA analysis to a different file separate to probe filtering?
#Plot log2 intensity of regular probes normalised:
svg('regular_probes_normalised.svg')
boxplot(normalised_filtered,range= 0, ylab =expression(log[2](intensity)), 
        las = 2, xlab = "", main = "Regular probes, normalised")
dev.off()

#Plot expressed probes in a multi-dimensional scaling plot:
svg('MDS_normalised_filtered.svg')
plotMDS(normalised_filtered, pch=1)
dev.off()

svg('MDS_normalised_filtered_by_targets.svg')
plotMDS(normalised_filtered, pch=1, labels=normalised_filtered$targets)
dev.off()

# svg("MDS_normalised_filtered_by_sample.svg")
# plotMDS(normalised_filtered, pch=1, labels=normalised_filtered$targets) 
# dev.off()
###############################

###############################
# Plot PCA of normalised samples:
# Compute the PCs, first transpose the expression values and ignore the first column, then run the PCs:
pca_normalised_filtered <- prcomp(t(normalised_filtered), center=TRUE, scale=TRUE)

# Obtain values for all PCs output:
pc <- data.frame(round(pca_normalised_filtered$x, 2))
pc$sample_id <- rownames(pc) 
pc <- pc[, moveme(names(pc), 'sample_id first')]
names(pc)[1:10]
class(pc)
write.table(pc, 'principal_components_normalised_filtered.tsv', quote = FALSE,
            sep = '\t', row.names = FALSE)

# Explore dimensions and plot first 10 or so components:
dim(pc)
dim(normalised_filtered)
str(pca_normalised_filtered)
head(pc)
# pc[1:5, 1:5]
sum_pca <- summary(pca_normalised_filtered)
sum_pca$importance[, 1:10]
sum_pca
sum_pca_df <- as.data.frame(sum_pca$importance)
sum_pca_df <- t(sum_pca_df)
sum_pca_df <- as.data.frame(sum_pca_df)
# View(sum_pca_df)
sum_pca_df$percent_var <- round(100 * (sum_pca_df$`Proportion of Variance`), 3)
sum_pca_df$PC <- factor(row.names(sum_pca_df), levels = row.names(sum_pca_df),
                        labels = row.names(sum_pca_df))
head(sum_pca_df)
tail(sum_pca_df)
sum_pca_df[1:20, ]
names(sum_pca_df)
str(sum_pca_df)
summary(sum_pca_df[1:100, 'percent_var'])
###############################

###############################
# Plot proportion of variance of first x PCs:
plot_prop_vars <- ggplot(sum_pca_df[1:100, ], aes(y = percent_var, x = PC)) + 
  geom_bar(stat = 'identity') +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = 'Principal component', y = 'Proportion of variance (%)') +
  theme_Publication() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  theme(panel.grid = element_blank(), panel.border = element_blank())
ggsave('plot_PCA_normalised_filtered_prop_var.svg', plot_prop_vars)
###############################

###############################
# Plot PCA results:
svg('plot_PCA_normalised_filtered.svg')
par(mfrow=c(2,2))
# Histogram of first x PCs:
plot(pca_normalised_filtered, main='Normalised expression values')
# Scatterplot of PC1 and PC2:
biplot(pca_normalised_filtered, main='Normalised expression values')
par(mfrow=c(1,1))
dev.off()
###############################


###############################
# Will error if raw values aren't loaded.
# Plot PCA for non-normalised expression values:
# Plot PCA of normalised samples:
pca_raw <- prcomp(t(read_files_cleaned_QC$E), center=TRUE, scale=TRUE)
pc_raw <- data.frame(round(pca_raw$x, 2))
pc_raw$sample_id <- rownames(pc_raw) 
pc_raw <- pc_raw[, moveme(names(pc_raw), 'sample_id first')]
names(pc_raw)[1:10]
class(pc_raw)
dim(pc_raw)
dim(read_files_cleaned_QC$E)
str(pca_raw)
head(pc_raw)
sum_pca_raw <- summary(pca_raw)
sum_pca_raw$importance[, 1:10]

svg('plot_PCA_raw.svg')
ggplot(data = pc_raw, aes(x = PC1, y = PC2)) + geom_point(alpha = 0.7) +
  theme_Publication()
dev.off()
###############################

###############################
# TO DO: clean up
# PCA for normalised values:
# Run PCA analysis by groups of interest: TO DO: Cross files and IDs first
head(membership_file_cleaned)
tail(membership_file_cleaned)
dim(membership_file_cleaned)
head(normalised_filtered)
dim(normalised_filtered)
dim(normalised_filtered_annotated)

pc_data <- data.frame(pca_normalised_filtered$x[,1:13])
str(pc_data)
head(pc_data)

pca_by_groups <- data.frame(merge(membership_file_cleaned, pc_data, by='row.names'))
head(pca_by_groups)
dim(pca_by_groups)
dim(pc_data)
# View(pc_data)
# View(pca_by_groups)
#pc_data['120000222',]
head(arrange(pc_data, PC1), 10)
head(arrange(pca_by_groups, PC1), 10)
###############################

###############################
# Plot by treatment and time point
# Set up factors:
time_points <- factor(membership_file_cleaned$visit_type, levels = c('Randomisation', 'FinalVisit'), 
                      labels = c('Randomisation', 'Final Visit'))
treatments <- factor(membership_file_cleaned$treatment, levels = c("untreated", "treated_2000", "treated_4000"), 
                     labels = c("Untreated", "2000 IU", "4000 IU"))

svg('plot_PCA_normalised_filtered_by_groups_1a.svg')
# Panel A:
# unfilled circle = 1, black randomisation
# filled circles = 16, black final visit
p2 <- ggplot(data=pca_by_groups, aes(x=PC1, y=PC2, colour=time_points, shape=time_points)) + 
  theme_bw() + geom_point(alpha = 0.8) + 
  # scale_shape_manual(values=c(19, 7)) +
  scale_shape_manual(values=c(1, 16)) +
  scale_colour_manual(values = c("black", "black")) + #c('red', 'blue')
  theme_Publication() +
  theme(legend.position ="bottom", legend.title = element_blank()) +
  guides(col = guide_legend(title.position = 'top'))

# Panel B: 
# red, filled circle untreated = 1 
# green, filled circle 2000IU = 16
# blue, filled circle 400IU = 16
p3 <- ggplot(data=pca_by_groups, aes(x=PC1, y=PC2, colour=treatments, shape=treatments)) + 
  theme_bw() + geom_point(alpha = 0.8) + scale_shape_manual(values = c(19, 15, 17)) + #c(1, 16, 16)
  # scale_color_brewer(palette="Dark2") + # used for genotype boxplots but too pale for small symbols
  scale_colour_manual(values = c("red","darkgreen", "darkblue")) + #c("purple","darkgreen", "black")
  theme_Publication() +
  theme(legend.position ="bottom", legend.title = element_blank()) +
  guides(col = guide_legend(title.position = 'top'))
plot_grid(p2, p3, nrow = 1, labels = c('A', 'B'))
dev.off()

svg('plot_PCA_normalised_filtered_by_groups_1.svg')
p1 <- qplot(x=PC1, y=PC2, data=pca_by_groups, colour=factor(membership_file_cleaned$arm)) + 
  theme_Publication() +
  theme(legend.position = "bottom") +
  guides(col = guide_legend(title.position = 'top'))
p2 <- qplot(x=PC1, y=PC2, data=pca_by_groups, colour=factor(membership_file_cleaned$visit_type)) + 
  theme_Publication() +
  theme(legend.position="bottom") +
  guides(col = guide_legend(title.position = 'top'))
p3 <- qplot(x=PC1, y=PC2, data=pca_by_groups, colour=factor(membership_file_cleaned$treatment)) + 
  theme_Publication() +
  theme(legend.position="bottom") +
  guides(col = guide_legend(title.position = 'top'))
p4 <- qplot(x=PC2, y=PC3, data=pca_by_groups, colour=factor(membership_file_cleaned$treatment)) + 
  theme_Publication() +
  theme(legend.position="bottom") +
  guides(col = guide_legend(title.position = 'top'))
plot_grid(p1, p2, p3, p4, ncol=2, labels = c('A', 'B', 'C', 'D'))
dev.off()

svg('plot_PCA_normalised_filtered_by_groups_2.svg')
p5 <- qplot(x=PC3, y=PC4, data=pca_by_groups, colour=factor(membership_file_cleaned$treatment)) +
  theme_Publication() +
  theme(legend.position="bottom") +
  guides(col = guide_legend(title.position = 'top'))
p6 <- qplot(x=PC4, y=PC5, data=pca_by_groups, colour=factor(membership_file_cleaned$treatment)) +
  theme_Publication() +
  theme(legend.position="bottom") +
  guides(col = guide_legend(title.position = 'top'))
p7 <- qplot(x=PC5, y=PC6, data=pca_by_groups, colour=factor(membership_file_cleaned$treatment)) +
  theme_Publication() +
  theme(legend.position="bottom") +
  guides(col = guide_legend(title.position = 'top'))
p8 <- qplot(x=PC6, y=PC7, data=pca_by_groups, colour=factor(membership_file_cleaned$treatment)) +
  theme_Publication() +
  theme(legend.position="bottom") +
  guides(col = guide_legend(title.position = 'top'))
plot_grid(p5, p6, p7, p8, ncol=2, labels = c('A', 'B', 'C', 'D'))
dev.off()

svg('plot_PCA_normalised_filtered_by_groups_3.svg')
p9 <- qplot(x=PC7, y=PC8, data=pca_by_groups, colour=factor(membership_file_cleaned$treatment)) +
  theme_Publication() +
  theme(legend.position="bottom") +
  guides(col = guide_legend(title.position = 'top'))
p10 <- qplot(x=PC8, y=PC9, data=pca_by_groups, colour=factor(membership_file_cleaned$treatment)) +
  theme_Publication() +
  theme(legend.position="bottom") +
  guides(col = guide_legend(title.position = 'top'))
p11 <- qplot(x=PC9, y=PC10, data=pca_by_groups, colour=factor(membership_file_cleaned$treatment)) +
  theme_Publication() +
  theme(legend.position="bottom") +
  guides(col = guide_legend(title.position = 'top'))
p12 <- qplot(x=PC10, y=PC11, data=pca_by_groups, colour=factor(membership_file_cleaned$treatment)) +
  theme_Publication() +
  theme(legend.position = "bottom") +
  guides(col = guide_legend(title.position = 'top'))
plot_grid(p9, p10, p11, p12, ncol=2, labels = c('A', 'B', 'C', 'D'))
dev.off()
###############################

###############################
# Plot PCA by treatment group individually at baseline and 12 months:
head(pca_by_groups)
tail(pca_by_groups)
dim(pca_by_groups)
dim(pc_data)
# View(pc_data)
# View(pca_by_groups)
#pc_data['120000222',]
head(arrange(pc_data, PC1), 10)
head(arrange(pca_by_groups, PC1), 10)

# Add group membership (eg baseline and 12 months for 4000IU only):
count(pca_by_groups$arm)
count(pca_by_groups$visit_type)
pca_by_groups$arm_treatment <- ifelse(pca_by_groups$arm == 0 & 
                                        pca_by_groups$visit_type == 'Randomisation', 'baseline_4000IU',
                                      ifelse(pca_by_groups$arm == 0 & 
                                               pca_by_groups$visit_type == 'FinalVisit', '12months_4000IU',
                                             ifelse(pca_by_groups$arm == 1 & 
                                                      pca_by_groups$visit_type == 'Randomisation', 'baseline_2000IU',
                                                    ifelse(pca_by_groups$arm == 1 & 
                                                             pca_by_groups$visit_type == 'FinalVisit', '12months_2000IU',
                                                           ifelse(pca_by_groups$arm == 2 & 
                                                                    pca_by_groups$visit_type == 'Randomisation', 'baseline_placebo',
                                                                  ifelse(pca_by_groups$arm == 2 & 
                                                                           pca_by_groups$visit_type == 'FinalVisit', '12months_placebo', 'NA'
                                                                  ))))))
count(pca_by_groups$arm_treatment)
head(pca_by_groups)
# Setup labels and plotting for each allocation group:
svg('plot_PCA_normalised_filtered_by_allocation.svg')
allocation_group_placebo <- factor(pca_by_groups$arm_treatment, levels = c('baseline_placebo', '12months_placebo'),
                                   labels = c('Baseline placebo', '12 months placebo'))
# allocation_group_placebo <- na.omit(allocation_group_placebo)
summary(allocation_group_placebo)
summary(pca_by_groups$PC1)
summary(pca_by_groups$PC2)

p1 <- ggplot(data = pca_by_groups, aes(x = PC1,
                                 y = PC2,
                                 colour = allocation_group_placebo,
                                 shape = allocation_group_placebo)) + 
  geom_point(alpha = 0.8) +
  scale_shape_manual(values=c(19, 7), na.translate = FALSE) +
  scale_colour_manual(values = c("red", "blue"), na.translate = FALSE) +
  theme_Publication() +
  theme(legend.position="bottom", legend.title = element_blank())

allocation_group_2000 <- factor(pca_by_groups$arm_treatment, levels = c('baseline_2000IU', '12months_2000IU'),
                                labels = c('Baseline 2000 IU', '12 months 2000 IU'))
p2 <- ggplot(data=pca_by_groups, aes(x=PC1,
                                     y=PC2,
                                     colour=allocation_group_2000,
                                     shape=allocation_group_2000)) + 
  geom_point(alpha = 0.8) +
  scale_shape_manual(values=c(19, 7), na.translate = FALSE) +
  scale_colour_manual(values = c("red", "blue"), na.translate = FALSE) +
  theme_Publication() +
  theme(legend.position="bottom", legend.title = element_blank())

allocation_group_4000 <- factor(pca_by_groups$arm_treatment, levels = c('baseline_4000IU', '12months_4000IU'),
                                labels = c('Baseline 4000 IU', '12 months 4000 IU'))
p3 <- ggplot(data=pca_by_groups, aes(x=PC1,
                                     y=PC2,
                                     colour=allocation_group_4000,
                                     shape=allocation_group_4000)) + 
  geom_point(alpha = 0.8) +
  scale_shape_manual(values=c(19, 7), na.translate = FALSE) +
  scale_colour_manual(values = c("red", "blue"), na.translate = FALSE) +
  theme_Publication() +
  theme(legend.position="bottom", legend.title = element_blank())
# Plot together with proportion of variance:
# svg('plot_PCA_normalised_filtered_prop_var.svg', width = 12, height = 12, units = 'in', res = 300)
# p4 <- ggplot(sum_pca_df[1:100, ], aes(y = percent_var, x = PC)) + 
#   geom_bar(stat = 'identity') +
#   theme_classic() +
#   labs(x = 'Principal component', y = 'Proportion of variance (%)') +
#   theme_Publication() +
#   theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
plot_grid(plot_prop_vars, p1, p2, p3,
          # align = 'vh',
          align = 'h', axis = 'b', # aligning vertically along the left axis
          nrow = 2,
          labels = c('A', 'B', 'C', 'D'))
dev.off()
###############################

###############################
# Plot all allocations groups in one figure:
allocation_group <- factor(pca_by_groups$arm_treatment, levels = c('baseline_placebo', '12months_placebo',
                                                                   'baseline_2000IU', '12months_2000IU',
                                                                   'baseline_4000IU', '12months_4000IU'),
                           labels = c('Baseline placebo', '12 months placebo',
                                      'Baseline 2000 IU', '12 months 2000 IU',
                                      'Baseline 4000 IU', '12 months 4000 IU'))

svg('plot_PCA_normalised_filtered_by_allocation_PCs.svg')
# Just plot for legend, use cowplot to extract it:
p_legend <- ggplot(data=pca_by_groups, aes(x=PC1,
                                     y=PC2,
                                     colour=allocation_group,
                                     shape=allocation_group)) +
  geom_point(alpha = 0.8) +
  scale_colour_manual(values = c("red", "blue", "purple",
                                 "darkgreen","black", "orange"),
                      na.translate = FALSE) +
  theme_Publication() +
  theme(legend.position="bottom", legend.title = element_blank())
  # theme(legend.key.size = unit(0.5, "cm"))
legend_PCs <- cowplot::get_legend(p_legend)
legend_PCs

# Actual plots:
p1 <- ggplot(data=pca_by_groups, aes(x=PC1,
                                     y=PC2,
                                     colour=allocation_group,
                                     shape=allocation_group)) +
  geom_point(alpha = 0.8) +
  scale_colour_manual(values = c("red", "blue", "purple",
                                 "darkgreen","black", "orange"),
                      na.translate = FALSE) +
  theme_Publication() +
  theme(legend.position = 'none')
  # theme(legend.position = "bottom", legend.title = element_blank()) +
  # theme(legend.key.size = unit(0.5, "cm"))

p2 <- ggplot(data=pca_by_groups, aes(x=PC2,
                                     y=PC3,
                                     colour=allocation_group,
                                     shape=allocation_group)) + 
  geom_point(alpha = 0.8) +
  scale_colour_manual(values = c("red", "blue", "purple",
                                 "darkgreen","black", "orange"),
                      na.translate = FALSE) +
  theme_Publication() +
  theme(legend.position = 'none')

p3 <- ggplot(data=pca_by_groups, aes(x=PC3,
                                     y=PC4,
                                     colour=allocation_group,
                                     shape=allocation_group)) + 
  geom_point(alpha = 0.8) +
  scale_colour_manual(values = c("red", "blue", "purple",
                                 "darkgreen","black", "orange"),
                      na.translate = FALSE) +
  theme_Publication() +
  theme(legend.position = 'none')

p4 <- ggplot(data=pca_by_groups, aes(x=PC4,
                                     y=PC5,
                                     colour=allocation_group,
                                     shape=allocation_group)) + 
  geom_point(alpha = 0.8) +
  scale_colour_manual(values = c("red", "blue", "purple",
                                 "darkgreen","black", "orange"),
                      na.translate = FALSE) +
  theme_Publication() +
  theme(legend.position = 'none')

# Use cowplot to arrange plots in one and share a legend:
plot_1 <- plot_grid(p1, p2, p3, p4,
                    nrow = 2,
                    align = 'vh',
                    labels = c("A", "B", "C", "D")
                    )
plot_grid(plot_1, legend_PCs, ncol = 1, rel_heights = c(1, 0.1))
dev.off()
###############################

###############################
# Plot more PCs from allocations groups in one figure:
svg('plot_PCA_normalised_filtered_by_allocation_PCs_2.svg')
# Just plot for legend, use cowplot to extract it:
p_legend <- ggplot(data=pca_by_groups, aes(x=PC5, y=PC6, colour=allocation_group, shape=allocation_group)) + 
  theme_bw() + geom_point(alpha = 0.8) + 
  scale_colour_manual(values = c("red", "blue", "purple", "darkgreen","black", "orange")) +
  theme_Publication() +
  theme(legend.position="bottom", legend.title = element_blank())
  # theme(legend.key.size = unit(0.5, "cm"))
legend_PCs <- cowplot::get_legend(p_legend)
legend_PCs
# Actual plots:
p1 <- ggplot(data=pca_by_groups, aes(x=PC5, y=PC6, colour=allocation_group, shape=allocation_group)) + 
  theme_bw() + geom_point(alpha = 0.8) + 
  scale_colour_manual(values = c("red", "blue", "purple", "darkgreen","black", "orange")) +
  theme_Publication() +
  theme(legend.position = 'none')
  # theme(legend.position="bottom", legend.title = element_blank())
p2 <- ggplot(data=pca_by_groups, aes(x=PC6, y=PC7, colour=allocation_group, shape=allocation_group)) + 
  theme_bw() + geom_point(alpha = 0.8) + 
  scale_colour_manual(values = c("red", "blue", "purple", "darkgreen","black", "orange")) +
  theme_Publication() +
  theme(legend.position = 'none')
p3 <- ggplot(data=pca_by_groups, aes(x=PC7, y=PC8, colour=allocation_group, shape=allocation_group)) + 
  theme_bw() + geom_point(alpha = 0.8) + 
  scale_colour_manual(values = c("red", "blue", "purple", "darkgreen","black", "orange")) +
  theme_Publication() +
  theme(legend.position = 'none')
p4 <- ggplot(data=pca_by_groups, aes(x=PC8, y=PC9, colour=allocation_group, shape=allocation_group)) + 
  theme_bw() + geom_point(alpha = 0.8) + 
  scale_colour_manual(values = c("red", "blue", "purple", "darkgreen","black", "orange")) +
  theme_Publication() +
  theme(legend.position = 'none')
# Use cowplot to arrange plots in one and share a legend:
plot_1 <- plot_grid(p1, p2, p3, p4,
                    nrow = 2,
                    align = 'vh',
                    labels = c("A", "B", "C", "D")
                    )
plot_grid(plot_1, legend_PCs, ncol = 1, rel_heights = c(1, 0.1))
dev.off()
###############################

###############################
# Even more:
svg('plot_PCA_normalised_filtered_by_allocation_PCs_3.svg')
# Just plot for legend, use cowplot to extract it:
p_legend <- ggplot(data=pca_by_groups, aes(x=PC9, y=PC10, colour=allocation_group, shape=allocation_group)) + 
  theme_bw() + geom_point(alpha = 0.8) + 
  scale_colour_manual(values = c("red", "blue", "purple", "darkgreen","black", "orange")) +
  theme_Publication() +
  theme(legend.position="bottom", legend.title = element_blank())
  # theme(legend.key.size = unit(0.5, "cm"))
legend_PCs <- cowplot::get_legend(p_legend)
legend_PCs
# Actual plots:
p1 <- ggplot(data=pca_by_groups, aes(x=PC9, y=PC10, colour=allocation_group, shape=allocation_group)) + 
  theme_bw() + geom_point(alpha = 0.8) + 
  scale_colour_manual(values = c("red", "blue", "purple", "darkgreen","black", "orange")) +
  theme_Publication() +
  theme(legend.position = 'none')
p2 <- ggplot(data=pca_by_groups, aes(x=PC10, y=PC11, colour=allocation_group, shape=allocation_group)) + 
  theme_bw() + geom_point(alpha = 0.8) + 
  scale_colour_manual(values = c("red", "blue", "purple", "darkgreen","black", "orange")) +
  theme_Publication() +
  theme(legend.position = 'none')
p3 <- ggplot(data=pca_by_groups, aes(x=PC11, y=PC12, colour=allocation_group, shape=allocation_group)) + 
  theme_bw() + geom_point(alpha = 0.8) + 
  scale_colour_manual(values = c("red", "blue", "purple", "darkgreen","black", "orange")) +
  theme_Publication() +
  theme(legend.position = 'none')
p4 <- ggplot(data=pca_by_groups, aes(x=PC12, y=PC13, colour=allocation_group, shape=allocation_group)) + 
  theme_bw() + geom_point(alpha = 0.8) + 
  scale_colour_manual(values = c("red", "blue", "purple", "darkgreen","black", "orange")) +
  theme_Publication() +
  theme(legend.position = 'none')
# Use cowplot to arrange plots in one and share a legend:
plot_1 <- plot_grid(p1, p2, p3, p4,
                    nrow = 2,
                    align = 'vh',
                    labels = c("A", "B", "C", "D")
)
plot_grid(plot_1, legend_PCs, ncol = 1, rel_heights = c(1, 0.1))
dev.off()
###############################

###############################
# # TO DO: Next time use these functions...
# # Plot first 10 PCs:
# svg('plot_PCA_normalised_filtered_allocation_PCs_10.svg')
# grid_list <- list()
# nums <- seq(3, 12)
# for(i in nums) {
#   # plot_PCs(i)
#   p1 <- loop_ggplot_PCs(pca_by_groups, i, i+1, allocation_group)
#   # Get one of the legends:
#   legend <- get_legend(p1)
#   # Remove the legend:
#   p1 <- p1 + theme(legend.position="none")
#   grid_list[[i]] <- p1
# }
# # Arrange ggplot2 graphs with a specific width, this errors:
# # plot_grid(grid_list[3:12], ncols = 3)
# # , legend_PCs, ncol = 1, rel_heights = c(1, 0.1))
# # grid.arrange(multiplot(plotlist = grid_list[3:12], cols = 3), legend)
# dev.off()
###############################


###############################
# TO DO:
# Other plots
# Plot dendrogram
# Prepare hierarchical cluster:
# correlation <- cor(normalised_filtered, method='pearson')
# head(correlation)
# # View(correlation)
# 
# hc <- flashClust(dist(correlation))
# str(hc)
# head(hc)
# 
# # Convert to a dendogram object:
# dendro <- as.dendrogram(hc)
# 
# # Add colours by groups
# colourCodes <- c(treated_4000="red", treated_2000="green", untreated="blue")
# 
# # Assign labels of dendrogram object with new colors:
# labels_colors(dendro) <- colourCodes[membership_file_cleaned$treatment][order.dendrogram(dendro)]
# head(colourCodes[membership_file_cleaned$treatment][order.dendrogram(dendro)])
# head(order.dendrogram(dendro))
# head(membership_file_cleaned$treatment)
# 
# # Plot simple dendrogram, labels at the same level:
# plot_dendrogram_normalised_filtered <- ('plot_dendrogram_normalised_filtered.svg')
# svg(plot_dendrogram_normalised_filtered)
# par(mfrow=c(1,2), cex = 1)
# plot(dendro, hang = -1)
# # Zoom in to a sub tree:
# plot(dendro[[1]], horiz = TRUE)
# par(mfrow=c(1,1))
# dev.off()

# # Plot correlation between samples in a heatmap:
# cor_normalised_filtered <- melt(correlation)
# head(cor_normalised_filtered)
# summary(cor_normalised_filtered)
# plot_heatmap_correlations <- ('plot_heatmap_correlations_normalised_filtered.svg')
# svg(plot_heatmap_correlations)
# p1 <- ggplot(cor_normalised_filtered, aes(Var1, Var2, fill = value)) + geom_tile() + scale_fill_gradient2(low = 'blue', high = 'yellow') +
#   theme_Publication()
# grid.arrange(p1, ncol=1)
# dev.off()
###############################


###############################
# TO DO: 
# Plot expression values in a heatmap:
#plot_heatmap_correlations <- ('plot_heatmap_correlations_normalised_filtered.svg')
#svg(plot_heatmap_correlations, width = 12, height = 12, units = 'in', res = 300)
#par(mfrow=c(1,2))
#heatmap(normalised_filtered[,1:10])
#heatmap(normalised_filtered[1:100,1:100], ColSideColors=colourCodes)
# p1 <- ggplot(cor_normalised_filtered, aes(Var1, Var2, fill = value)) + geom_tile() + scale_fill_gradient2(low = 'blue', high = 'yellow')
# grid.arrange(p1, ncol=1)
#par(mfrow=c(1,1))
#dev.off()

# # Try again with heatmap.2:
# # http://sebastianraschka.com/Articles/heatmaps_in_r.html
# 
# # Create colour palette from red to green:
# my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)
# 
# # (optional) defines the color breaks manually for a "skewed" color transition:
# col_breaks = c(seq(-1,0,length=100),  # for red
#                seq(0,0.8,length=100),              # for yellow
#                seq(0.8,1,length=100))              # for green
# 
# # creates a 5 x 5 inch image
# svg('heatmap_normalised_filtered_heatmap2.svg',    # create svg for the heat map        
#     width = 5*300,        # 5 x 300 pixels
#     height = 5*300,
#     res = 300,            # 300 pixels per inch
#     pointsize = 8)        # smaller font size
# 
# heatmap.2(normalised_filtered[1:100, 1:100], 
#           cellnote = normalised_filtered[1:100, 1:100],  # same data set for cell labels
#           main = "Correlation", # heat map title
#           notecol="black",      # change font color of cell labels to black
#           density.info="none",  # turns off density plot inside color legend
#           #          trace="none",         # turns off trace lines inside the heat map
#           margins =c(12,9),     # widens margins around plot
#           col=my_palette,       # use on color palette defined earlier 
#           #          breaks=col_breaks,    # enable color transition at specified limits
#           dendrogram="column",     # only draw a column dendrogram
#           #          Colv="NA")            # turn off column clustering
# )
# dev.off()               # close the svg device


# hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
# t <- ExpressionSet(e, AnnotatedDataFrame(tab))
# rv <- rowVars(exprs(t))
# idx <- order(-rv)[1:40]
# heatmap(exprs(t)[idx, ], col = hmcol)
#############################


#############################
# The end:
# Remove objects that are not necessary to save:
#ls()
#object_sizes <- sapply(ls(), function(x) object.size(get(x)))
#as.matrix(rev(sort(object_sizes))[1:10])
#rm(list=ls(arrayQualityMetrics_raw_cleaned, arrayQualityMetrics_preprocessed))
# objects_to_save <- (c('normalised_filtered_annotated', 'normalised_filtered', 'membership_file_cleaned', 'FAILED_QC_unique'))
# save(list=objects_to_save, file=R_session_saved_image, compress='gzip')

# To save R workspace with all objects to use at a later time:
save.image(file = R_session_saved_image_full, compress = 'gzip')

# Print session information:
sessionInfo()

q()

# Next: run the script for differential gene expression.
#############################
