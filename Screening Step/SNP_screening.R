# Author: Edwina Rossi 
# Date: 14.03.22
# Title: Screening Step for T2D SNPs 

# This is a script designed to filter summary statistics from GWAS meta-analysis
# The aim of this procedure is to create a new summary-statistics file which 
# only includes SNPs of genome-wide signifiance, in the hope of increasing power
# for downstream GxE testing

library(tidyverse)
# Read in the T2D summary statistics file
sum_stat_T2D <- read.table("t2d_sum_stat.txt", header=TRUE)
head(sum_stat_T2D)

# now need to order the SNPs according to their p-value
# sorting by p_value
sorted_by_p_value_T2D <- sum_stat_T2D[order(sum_stat_T2D[, "Pvalue"]), , drop = FALSE]
# check that the above command has sorted by p value
head(sorted_by_p_value_T2D)

# filtering to only include the SNPs that reached genome-wide significance
# i.e. only including SNPs with P-value < 5x10^-8 
filter(sorted_by_p_value_T2D, Pvalue < 5*10^-8) -> filtered_t2d_sum_stat 
head(filtered_t2d_sum_stat)

# output filtered data set as new txt file 
write.table(filtered_t2d_sum_stat, file = "filtered_t2d_sum_stat.txt", quote = FALSE)
