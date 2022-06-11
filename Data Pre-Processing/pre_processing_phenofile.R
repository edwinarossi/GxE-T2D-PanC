# Author: Edwina Rossi
# Date: 14.03.22
# Title: Pre-Processing UKB Phenotype File

# This script is written to pre-process the UKB file (phenotype file)
# The phenotype file contains a set of variables extracted from the UKB (UK Biobank)
# Not all variables will be used in downstream analyses, and for ease of use
# certain variables are renamed or re-coded 
# Filtering is also done in order to exclude certain participants 

# importing libraries
library(dplyr)
library(data.table)

# setting the working directory
setwd("~/Desktop/Masters Project/UK Biobank Sample File")

# read in the raw UK BB file
df <- read.delim("final_UK_BB_variables_file_for_processing.txt", header=TRUE, sep = " ")

#Pre-processing variables-------------------------------------------------------------------------

# function to remove columns (<<- is used in order to assign as global variable) 
remove_columns = function(col_name) {
  stopifnot(!is.null(col_name)) # missing input
  df <<- df[,!grepl(col_name, colnames(df))]
  return(df)
}

# removing unnecessary columns, i.e. irrelevant fields or further assessment center visits
# list of columns to remove
cols_to_remove = c("^X30740.1.0", #second instance of blood glucose
                   "X54.1.0",    #second instance assessment center
                   "X54.2.0",    #third instance assessment center
                   "X54.3.0",    #fourth instance assessment center
                   "^X5408",     #irrelevant field - eyes
                   "^X5419",     #irrelevant field - eyes
                   "^X5430",     #irrelevant field - eyes
                   "^X5441",     #irrelevant field - eyes
                   "^X5452",     #irrelevant field - leg pain 
                   "^X5463",     #irrelevant field - leg pain 
                   "^X5474",     #irrelevant field - leg pain 
                   "^X5485",     #irrelevant field - leg pain 
                   "^X5496",     #irrelevant field - leg pain 
                   "^X41202",    #primary diagnoses ICD10 columns
                   "^X41204",    #secondary diagnoses ICD10 columns
                   "^X40006",    #cancer type ICD10 columns
                   "^X40008",    #age at cancer diagnosis 
                   "^X20001",    #cancer code self-reported columns
                   "^X84",       #cancer year/age first occurrence (linked to field 20001)
                   "^X87",       #age at non-cancer illness  
                   "^X20002",    #non-cancer illness code self-reported columns
                   "^X20009",    #interpolated age of participant, non-cancer illness first diagnosed
                   "^X40002",    #secondary causes of death ICD10
                   "^X40001.1",  #second instance of primary cause of death ICD10
                   "^X40001.2",  #third instance of primary cause of death ICD10
                   "^X22007",    #genotype measurement plate (not used)
                   "^X30750",    #HbA1c measurements (not used)
                   "sex.y")      #other sex column from merging of data frames

# using lapply to remove the columns listed in cols_to_remove
lapply(cols_to_remove, remove_columns)

# renaming some of the single column variables
df %>% rename(
  Sex = X22001.0.0, 
  Biological.sex = X31.0.0,  
  GeneticEthnicGrouping = X22006.0.0,  
  Array = X22000.0.0,
  AssessmentCentre = X54.0.0, 
  Fasting.time = X74.0.0,
  Random.glucose = X30740.0.0,
  AgeBaseline = X21022.0.0,
  Date.E11 = Date.E11.first.reported
) -> df

# Making UK BB assessment centre a factor 
df$AssessmentCentre <- as.factor(df$AssessmentCentre)

##Redefining the coding of genotype measurement batch (array)
# Currently, the genotype measurement batch variable is coded according to which
# batch a participant's sample was analysed from 
# Need to recode the variable as a categorical variable
# which indicates which one of two arrays was used for genotyping participant
# 0 = BiLEVE, 1 = Axiom 
df$Array <- ifelse(df$Array < 0, 0, 1)

#Renaming and removing PCs---------------------------------------------------------------
# I am keeping the top 10 genetic principal components to adjust for population stratificaiton
# There are a total of 40 genetic PCs in the UK BB (Field: 22009)
# Keeping 22009.1-22009.10
# Removing 22009.11-22009.40

# Removing the PCs I am not using
PCs_to_remove = c("X22009.0.11", "X22009.0.12", "X22009.0.13", "X22009.0.14", "^X22009.0.15",
                  "X22009.0.16", "X22009.0.17", "X22009.0.18", "X22009.0.19", "X22009.0.20", 
                  "X22009.0.21", "X22009.0.22", "X22009.0.23", "X22009.0.24", "X22009.0.25", 
                  "X22009.0.26", "X22009.0.27", "X22009.0.28", "X22009.0.29", "X22009.0.30", 
                  "X22009.0.31", "X22009.0.32", "X22009.0.33", "X22009.0.34", "X22009.0.35", 
                  "X22009.0.36", "X22009.0.37", "X22009.0.38", "X22009.0.39", "X22009.0.40")

lapply(PCs_to_remove, remove_columns)

# Renaming the PCs I am using
df %>% rename(
  PC1 = X22009.0.1,
  PC2 = X22009.0.2,
  PC3 = X22009.0.3,
  PC4 = X22009.0.4,
  PC5 = X22009.0.5,
  PC6 = X22009.0.6, 
  PC7 = X22009.0.7, 
  PC8 = X22009.0.8, 
  PC9 = X22009.0.9, 
  PC10 = X22009.0.10
) -> df

#Date-format-----------------------------------------------------------------------------

# Changing the date-format of "Baseline.Date" as to match the date format in 1CD10
df$BaselineDate = as.Date(df$BaselineDate, format = "%d/%m/%Y")

#Filtering-------------------------------------------------------------------------------

# filtering to exclude participants with sex chromosome aneuploidy
# in other words, only include participants with an "NA" value in this field
filter(df, is.na(Aneuploidy)) -> df

# filtering to exclude participants with a mismatch between genetic and biological sex
# in other words, only include participants which have same value in both cols
df %>% filter(df$Sex == df$Biological.sex)

# removing biological sex and aneuploidy columns
more_cols_to_remove = c("Aneuploidy", "Biological.sex")
lapply(more_cols_to_remove, remove_columns)

# filtering to exclude participant that are of non-White British ethnicity 
# in other words, only include participants that have a value "1" in this field
filter(df, GeneticEthnicGrouping == 1) -> df

#Saving file------------------------------------------------------------------------------

# saving the pre-processed UK BB file
write.table(df, file = "~/Desktop/Masters Project/UK Biobank Sample File/pre_processed_UK_BB_file.txt")
