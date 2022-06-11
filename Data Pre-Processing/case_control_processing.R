# Author: Edwina Rossi
# Date: 11.04.22
# Title: Pre-Processing UKB Phenotype File

# This script was written to define cases and controls from the phenotype file 
# The diagnoses from the UKB are organised either in separate fields (e.g. for T2D)
# or in variables with large array structures 
# PanC diagnoses must be extracted from the array variables and matched with their date
# of diagnosis in a corresponding array variable

# All diagnoses data are coded using ICD-10 codes
# T2D = E11 (all)
# PanC = C25 (all)

# importing libraries
library(dplyr)
library(data.table)

# creating a 'notin' operator
`%notin%` <- Negate(`%in%`)

# setting the working directory
setwd("~/Desktop/Masters Project/UK Biobank Sample File")

# read in the raw UK BB file
df <- read.delim("pre_processed_UK_BB_file.txt", header=TRUE, sep = " ")

# function to remove columns (<<- is used in order to assign as global variable) 
remove_columns = function(col_name) {
  stopifnot(!is.null(col_name)) # missing input
  df <<- df[,!grepl(col_name, colnames(df))]
  return(df)
}

# function to make new variable
make_variable = function(pattern = NULL, columns = NULL) {
  stopifnot(!is.null(pattern), !is.null(columns)) # missing input
  mask = apply(df[, columns], 1, function(x) {grepl(pattern, x)}) %>% as.data.frame() %>% transpose()
  sum(mask == T)
  new_column = ifelse(rowSums(mask) > 0, 1, 0)
  return(new_column)
}

# function for making new variable where diagnosis is matched with date of diagnosis for ICD10 coded variables 
make_variable_with_date_ICD10 = function(pattern = NULL, diagnosis_columns = NULL, date_columns = NULL) {
  stopifnot(!is.null(pattern), !is.null(diagnosis_columns), !is.null(date_columns)) # missing input
  mask = apply(df[, diagnosis_columns], 1, function(x) {grepl(pattern, x)}) %>% as.data.frame() %>% transpose()
  list_of_array_numbers = apply(mask, 1, function(x) {which(x == T)})
  
  disease_array_numbers_col = matrix(NA, nrow = length(list_of_array_numbers), ncol = 1) %>% as.data.frame() #10 is max number of t2d diagnoses
  for (i in 1:length(list_of_array_numbers)) {
    
    # if the patient has no diagnosis 
    if (length(list_of_array_numbers[[i]]) == 0) {
      
      disease_array_numbers_col[i,1] = 0
      
    }
    # if the patient a has exactly 1 T2D diagnosis 
    if (length(list_of_array_numbers[[i]]) == 1) {
      
      disease_array_numbers_col[i,1] = list_of_array_numbers[[i]]
      
    }
    
    # this is the case if the patient has > 1 diagnosis
    if (length(list_of_array_numbers[[i]]) > 1) {
      
      disease_array_numbers_col[i,1] = list_of_array_numbers[[i]][1]
      
    }
    
  }
  
  diagnosis_dates = df[, grep(date_columns, colnames(df))]
  
  disease_case_dates = matrix(NA, nrow = length(list_of_array_numbers), ncol = 1) %>% as.data.frame() %>% mutate(V1 = as.Date(V1, format= "%Y.%m.%d"))
  for (i in 1:nrow(df)) {
    
    if (disease_array_numbers_col[i, 1] == 0) {
      
      disease_case_dates[i, 1] = NA
      
    } else {
      
      disease_case_dates[i, 1] = as.Date(diagnosis_dates[i, disease_array_numbers_col[i, 1]])
      
    }
  }
  
  
  return(as.vector(disease_case_dates))
  
  
} 

#Defining T2D Cases----------------------------------------------------------------------
## Code to extract T2D cases (E11) and their dates, as well as CRC diagnoses and their dates 

# IDENTIFYING T2D patients
# Identifying cases of T2D and linking them to the dates of diagnosis 
# Fields 130709 contains dates of E11 first-reported
# Thus, there is no need to match E11 diagnoses to their date of diagnosis, as 
# we already have a variable containing date of diagnosis 

# counting number of ICD0 T2D cases
no_T2D_cases = sum(!is.na(df$Date.E11))

# making new binary variable of T2D cases
df$T2D = ifelse(!is.na(df$Date.E11), 1, 0) 

# IDENTIFYING PanC patients
# Identifying cases of PanC and linking them to the dates of diagnosis 
# Fields 41270 (diagnoses - ICD10) and 41280 (date of first in-patient diagnosis) are 
# related through the array structure detailed on the UKB Showcase website
# Using make_variable_with_date_ICD10 function to link the two 

# Field 41270 = diagnoses (ICD10)
# Field 41280 = date of diagnoses 

# making a new variable for ICD10 diagnoses of PanC (C25)
df$C25.date <- make_variable_with_date_ICD10("^C25", grep("^X41270", colnames(df)), "^X41280")

# counting number of ICD10 PanC cases for C25
no_PanC_cases = sum(!is.na(df$C25.date))

# making new column containing ICD10 PanC diagnoses
df$C25.diagnosis = ifelse(!is.na(df$C25.date), 1, 0)

#Defining PanC Death Cases------------------------------------------------------------
# Defining more cases from patients who have a PanC listed as their cause of death
# Grab all PanC (C25) death cause from death cause fields

# Field 40001 = Underlying (primary) cause of death: ICD10

# new variable for patients with primary death C25
df$Death.primary.C25 = ifelse(grepl("^C25", df$X40001.0.0) == T, 1, 0)

# new variable which only counts the PanC deaths
df$C25.death.only = ifelse(df$C25.diagnosis == 0 & df$Death.primary.C25 == 1, 1, 0)
C25_death_only = sum(df$C25.death.only == 1)

# final variable combining all reports of PanC (death or diagnosis)
df$C25 = ifelse(df$C25.death.only == 1 | df$C25.diagnosis == 1, 1, 0)
no_C25_cases_all = sum(df$C25 == 1)


#Defining RG Cases----------------------------------------------------------------------
# First define people with an RG measurement available
# Excluding participants with a T2D diagnosis diagnosed prior to date of assessment centre

# counting number of participants with RG measurements available
RG_measurements = sum(!is.na(df$Random.glucose))

# defining number of people with RG measurements and no T2D diagnosis
df$RG.control = ifelse((!is.na(df$Random.glucose) & df$T2D == 0), df$Random.glucose, NA)

# counting number of RG measurements without T2D diagnosis
no_RG_measurements_no_T2D = sum(!is.na(df$RG.control), na.rm = T)


# Now defining PanC Incident with RG measurements 
# Including only participants where RG measurement > PanC diagnosis 
# Compare values in PanC date of diagnosis column with values in column BaselineDate
# If BaselineDate > PanC date of diagnosis, include participant as PanC incident
# Logic: if "value" of BaselineDate < date of PanC of diagnosis, include PanC incident

# adding column which represents all PanC patients with an RG measurement and without a T2D diagnosis 
# (where RG measurement was taken before the PanC diagnosis)
# Make new column containing PANC Incident's RG value
df$PANC.incident.RG.measurements = ifelse((df$T2D == 0) & 
                                                (df$BaselineDate < df$C25.date), 
                                                 df$Random.glucose, NA)

# Binary variable for PANC incident RG 
df$PANC.incident.RG = ifelse(!is.na(df$PANC.incident.RG.measurements), 1, 0)

# counting number of PANC Incident RG cases
PANC_incident_RG = sum(df$PANC.incident.RG == 1, na.rm = T)


#Final Processing and Export--------------------------------------------------------------------------------------------------------
# removing now redundant columns
old_cols_to_remove = c("^X40001",                       #primary death ICD10 columns
                       "^X40002",                       #secondary death ICD10 columns
                       "^X41270",                       #all diagnoses ICD10 columns
                       "^X41280",                       #all diagnoses date ICD10 columns
                       "^Death.primary.C25",            #redundant primary death PanC columns
                       "C25.death.only",                #redundant PanC death variable
                       "C25.diagnosis")                 #redundant C25 diagnosis column
  
lapply(old_cols_to_remove, remove_columns)

# saving 'processed' UK BB file to be further used for summary statistics and EDA
write.table(df, file = "~/Desktop/Masters Project/UK Biobank Sample File/Pheno_file_RG_T2D_PANC", row.names = F)
