################################################
################# SATB2 PORTAL ################
################################################
################## DATA PREPARATION ############

library(readxl)
library(readr)
library(tidyverse)


# Patient variants table
## read in from xlsx file
pat.df <- read_excel("prep_data/SATB2_Patient_variants_v1.xlsx")

## write to .txt
write_delim(pat.df, "data/SATB2_Patient_variants_v1.txt", # write to data folder
            delim = "\t")
  
  # "data/Patient_variants_SLC6A1_v5.txt", 
  #                         "\t", escape_double = FALSE, col_types = cols(Chr = col_character(), 
  #                                                                       X27 = col_character(), `Cognitive Impairment` = col_character()), 
  #                         trim_ws = TRUE) %>% 
  # mutate(Genomic_pos = as.numeric(Genomic_pos), `Age at Inclusion (years)` = as.numeric(`Age at Inclusion (years)`))

# Patient variants table with SATB2 registry IDs
## read in from xlsx file
pat.df.id <- read_excel("prep_data/SATB2_Patient_variants_v1_with_IDs.xlsx")

# ## add new column with registry IDs only: remove first 6 characters 'SATB2-'
# pat.df.id$`SATB2 ID` <- gsub("^.{0,6}", "", pat.df.id$`SATB2 ID`)

## rename registry id column and make it NUM
pat.df.id <- pat.df.id %>% 
  rename(registry_id = `SATB2 ID`)

pat.df.id %>%
  mutate(registry_id = str_sub(`SATB2 ID`, 6, -1))

## write to .txt
write_delim(pat.df.id, "data/SATB2_Patient_variants_v1.txt", # write to data folder
            delim = "\t")