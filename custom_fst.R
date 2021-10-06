# custom_pairwise_fst.R
#
# Author: Mouhamadou Fadel DIOP
# Date: 2021-10-06
#
# Purpose:
# Produce per locus Fst between Gambia and Cambodia.
#
# ------------------------------------------------------------------

# NOTE - uncomment these lines to install packages as needed
# install.packages("data.table")
# install.packages("tidyverse")

# Load packages
library(data.table)
library(tidyverse)

source('craig_functions.r')

# read raw data

Imputation <- read_tsv("../Data/reference/Imputed.txt")

datalegend <- Imputation %>% select(c(1:4))
datagenotypes <- as.data.frame(Imputation %>% select(-c(1:4)))
datagenotypes <- apply(datagenotypes, 2, as.numeric)

## Get Isolates for each population
gambia <- datagenotypes[, -(grep("PH", colnames(datagenotypes)))]
cambodia <- datagenotypes[, grep("PH", colnames(datagenotypes))]

ncol(gambia)
ncol(cambodia)

## Pairwise Fst Gambia and Cambodia
fst <- fstsingle(cbind(gambia, cambodia), datalegend, ncol(gambia), ncol(cambodia))
fst$fst <- round(fst$fst, digits=4)

# Save into a excel file
fst %>% 
    distinct(POS, .keep_all = TRUE) %>% 
    write.table('pairwise_fst.xlsx', col.names = TRUE, 
                row.names = FALSE, quote = FALSE, sep = '\t')
