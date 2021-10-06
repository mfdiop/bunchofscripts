# filter.R
#
# Author: Mouhamadou Fadel DIOP
# Date: 2021-10-06
#
# Purpose:
# Remove missing data on sample and locus
# Produce filtered VCF file.
#
# ------------------------------------------------------------------

# NOTE - uncomment these lines to install packages as needed
# install.packages("data.table")
# install.packages("tictoc")

arguments <- commandArgs(trailingOnly = TRUE)

Dir <- arguments[1]

# Load packages
library(data.table)
library(tictoc)
options(scipen=999)

source("Functions.R")

#=============================================
#=== Extracting genotypes from the raw data
#=============================================
setwd(Dir)
vcf <- arguments[2]
Genotypes <- '../process/Genotypes.txt'

expression <- '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n'

system(sprintf("bcftools query -f'%s' %s > %s", expression, vcf, Genotypes))

#=============================================
#-- Computing the missingness on SNPs and samples
#-- Get the sample names from the VCF file
#=============================================

sampleList <- '../process/SampleList.txt'
system(sprintf("bcftools query -l %s > %s", vcf, sampleList))
AfricanSamplesList <- fread(sampleList, header = FALSE)

#=============================================
#---- put the first four columns in a variable
#=============================================
Genotype <- fread(Genotypes, header = FALSE)

first4Column <- subset(Genotype, select=c(1:4))
Genotype <- subset(Genotype, select=-c(1:4))

#===================
## Remove invariants
#==================
varsnp <- apply(Genotype, 1, triquad)
keepvar <- which(varsnp != 1)

data2 <- cbind(first4Column, Genotype)
data2 <- data2[keepvar,]
geno <- data2[,5:ncol(data2)]
first4Column <- data2[,1:4]
nrow(geno)
ncol(geno)

#======================================
#--- computing missingness on SNPs
#======================================

tic();snpMissingness <- computeMissingnessOnSNPs(geno);toc()

geno <- cbind(first4Column, geno, snpMissingness)

geno <- geno[which(geno$snpMissingness <= 0.2), ]
first4Column <- subset(geno, select = c(1:4))
geno <- subset(geno, select = -c(1:4, ncol(geno)))

#================================
#=== Compute Missingness on samples
#================================
tic(); sampleMissingness <- computeMissingness(geno); toc()

geno <- t(geno)
rownames(geno) <- AfricanSamplesList$V1
geno <- as.data.frame(geno)
geno$Missingness <- sampleMissingness

index <- which(geno$Missingness > 0.15)

geno <- subset(geno, select = -c(ncol(geno)))
geno <- as.data.frame(t(geno))

if(!is_empty(index)) geno <- geno[, -index]

#===================
## Remove invariants
#==================
varsnp <- apply(geno, 1, triquad)
keepvar <- which(varsnp != 1)

data2 <- cbind(first4Column, geno)
data2 <- data2[keepvar,]
geno <- data2[,5:ncol(data2)]
first4Column <- data2[,1:4]
nrow(geno)
ncol(geno)

## Save positions and isolates to keep for downstream analysis
write.table(colnames(geno), '../process/samplesTokeep.txt', 
            col.names = FALSE, row.names = FALSE, quote = FALSE)

write.table(first4Column[,1:2], "../process/snpsTokeep.txt", 
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = '\t')

#=========================================================
#---- Removing SNPs and Isolates to discard from vcf file
#========================================================
samplesToKeep <- "../process/samplesTokeep.txt"
snpsToKeep <- "../process/snpsTokeep.txt"

filtered.vcf <- paste0("../process/", gsub(".vcf.gz", ".filtered", vcf))

system(paste0("vcftools --gzvcf ", vcf,
              " --keep ", samplesToKeep, 
              " --positions ", snpsToKeep,
              " --not-chr Pf3D7_API_v3", 
              " --recode --recode-INFO-all --out ", 
			  filtered.vcf))

file.remove(Genotypes, samplesToKeep, snpsToKeep)

system(paste0("mv ", filtered.vcf, ".recode.vcf ", filtered.vcf, ".vcf"))
system(paste0("bgzip ", filtered.vcf, ".vcf"))
system(paste0("tabix ", filtered.vcf, ".vcf.gz"))