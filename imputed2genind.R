
###############################
## Convert vcf Column Format
## To freebayes format 
###############################

rm(list = ls())

args = (commandArgs(TRUE))

if(length(args) == 0)
{
    print("Convert GATK vcf format to Freebayes")
    print("Usage: gatk2freebayesformat.R imputed_data.txt")
    
    stop("Not enough arguments supplied ...")
}

imputed_data <- args[1]

library(data.table)
library(tidyverse)

Replace <- function(x){
    y <- as.matrix(x)
    y[y == 0] <- strrep(y[3], 2) 
    y[y == 1] <- strrep(y[4], 2)
    y
}

vcf <- as_tibble(read.table(imputed_data, header=TRUE, comment.char="", 
                            stringsAsFactors = FALSE, check.names= FALSE))

Data <- as_tibble(t(apply(vcf, 1, Replace)))

# Subset genotypes
snpgeno = as_tibble(t(Data[ , 5:ncol(Data)]))
loci <- vcf %>% select(POS) %>% pull()
names(snpgeno) <- loci

indv <- names(vcf)[-c(1:4)]
populations <- c(rep("Gambia", length(grep("PA", indv))), rep("Cambodia", length(grep("PH", indv))))

gam_cam_gen = df2genind(snpgeno, ploidy = 2, ind.names = indv, pop = populations, sep = "")
gam_cam_gen$tab[1:5, 1:10]

is.genind(gam_cam_gen)
bs.nc <- basic.stats(gam_cam_gen)
bs.nc
