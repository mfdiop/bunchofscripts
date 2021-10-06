# haplodise.R
#
# Author: Mouhamadou Fadel DIOP
# Date: 2021-10-06
#
# Purpose:
# Produce tajima's D excel file.
#
# ------------------------------------------------------------------

# NOTE - uncomment these lines to install packages as needed
# install.packages("data.table")
# install.packages("tictoc")
# install.packages("statip")

# Load packages
library(data.table)
library(tictoc)
library(statip)

options(scipen=999)

SelectionData = "/media/Data/Data/Documents_Karim/Fadel/Alfred/Fst/"

vcf = paste0(SelectionData, 'Cam_Gam14_15.vcf.gz')
Genotypes = paste0(SelectionData, 'Genotypes.txt')
AllelicDeph = paste0(SelectionData, 'AllelicD.txt')
GTexpression = '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n'
ADexpression = '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n'

system(sprintf("bcftools query -f'%s' %s > %s", GTexpression, vcf, Genotypes))  
system(sprintf("bcftools query -f'%s' %s > %s", ADexpression, vcf, AllelicDeph))

# ---------- Phased the final vcf file: FinalSelectionData.vcg.gz

genotypeData = fread(Genotypes, header = F)
allelicDepthData = fread(AllelicDeph, header = F)

first4Columns= subset(genotypeData, select=c(1:4))
genotypeData= as.data.frame(subset(genotypeData, select=-c(1:4)))
allelicDepthData= as.data.frame(subset(allelicDepthData, select=-c(1:4)))

write.table(first4Columns, paste0(SelectionData,"first4Columns.txt"), quote = FALSE, row.names = FALSE, col.names = F)  
FileName =paste0(SelectionData,'PhasedData.txt')

numberOfSimulation=1

genotypePhasing <-function (genotypeData, allelicDepth, first4Columns, numberOfSimulation, FileName)
{
    i=1
    while(i <= numberOfSimulation)
    {
        print(paste0("i = ", i), quote = FALSE)
        phasedData = matrix(9, nrow=dim(genotypeData)[1], ncol=dim(genotypeData)[2])
        for(j in 1:nrow(genotypeData))
        {
            k=1
            while(k<=ncol(genotypeData))
            {
                if(genotypeData[j,k]=='0/0')
                    phasedData[j,k]=0   
                else if(genotypeData[j,k]=='./.')
                    phasedData[j,k]='.'
                else if(genotypeData[j,k]=='1/1')
                    phasedData[j,k]=1
                else if(genotypeData[j,k]=='0/1')
                {
                    target=as.integer(unlist(strsplit(allelicDepthData[j,k], ',')))
                    if(sum(target)>0)
                    {
                        if(target[1]<target[2])
                        {
                            minor=target[1]
                            maf=minor/sum(target)
                            phasedData[j,k]=rbern(1,maf)
                        }
                        else
                        {
                            minor=target[2] 
                            maf=minor/sum(target)
                            phasedData[j,k]=rbern(1,maf)
                        }
                    }
                    else
                    {
                        maf=target[1]
                        phasedData[j,k]=rbern(1,maf)
                    }
                    
                }
                k=k+1
            }
        }
        
        i=i+1
    }
    
    phasedData = as.data.frame(cbind(first4Columns,phasedData))
    
    write.table(phasedData, FileName, quote = FALSE, row.names = FALSE, col.names = F)
}

system.time(genotypePhasing(genotypeData, allelicDepthData, first4Columns,numberOfSimulation,FileName))


#----------- Impute missing data ------------
PhasedData = fread(FileName, header = FALSE)
PhasedData = subset(PhasedData, select = -c(1:4))

Name = paste0(SelectionData, 'Cam_Gam14_15_ImputedData.txt')

system(paste0('rm -rf ', SelectionData, 'Imputation_Folder/'))
system(paste0('mkdir -p ', SelectionData, 'Imputation_Folder/'))

imputeMissingGenotypes = function(PhasedData, first4Columns, numberOfSimulation, Name)
{
    maf = vector(mode = "numeric", length = dim(PhasedData)[1])
    PhasedData =as.data.frame(PhasedData)
    print(paste0("computing the MAF before simulation"), quote = FALSE)
    for(i in 1:dim(PhasedData)[1])
    {
        print(paste0("i=",i), quote = FALSE)
        count1=0
        count2=0
        for(j in 1:dim(PhasedData)[2])
        {
            if(PhasedData[i,j]=='.')
                count1=count1+1
            if(PhasedData[i,j]=='1')
                count2=count2+1
        }
        validCount = dim(PhasedData)[2] - count1
        count3 = validCount-count2
        if(count2 < count3)
            maf[i]=count2/validCount
        else if(count2 > count3)
            maf[i]=count3/validCount
        else
            maf[i]=count2/validCount
    }
    
    print(paste0("Imputing the missing genotypes"), quote = FALSE)
    s=1
    correlationCoef = vector(mode = "numeric", length = numberOfSimulation)
    while(s <= numberOfSimulation)
    {
        imputedGenotypeMatrix = matrix(-9, nrow = dim(PhasedData)[1], ncol = dim(PhasedData)[2])
        for(j in 1:dim(PhasedData)[1])
        {
            k=1
            while(k<=dim(PhasedData)[2])
            {
                if(PhasedData[j,k]=='0')
                    imputedGenotypeMatrix[j,k]=0   
                else if(PhasedData[j,k]=='1')
                    imputedGenotypeMatrix[j,k]=1
                else if(PhasedData[j,k]=='.')
                    imputedGenotypeMatrix[j,k]=rbern(1,maf[j])
                k=k+1
            }
        }
        
        print(paste0("computing the MAF after simulation"), quote = FALSE)
        mafAfter=vector(mode = "numeric", length = dim(imputedGenotypeMatrix)[1])
        count4=rowSums(imputedGenotypeMatrix)
        count5=dim(imputedGenotypeMatrix)[2]-count4
        for(l in 1:length(mafAfter))
         {
          if(count4[l] < count5[l])
             mafAfter[l] = count4[l]/dim(imputedGenotypeMatrix)[2]
           else
             mafAfter[l] = count5[l]/dim(imputedGenotypeMatrix)[2]
           
        }
         correlationCoef[s] = cor(maf, mafAfter)
        
        imputedGenotypeMatrix = cbind(first4Columns, as.data.frame(imputedGenotypeMatrix))
        outFileName = paste0(SelectionData, 'Imputation_Folder/Imputed_', s, ".txt")
        write.table(imputedGenotypeMatrix, outFileName, quote = FALSE, row.names = FALSE, col.names = FALSE)
        print(paste0("simulation = ", s), quote = FALSE)
        s=s+1
    }
    
    which(correlationCoef==max(correlationCoef))
    print('----- Save The Imputed File ------')

    imputedGenotypeMatrix = cbind(first4Columns, as.data.frame(imputedGenotypeMatrix))
    write.table(imputedGenotypeMatrix, Name, quote = FALSE, row.names = FALSE, col.names = FALSE)
}

system.time(imputeMissingGenotypes(PhasedData, first4Columns, numberOfSimulation, Name))

Imputation <- fread(paste0(SelectionData, 'Imputation_Folder/Imputed_75.txt'), header = FALSE)

system(paste0("bcftools query -l ", paste0(SelectionData, 'Cam_Gam14_15.vcf.gz'), " > ", paste0(SelectionData, 'SampleIDs.txt')))
samples <- scan(paste0(SelectionData, 'SampleIDs.txt'), what = 'character')

names(Imputation) <- c("CHROM", "POS", "REF", "ALT", samples)

write.table(Imputation, Name, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')


