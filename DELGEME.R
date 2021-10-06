
library(data.table)
library(tictoc)
library(statip)

options(scipen=999)


SelectionData = "/home/fadel/DELGEME/"

#-------------------- Splitting the data by chromosome------------------
splitVCF = function(vcfileName)
{
  outFolder = gsub('.vcf.gz','',basename(vcfileName))
  infolder = dirname(vcfileName)
  setwd(infolder)
  system(sprintf("mkdir %s", paste0(outFolder,'/')))
  system(sprintf("mv %s %s", basename(vcfileName), paste0(outFolder,'/')))
  contigs = c("Pf3D7_01_v3", "Pf3D7_02_v3", "Pf3D7_03_v3", "Pf3D7_04_v3", "Pf3D7_05_v3", "Pf3D7_06_v3", "Pf3D7_07_v3",
              "Pf3D7_08_v3", "Pf3D7_09_v3", "Pf3D7_10_v3", "Pf3D7_11_v3", "Pf3D7_12_v3", "Pf3D7_13_v3", "Pf3D7_14_v3")
  outdir = paste0(infolder, '/', outFolder)
  rr = paste0(outdir, '/', basename(vcfileName))
  system(sprintf("bcftools index -t %s", rr))
  for(i in 1:length(contigs))
  {
    r=sprintf(contigs[i])
    out = paste0(outdir, '/', "Chromosome",i,'.vcf.gz')
    system(sprintf("bcftools view -r %s -o %s -O z %s", r, out, rr))
  }
}

#------------ Extracting genotypes from the raw data: AfricanPop_M001_CDS.vcf.gz

vcf="/home/fadel/DELGEME/Candidates.vcf.gz"

vcf="/home/fadel/DELGEME/GenomeSNP-maf001.vcf.gz"
SelectionData='/home/fadel/DELGEME/Target/PopGenStats/IBD'

Genotypes = paste0(SelectionData, '/', 'Genotypes.txt')

expression = '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n'

system(sprintf("bcftools query -f'%s' %s > %s",expression, vcf, Genotypes))

#------- Computing the missingness on SNPs and samples from Candidates.vcf.gz
#---- --------get the sample names from the VCF file

#sampleList = "/home/fadel/DELGEME/SampleList.txt"
sampleList = paste0(SelectionData, '/', 'SampleList.txt')
system(sprintf("bcftools query -l %s > %s", vcf, sampleList))
AfricanSamplesList = fread(sampleList, header = FALSE)

#---- put the first four columns in a variable
Genotype = fread(paste0(SelectionData, '/', 'Genotypes.txt'), header = FALSE)

first4Column = subset(Genotype, select=c(1:4))

Genotype = subset(Genotype, select=-c(1:4))

#-------- computing missingness on samples
computeMissingness = function(dataFrame)
{
  dataFrame = t(dataFrame)
  #dataFrame = as.matrix(dataFrame)
  missingness = vector(mode = "numeric", length = dim(dataFrame)[1])
  for(i in 1:dim(dataFrame)[1])
  {
    m=0
    for(j in 1:dim(dataFrame)[2])
    {
      if(dataFrame[i,j] == './.')
        m=m+1
    }
    missingness[i] = m/(dim(dataFrame)[2])
  }
  return(missingness)
}
tic();sampleMissingness = computeMissingness(Genotype);toc()

write.table(sampleMissingness, paste0(SelectionData, '/', 'sampleMissingness.txt'), col.names = FALSE, row.names = FALSE, quote = FALSE)

Genotype = t(Genotype)

rownames(Genotype) = AfricanSamplesList$V1

Genotype = as.data.frame(Genotype)

Genotype$Missingness = sampleMissingness

samplesToBeDiscarded = rownames(Genotype[which(Genotype$Missingness > 0.2),])

write.table(samplesToBeDiscarded, paste0(SelectionData, '/', 'samplesToBeDiscarded.txt'), col.names = FALSE, row.names = FALSE, quote = FALSE)

Genotype = subset(Genotype, select = -c(ncol(Genotype)))

Genotype = as.data.frame(t(Genotype))

#-------- computing missingness on SNPs
computeMissingnessOnSNPs = function(dataFrame)
{
  #dataFrame = t(as.matrix(subset(dataFrame, select=-c(1))))
  dataFrame = as.matrix(dataFrame)
  missingness = vector(mode = "numeric", length = dim(dataFrame)[1])
  for(i in 1:dim(dataFrame)[1])
  {
    m=0
    for(j in 1:dim(dataFrame)[2])
    {
      if(dataFrame[i,j] == './.')
        m=m+1
    }
    missingness[i] = m/(dim(dataFrame)[2]-1)
  }
  return(missingness)
}
tic();snpMissingness = computeMissingnessOnSNPs(Genotype);toc()

Genotype = cbind(first4Column, Genotype, snpMissingness)

snpsToBeDiscarded = Genotype[which(Genotype$snpMissingness > 0.2),]

snpsToBeDiscarded = subset(snpsToBeDiscarded, select = c(1:2))

write.table(snpsToBeDiscarded, paste0(SelectionData,"/snpsToBeDiscarded.txt"), col.names = FALSE, row.names = FALSE, quote = FALSE)

#-------- Plotting the missingness
pdf("/home/fadel/DELGEME/Missingness.pdf")
par(mfrow=c(2,1))
hist(snpMissingness,100,ylim = c(0,500))
hist(sampleMissingness,100, ylim = c(0,500))
dev.off()

#-------------------- Removing the samples to discard from the vcf file
samplesToBeDiscarded = paste0(SelectionData,"/samplesToBeDiscarded.txt")
vcf1 = paste0(SelectionData,"/WholeGenome-maf0.01_filteredMissingness.vcf.gz")
system(sprintf("bcftools view -S^%s --threads 15 -o %s -O z %s", samplesToBeDiscarded, vcf1, vcf))

#-------------------- Removing SNPs to discard after removing uneeded samples from Target.vcf.gz
#-------------------- Spliting the vcf file to avoid duplicate loci position on chromosomes
tic();splitVCF(vcfileName=vcf1);toc()

#-------------------- Getting the positions to remove by chromosome
d = paste0(dirname(vcf1),'/',gsub('.vcf.gz','', basename(vcf1)))
for(i in 1:length(contigs))
{
  c = snpsToBeDiscarded[which(snpsToBeDiscarded$V1 == contigs[i]),]
  write.table(t(c$V2), paste0(d,'/','Chrom',i,'.txt'), quote = FALSE, row.names = FALSE, col.names = FALSE)
}

#-------------------- Concatenate the obtained VCF files 
#path = "/home/fadel/DELGEME/Target/"
#setwd(path)
#files = paste0(path,'/','list.txt')
files='/home/fadel/DELGEME/Target/PopGenStats/IBD/WholeGenome-maf0.01_filteredMissingness/list.txt'
#out = paste0(SelectionData,"/WholeGenome-maf0.01_filteredMissingness.vcf.gz")
#system(sprintf("bcftools concat -f %s -o %s -O z --threads 15", files, out))
#system(sprintf("tabix %s", out))

out="/home/fadel/DELGEME/Target/PopGenStats/IBD/WholeGenome-maf0.01_filteredMissingness/Target.vcf.gz"
system(paste0("bcftools concat -f ",files," --threads 15 -Oz -o ",out))
#------------------- Extract genotypes and allelic deph from the final vcf file
#path="/home/fadel/DELGEME/Target"
vcf= paste0(SelectionData,'/WholeGenome-maf0.01_filteredMissingness.vcf.gz')
Genotypes = paste0(SelectionData,'/Genotypes.txt')
AllelicDeph = paste0(SelectionData,'/AllelicD.txt')
GTexpression = '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n'
ADexpression = '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n'

system(sprintf("bcftools query -f'%s' %s > %s",GTexpression, vcf, Genotypes))  
system(sprintf("bcftools query -f'%s' %s > %s",ADexpression, vcf, AllelicDeph))

# ---------- Phased the final vcf file: FinalSelectionData.vcg.gz

genotypeData = fread(Genotypes, header = F)
allelicDepthData = fread(AllelicDeph, header = F)


first4Columns= subset(genotypeData, select=c(1:4))
genotypeData= as.data.frame(subset(genotypeData, select=-c(1:4)))
allelicDepthData= as.data.frame(subset(allelicDepthData, select=-c(1:4)))


write.table(first4Columns, paste0(SelectionData,"/first4Columns.txt"), quote = FALSE, row.names = FALSE, col.names = F)  
FileName =paste0(SelectionData,'/PhasedData.txt')

numberOfSimulation=100

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

tic();genotypePhasing(genotypeData, allelicDepthData, first4Columns,numberOfSimulation,FileName);toc()


#---- Create my metadata with samples ID, country and sites
listSamples= scan('SampleList.txt', what = character())

meta=fread("/home/fadel/mrcuser/Fadel/New_process/SNPs_files/working_file/DELGEME-ANALYSIS/BALANCING-SELECTION-STUDY/Analysis/SelectionData/FinalData/RawData/MetaData/AfricanSamplesMetaData.txt", fill = T)

newcol1 = NULL
for (i in 1:length(listSamples))
{
  for (j in 1:nrow(meta))
  {
    if (listSamples[i] == meta$sample[j])
    {
      if(i==1)
      {
        newcol1 = cbind(newcol1, meta[j, c(1,9)])
        
      }
      newcol1 = rbind(newcol1,meta[j, c(1,9)])

    }  
  }
}

names(newcol1)=c("SampleID", "Country")

mau= as.data.frame(grep('FP', listSamples, value = T)) 
mau$country='mauritanie'
colnames(mau)=c('sample', 'country')

write.table(newcol1, 'MetaAfricanSamples.txt', col.names = T, row.names = F, quote = F )
write.table(mau, 'Mauritanie.txt', col.names = T, row.names = F, quote = F )

#---------- Get missing sample ID in meta data
#---------- Samples which are in the list and not in the meta data
listSamples= fread('SampleList.txt', header = F)
Samples=newcol1$SampleID

MissingSamples = listSamples[!which(listSamples$V1 %in% Samples),]
write.table(MissingSamples, 'missing.txt', col.names = F, row.names = F, quote = F )

#----------- Impute missing data ------------
PhasedData=fread('PhasedData.txt', header = F)
PhasedData=subset(PhasedData, select = -c(1:4))

Name=paste0(path,'/','ImputedData.txt')

system(sprintf('mkdir Imputation_Folder/'))

imputeMissingGenotypes = function(PhasedData, first4Columns, numberOfSimulation, Name,path)
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
  #correlationCoef = vector(mode = "numeric", length = numberOfSimulation)
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
    
    #print(paste0("computing the MAF after simulation"), quote = FALSE)
    #mafAfter=vector(mode = "numeric", length = dim(imputedGenotypeMatrix)[1])
    #count4=rowSums(imputedGenotypeMatrix)
    #count5=dim(imputedGenotypeMatrix)[2]-count4
    #for(l in 1:length(mafAfter))
   # {
    #  if(count4[l] < count5[l])
   #     mafAfter[l] = count4[l]/dim(imputedGenotypeMatrix)[2]
   #   else
   #     mafAfter[l] = count5[l]/dim(imputedGenotypeMatrix)[2]
   #   
    #}
   # correlationCoef[s] = cor(maf, mafAfter)
    
    #imputedGenotypeMatrix = cbind(first4Columns, as.data.frame(imputedGenotypeMatrix))
    #outFileName = paste0(path, "/Imputation_Folder/Imputed_",s,".txt")
   # write.table(imputedGenotypeMatrix, outFileName, quote = FALSE, row.names = FALSE, col.names = FALSE)
   # print(paste0("simulation = ", s), quote = FALSE)
    s=s+1
  }
  #return(which(correlationCoef==max(correlationCoef)))
  print('----- Save The Imputed File ------')
  
  Name=paste0(path,'/',Name)
  imputedGenotypeMatrix = cbind(first4Columns, as.data.frame(imputedGenotypeMatrix))
  
  write.table(imputedGenotypeMatrix, Name, quote = FALSE, row.names = FALSE, col.names = FALSE)
}

tic();imputeMissingGenotypes(PhasedData, first4Columns,numberOfSimulation,Name, path);toc()

Imputation=fread('ImputedData.txt', header = F)
header=fread('header.txt', header = T)

names(Imputation)=colnames(header)
names(PhasedData)=colnames(header)
write.table(Imputation, 'ImputedData.txt', col.names = T, row.names = F, quote = F)
write.table(PhasedData, 'PhasedData.txt', col.names = T, row.names = F, quote = F)

#---------- Split the data per Chromosome (gene) ----------------
SplitDataByChromosome <- function(targetDataSet, rawDataPath)
{
  first4Columns= subset(targetDataSet, select=c(1:4))
  chrom = unique(targetDataSet$`#CHROM`)
  for(i in 1:length(chrom))
  {
    cat(paste("i= ", i))
    d = targetDataSet[which(targetDataSet$`#CHROM`== chrom[i]),]
    dd = as.data.frame(subset(d, select=-c(1:4)))
    dd[, ] = lapply(dd[, ], as.character)
    write.table(dd, paste0(rawDataPath, "chrom",i,".txt"), col.names = T, row.names = FALSE, quote = FALSE)
  }
}

#-------------------------------------------------------------------------------
#------ Estimate the selection parameters using PopGenome
#--- Tajima D, Fu & Li, Theta Watterson ect... for each gene and each population
library(data.table)
library(PopGenome)
#--- vcf files should be zipped and Indexed
Outpout='/home/fadel/DELGEME/Target/PopGenStats/TestOfSelection/PopGenome-Output'

chroms=c("Pf3D7_11_v3","Pf3D7_12_v3","Pf3D7_01_v3","Pf3D7_03_v3","Pf3D7_03_v3", "Pf3D7_07_v3","Pf3D7_10_v3", "Pf3D7_08_v3",
         "Pf3D7_02_v3","Pf3D7_09_v3","Pf3D7_10_v3","Pf3D7_14_v3","Pf3D7_04_v3","Pf3D7_08_v3","Pf3D7_14_v3","Pf3D7_13_v3")

frompos= c(1293856, 656738, 558201, 224160, 221323, 1358055, 1399195, 1307053, 796752, 1201812, 896865, 3188888, 1082363, 1310655, 3193285, 1464895)
topos= c(1295724, 660286, 558708, 225793, 222516, 1362929, 1402895, 1307924, 801583, 1206974, 903608, 3190548, 1084150, 1316883, 3199478, 1466819)

setwd('/home/fadel/DELGEME/Target/PopGenStats/TestOfSelection/PopGenome-Input/TRAP/')
NewDF=NULL
files=list.files('.', pattern = '.vcf.gz$', full.names = T)
i=16
for (j in 1:length(files)) 
{
  Nom=gsub('_GT.vcf.gz', '',basename(files[j]))
  print(paste0('File= ', Nom), quote = F)
  
  if(j==1)
  {
    chrom1 = readVCF(files[j], 10000, chroms[i], frompos[i], topos[i], include.unknown=T)
    
    GENOME.class <- neutrality.stats(chrom1, FAST=TRUE)
    
    TestNeutrality=as.data.frame(cbind(subset(round(get.neutrality(GENOME.class)[[1]],3), select=c(1:2,4,5))))
    
    f= as.data.frame(round(GENOME.class@theta_Watterson,3))
    
    NewDF= as.data.frame(cbind(TestNeutrality$n.segregating.sites,TestNeutrality$Tajima.D,TestNeutrality$Fu.Li.F,TestNeutrality$Fu.Li.D, f))
    rownames(NewDF)=Nom
  }
  else
  {
    chrom1 = readVCF(files[j], 10000, chroms[i], frompos[i], topos[i], include.unknown=T)
    
    GENOME.class <- neutrality.stats(chrom1, FAST=TRUE)
    
    TestNeutrality=as.data.frame(cbind(subset(round(get.neutrality(GENOME.class)[[1]],3), select=c(1:2,4,5))))
    
    f= as.data.frame(round(GENOME.class@theta_Watterson,3))
    
    TestNeutrality= as.data.frame(cbind(TestNeutrality$n.segregating.sites,TestNeutrality$Tajima.D,TestNeutrality$Fu.Li.F,TestNeutrality$Fu.Li.D, f))
    rownames(TestNeutrality)=Nom
    NewDF=rbind(NewDF,TestNeutrality)
  }
  
}

names(NewDF)= c('Segregating_Sites','Tajima_D', 'Fu&Li_F', 'Fu&Li_D', 'Theta_Watterson')

write.table(NewDF, paste0(Outpout, 'trap.xls'), col.names = T, row.names = T, quote = F, sep = '\t')

#-------------------------------------------------------------------------------
#------ Estimate the Nucleotide and Haplotype Diversities using PopGenome
#--- for each gene and each population: see the R script Fst-Popgenome
#*********************************************************************
#************* Compute pairwise Fst (Fst between populations)
#****** using hierfstat package. Takes as input the imputed data
#***** convert to diploid genotypes 0/0 or 1/1

library(data.table)
library(hierfstat)
library(pegas)
library(tictoc)

ImputedFile=fread('/home/fadel/DELGEME/Target/ImputedData.txt',header = T)

first4Columns = subset(ImputedFile, select=c(1:4))
Genotype= as.matrix(subset(ImputedFile, select=-c(1:4)))
tic()
Data = matrix('z', nrow=dim(Genotype)[1], ncol=dim(Genotype)[2])
tic()
for(l in 1:nrow(Genotype))
{
  k=1
  while(k<=ncol(Genotype))
  {
    if(Genotype[l,k] ==0)
      Data[l,k]='0/0' 
    else
      Data[l,k]='1/1'
    k=k+1
  }
}
toc()

Data = cbind(first4Columns, as.data.frame(Data))

write.table(Data, file = '/home/fadel/DELGEME/Target/ImputedGenotypes.txt', sep = "\t", row.names = FALSE, col.names = F,quote = FALSE)
#------------------------------------------------------------------------------------------


files=list.files('/home/fadel/DELGEME/Target/PopGenStats/Fst/Fst-hierfstat', pattern = '.txt', full.names = T)

out='/home/fadel/DELGEME/Target/PopGenStats/Fst/Fst-hierfstat/'

Samples=fread('/home/fadel/DELGEME/Target/MetaAfricanSamples.txt', header = F)

for (i in 1:length(files)) 
{
  Nom=gsub('.txt', '', basename(files[i]))
  cat(paste0('File = ', Nom, '\n'))
  
  d <- read.table(files[i], header = F)
  f<- d[,-c(1:4)]
  
  loci= as.character(d$V2)
  
  population= as.character(Samples$V2)
  
  indv.names=as.character(Samples$V1)
  indv.names <- gsub("\\.", "_", indv.names)
  
  ff=t(f)
  
  #---------------- WE SHOULD CONVERT THE DATAFRAME OR THE MATRIX TO A GENIND OBJECT
  
  Mydata <- df2genind(ff, ncode = 2, ind.names = indv.names, loc.names = loci, pop = population, ploidy = 2, type ="codom")
  
  #-------- Compute Pairwise Fst using hierfstat package
  #---- Error in fstat(Mydata, fstonly = FALSE) : not implemented for non-diploid genotypes
  #---- Works on diploid data
  
  mat.pairwise.fst= pairwise.fst(Mydata, pop = population, res.type ="matrix")
  
  write.table(mat.pairwise.fst, paste0(out,Nom,'.xls'), row.names = T, col.names = T, quote = F)
  
  print('*******************************************')
}


#------------------------------------------------------------------------------------------------------------
#---- use isorelate to infer pairwise IBD: Use the genotype data containing all chromosomes and all populations
#-------- First lets design the input files for RealMcCOIL server
#-------- Got the snp list from Fadel
#setwd('/home/fadel/DELGEME/Target')
snpToBeUsed = fread("/home/fadel/DELGEME/Target/RealMcCoil/WTSI_BarcodeSNPlist.txt", header = FALSE)
delgemeGenotype = fread(paste0(SelectionData,'/Genotypes.txt'), header = F)
contigs = c("Pf3D7_01_v3", "Pf3D7_02_v3", "Pf3D7_03_v3", "Pf3D7_04_v3", "Pf3D7_05_v3", "Pf3D7_06_v3", "Pf3D7_07_v3",
            "Pf3D7_08_v3", "Pf3D7_09_v3", "Pf3D7_10_v3", "Pf3D7_11_v3", "Pf3D7_12_v3", "Pf3D7_13_v3", "Pf3D7_14_v3")

for(i in 1:length(contigs))
{
  data = delgemeGenotype[which(delgemeGenotype$V1==contigs[i]),]
  snp = snpToBeUsed[which(snpToBeUsed$V1==contigs[i]),]
  if(i==1)
    d = data[which(data$V2 %in% snp$V2),]
  else
    d = rbind(d, data[which(data$V2 %in% snp$V2),])
}
write.table(subset(d, select = c(1:2)), paste0(SelectionData,"/RealMcCoil_SnpsList.txt"), col.names = FALSE, row.names = FALSE, quote = FALSE)

dd = as.matrix(subset(d, select = -c(1:4)))
first4Column = subset(d, select = c(1:4))
modifiedD = matrix(data='z', nrow = dim(dd)[1], ncol = dim(dd)[2])
for(i in 1:dim(dd)[1])
{
  for(j in 1:dim(dd)[2])
  {
    if(dd[i,j]=='0/0')
      modifiedD[i,j]=first4Column$V3[i]
    if(dd[i,j]=='1/1')
      modifiedD[i,j]=first4Column$V4[i]
    if(dd[i,j]=='0/1')
      modifiedD[i,j]='X'
    if(dd[i,j]=='./.')
      modifiedD[i,j]='N'
  }
}

modifiedD = t(modifiedD)
delgemeRealMcCoilSNP = data.frame(sample=character(), snps=character(), stringsAsFactors = FALSE)
samples=fread(paste0(SelectionData,'/samples_list.txt'), header = F)
for(i in 1:nrow(modifiedD))
{
  delgemeRealMcCoilSNP[i,1] = samples$V1[i]
  delgemeRealMcCoilSNP[i,2] = paste(modifiedD[i,], collapse='')
}
write.table(delgemeRealMcCoilSNP, paste0(SelectionData,"/delgemeRealMcCoilSNP.txt"), col.names = TRUE, row.names = FALSE, quote = FALSE)

#-------- making ped and map files for isoRelate
#--- That ped file is created from the phased data (on which the missing genotypes have not been imputed yet)
phasedData = fread(paste0(SelectionData,"/PhasedData.txt"), header = FALSE)
first4Columns = subset(phasedData, select = c(1:4))
phasedData = subset(phasedData, select = -c(1:4))
#outDir = "/home/fadel/DELGEME/Target/PopGenStats/IBD"
phasedData = as.matrix(phasedData)
phasedData[phasedData=='1']=2
phasedData[phasedData=='0']=1
phasedData[phasedData=='.']=0
phasedData = cbind(first4Columns,as.data.frame(phasedData))
write.table(phasedData, paste0(SelectionData,'/','PhasedData.txt'), row.names = FALSE, col.names = FALSE, quote = FALSE)

#------- Create Ped file ---------
createPed = function(geno, phenotype, outputPath)
{
  library(data.table)
  myRaw = fread(geno, header = FALSE)
  contigs = c("Pf3D7_01_v3", "Pf3D7_02_v3", "Pf3D7_03_v3", "Pf3D7_04_v3", "Pf3D7_05_v3", "Pf3D7_06_v3","Pf3D7_07_v3",
              "Pf3D7_08_v3", "Pf3D7_09_v3", "Pf3D7_10_v3", "Pf3D7_11_v3", "Pf3D7_12_v3", "Pf3D7_13_v3", "Pf3D7_14_v3")
  for(k in 1:length(contigs))
  {
    print(paste0("Creating ped file for chromosome ",k), quote = FALSE)
    data = myRaw[which(myRaw$V1==contigs[k]),]
    genotypeData = subset(data, select=c(5:ncol(data)))
    genotypeData = as.data.frame(t(genotypeData))
    #rm(myRaw, data)
    
    mat1 = as.matrix(genotypeData)
    #mat1 = mat1+1
    genotype = matrix(-9 , nrow = dim(mat1)[1], ncol = (dim(mat1)[2]*2) )
    for(l in 1:dim(mat1)[1])
    {
      j = 1
      for(i in 1:dim(mat1)[2])
      {
        genotype[l,j] = mat1[l,i]
        genotype[l, j+1] = mat1[l,i]
        j = j+2
      }
    }
    genotype = as.data.frame(genotype)
    
    PED6 = data.frame(data = phenotype, FatherID = 0, MotherID = 0, Sex = 0, Phenotype = -9, stringsAsFactors = FALSE)
    ped = cbind(PED6, genotype)
    write.table(ped, paste0(outputPath,'/Delgeme_Chrom',k,'.ped'), row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
}

phasedData = paste0(SelectionData,'/PhasedData.txt')
phenotypeData = fread(paste0(SelectionData,'/Phenotypes.txt'), header = F)
outDir='/home/fadel/DELGEME/Target/PopGenStats/IBD/PedMapFiles'
tic("Time to create ped file"); ped = createPed(phasedData, phenotypeData, outDir);toc()

phasedData = fread(paste0(SelectionData,'/PhasedData.txt'), header = FALSE)
first4Columns = subset(phasedData, select = c(1:4))
#------- Create Map file ---------
create_map = function(first4Columns, outputPath)
{
  contigs = c("Pf3D7_01_v3", "Pf3D7_02_v3", "Pf3D7_03_v3", "Pf3D7_04_v3", "Pf3D7_05_v3", "Pf3D7_06_v3","Pf3D7_07_v3",
              "Pf3D7_08_v3", "Pf3D7_09_v3", "Pf3D7_10_v3", "Pf3D7_11_v3", "Pf3D7_12_v3", "Pf3D7_13_v3", "Pf3D7_14_v3")
  
  for(i in 1:length(contigs))
  {
    c1 = first4Columns[which(first4Columns$V1 == contigs[i]),]$V1
    c4 = first4Columns[which(first4Columns$V1 == contigs[i]),]$V2
    c2 = paste0(contigs[i],':',c4)
    c3 = c4/12000 # Get positions in centiMorgan
    map = data.frame(chrom = c1, post = c2, GD = c3, pos = c4, stringsAsFactors = FALSE)
    write.table(map, paste0(outputPath,'/','Delgeme_Chrom',i,'.map'), row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
}
tic("Time to create map file"); map = create_map(first4Columns, outDir);toc()

#--- combining ped and map files into one ped and map file
peds = list.files(outDir, pattern = '.ped', full.names = TRUE)
peds = c(peds[6:14], peds[1:5])
maps = list.files(outDir, pattern = '.map', full.names = TRUE)
maps = c(maps[6:14], maps[1:5])

for(i in 1:length(maps))
{
  if(i==1)
    delgemeMap = fread(maps[i], header = FALSE)
  else
    delgemeMap = rbind(delgemeMap, fread(maps[i], header = FALSE))
}
write.table(delgemeMap, paste0(outDir,'/','Delgeme.map'), row.names = FALSE, col.names = FALSE, quote = FALSE)

for(i in 1:length(peds))
{
  if(i==1)
    delgemePed = fread(peds[i], header = FALSE)
  else
    delgemePed = cbind(delgemePed, subset(fread(peds[i], header = FALSE), select = -c(1:6)))
}
write.table(delgemePed, paste0(outDir,'/','Delgeme.ped'), row.names = FALSE, col.names = FALSE, quote = FALSE)


#delgemePed = fread(paste0(outDir,'/','Delgeme.ped'), header = FALSE)
#delgemeMap = fread(paste0(outDir,'/','Delgeme.map'), header = FALSE)
delgemeCOI = fread("/home/fadel/DELGEME/Target/PopGenStats/IBD/MOI.txt", header = FALSE)
delgemePed$V5 = delgemeCOI$V2
#delgemePed2 = cbind(delgemePed$V2,delgemePed)
#delgemePed2 = subset(delgemePed2, select = -c(3))
write.table(delgemePed, paste0(outDir,'/','Delgeme.ped'), col.names = FALSE, row.names = FALSE, quote = FALSE)
ped.map = list(delgemePed, delgemeMap)

#------------ Running isoRelate
my_genotypes = getGenotypes(ped.map = ped.map, reference.ped.map = NULL, maf = 0.01, isolate.max.missing = 0.2, snp.max.missing = 0.2, chromosomes = NULL, input.map.distance = "cM", reference.map.distance = "cM")


#-------------------------------------------------------------
#---------- COMPUTE t-SNE FOR EACH GENE ACROSS ALL POPULATIONS
#----- Use the impute genotypes data for each gene

files=list.files('/home/fadel/DELGEME/Target/PopGenStats/TestOfSelection/PopGenome-Fst_Hap_Nuc_diversity/Input', 
                 pattern = '.vcf.gz$', full.names = T)

#------------------- Extract genotypes and allelic deph from the final vcf file
PathForGenotypeData='/home/fadel/DELGEME/Target/PopGenStats/Rtsne/Input/'

GTexpression = '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n'
ADexpression = '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n'

for (i in 1:length(files)) 
{
  print(paste0('File is = ', files[i]))
  Nom=gsub('.vcf.gz','', basename(files[i]))
  
  Genotypes = paste0(PathForGenotypeData, Nom,'_GT.txt')
  AllelicDeph = paste0(PathForGenotypeData, Nom,'_AD.txt')
  
  system(sprintf("bcftools query -f'%s' %s > %s",GTexpression, files[i], Genotypes))  
  system(sprintf("bcftools query -f'%s' %s > %s",ADexpression, files[i], AllelicDeph))
}

setwd('/home/fadel/DELGEME/Target/PopGenStats/Rtsne/Input/')
#---------- Phase and Impute the missing data ------------
GTfiles=list.files('.', pattern = 'GT.txt', full.names = T)
ADfiles=list.files('.', pattern = 'AD.txt', full.names = T)

for(i in 1:length(GTfiles))
{
  cat(paste0('i= ',i),'\n')
  GTfileName=gsub('_GT.txt','', basename(GTfiles[i]))
  for (j in 1:length(ADfiles)) 
  {
    cat(paste0('j= ',j),'\n')
    ADfileName=gsub('_AD.txt','', basename(ADfiles[j]))
    if(GTfileName == ADfileName)
    {
      genotypeData = fread(GTfiles[i], header = F)
      allelicDepthData = fread(ADfiles[j], header = F)
      
      
      first4Columns= subset(genotypeData, select=c(1:4))
      genotypeData= as.data.frame(subset(genotypeData, select=-c(1:4)))
      allelicDepthData= as.data.frame(subset(allelicDepthData, select=-c(1:4)))
      
      FileName =paste0(PathForGenotypeData, GTfileName,'_Phased.txt')
      
      tic();genotypePhasing(genotypeData, allelicDepthData, first4Columns,numberOfSimulation,FileName);toc()
    }
    else next
  }
}

PhasesFiles=list.files('.', pattern = '_Phased.txt', full.names = F)

path='/home/fadel/DELGEME/Target/PopGenStats/Rtsne/Input'

for (l in 1:length(PhasesFiles)) 
{
  cat(paste0('l= ',l),'\n')
  fileName=gsub('_Phased.txt','_Imputed.txt', basename(PhasesFiles[l]))
  
  PhasedData=fread(PhasesFiles[l], header = F)
  first4Columns= subset(PhasedData, select=c(1:4))
  PhasedData=subset(PhasedData, select = -c(1:4))
  
  tic();imputeMissingGenotypes(PhasedData, first4Columns,numberOfSimulation,fileName, path);toc()
}


#------------------------------------------------
#------ This script is to compute t-SNE on SNP data
#----- It takes as input the imputed genotypes data

library(data.table)
library(caret)
library(Rtsne)
library(ggplot2)
library(gridExtra)

options(scipen = 999)
setwd('/home/fadel/DELGEME/Target/PopGenStats/Rtsne/Input')

Output='/home/fadel/DELGEME/Target/PopGenStats/Rtsne/Output'

#chroms=c("Pf3D7_01_v3","Pf3D7_02_v3","Pf3D7_03_v3","Pf3D7_04_v3","Pf3D7_07_v3", "Pf3D7_08_v3","Pf3D7_09_v3",
#         "Pf3D7_10_v3","Pf3D7_11_v3","Pf3D7_12_v3","Pf3D7_13_v3","Pf3D7_14_v3")

files=list.files('.', pattern = '_Imputed.txt', full.names = T)

SamplesID <- fread('/home/fadel/DELGEME/Target/MetaAfricanSamples.txt', header = F)

couleurs <- c("red", "blue", "black","green3","yellow","magenta","forestgreen","blue4","firebrick","cyan","antiquewhite4","darkorange","deeppink2","gray","dodgerblue","darkviolet")

Populations=c('Guinea','Ghana','Mali','Cote_DIvoire','Kenya','Gabon','Gambia','Mauritanie','Cameroon','Tanzania','Ethiopia','Madagascar','Nigeria','Malawi','DRCongo','Senegal')

for(i in 1:length(files))
{
  cat(paste0('************ Running File= ', files[i]),'***********\n')
  data_tsne=read.delim(files[i], header =F, stringsAsFactors = F, sep = "")
  data=data_tsne[,-c(1:4)]
  data=as.matrix(t(data))
  FileName=gsub('_Imputed.txt', '', basename(files[i]))
  system(sprintf('mkdir %s', paste0(Output,FileName)))
  ## Rtsne function may take some minutes to complete...      
  set.seed(123)
  
  j=20
  while (j <=60) 
  {
    cat(paste0('****************** Perplexity= ',j),'****************\n')
    tsne_model = Rtsne(data, dims=3, check_duplicates=FALSE, pca=TRUE, perplexity=j, max_iter = 100000, theta=0.5)
    
    ## getting the two dimension matrix
    d_tsne =cbind(SamplesID, as.data.frame(tsne_model$Y))
    names(d_tsne)=c('SampleID',"Population","Dim1","Dim2","Dim3")
    
    write.table(d_tsne, paste0(Output,FileName,'/',FileName,'-', j, '.txt'), col.names = T, row.names = F, quote = F)
    
    j=j+10
  }
  
}

    
#-------- Plot the 3 plots in the same graph

p1 <- ggplot(lsa3,aes( x=lsa3$Dim1, y=lsa3$Dim2, colour=Population)) + 
  geom_point()+ 
  labs(title="LSA3", x = "", y = "Dim2", color = "Populations\n")  +
  scale_color_manual(labels = Populations, values = couleurs) + 
  theme_bw() + 
  theme_light(base_size=10) +
  theme(legend.direction = "horizontal", 
        legend.position = "bottom",
        legend.box = "horizontal")

p2 <- ggplot(csp,aes( x=csp$Dim1, y=csp$Dim2, colour=Population)) + 
  geom_point()+ 
  labs(title="CSP", x = "", y = "", color = "Populations\n")  +
  scale_color_manual(labels = Populations, values = couleurs) + 
  theme_bw() + 
  theme_light(base_size=10) +
  theme(legend.direction = "horizontal", 
        legend.position = "bottom",
        legend.box = "horizontal") #  + guides(fill = FALSE, color = FALSE, linetype = FALSE, shape = FALSE)


p3 <- ggplot(glurp,aes( x=glurp$Dim1, y=glurp$Dim2, colour=Population)) + 
  geom_point()+ 
  labs(title="GLURP", x = "Dim1", y = "Dim2", color = "Populations\n")  +
  scale_color_manual(labels = Populations, values = couleurs) + 
  theme_bw() + 
  theme_light(base_size=10) +
  theme(legend.direction = "horizontal", 
        legend.position = "bottom",
        legend.box = "horizontal") 

p4 <- ggplot(celtos,aes( x=celtos$Dim1, y=celtos$Dim2, colour=Population)) + 
  geom_point()+ 
  labs(title="CelTos", x = "Dim1", y = "", color = "Populations\n")  +
  scale_color_manual(labels = Populations, values = couleurs) + 
  theme_bw() + 
  theme_light(base_size=10) +
  theme(legend.direction = "horizontal", 
        legend.position = "bottom",
        legend.box = "horizontal") 

p5 <- ggplot(trap,aes( x=trap$Dim1, y=trap$Dim2, colour=Population)) + 
  geom_point()+ 
  labs(title="TRAP", x = "Dim1", y = "Dim2", color = "Populations\n")  +
  scale_color_manual(labels = Populations, values = couleurs) + 
  theme_bw() + 
  theme_light(base_size=10) +
  theme(legend.direction = "horizontal", 
        legend.position = "bottom",
        legend.box = "horizontal")
#************ Function *************
grid_arrange_shared_legend <- function(...) {
  plots <- list(...)
  g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  grid.arrange(
    do.call(arrangeGrob, lapply(plots, function(x)
      x + theme(legend.position="none"))),
    legend,
    ncol = 1,
    heights = unit.c(unit(1, "npc") - lheight, lheight))
}



grid_arrange_shared_legend(p1, p2, p3, p4, p5)

#******************** SECOND APPROACH *****************************
grid_arrange_shared_legend1 <- function(..., nrow = 1, ncol = length(list(...)), position = c("bottom", "right")) {
  
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position = "none"))
  gl <- c(gl, nrow = nrow, ncol = ncol)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  grid.newpage()
  grid.draw(combined)
  
}

grid_arrange_shared_legend1(p1, p2, p3, p4,p5, nrow = 1, ncol = 5)
grid_arrange_shared_legend1(p1, p2, p3, p4,p5, nrow = 2, ncol = 3)

# Creating the cluster models

# Next piece of code will create the k-means and hierarchical cluster models. 
# To then assign the cluster number (1, 2 or 3) to which each input case belongs.

## keeping original data
d_tsne_1_original=d_tsne_1

## Creating k-means clustering model, and assigning the result to the data used to create the tsne
fit_cluster_kmeans=kmeans(scale(d_tsne_1), 3)
d_tsne_1_original$cl_kmeans = factor(fit_cluster_kmeans$cluster)

## Creating hierarchical cluster model, and assigning the result to the data used to create the tsne
fit_cluster_hierarchical=hclust(dist(scale(d_tsne_1)))

## setting 3 clusters as output
d_tsne_1_original$cl_hierarchical = factor(cutree(fit_cluster_hierarchical, k=3))

# Plotting the cluster models onto t-SNE output

# Now time to plot the result of each cluster model, based on the t-SNE map.

plot_cluster=function(data, var_cluster, palette)
{
  ggplot(data, aes_string(x="Dim1", y="Dim2", color=var_cluster)) +
    geom_point(size=0.25) +
    guides(colour=guide_legend(override.aes=list(size=6))) +
    xlab("") + ylab("") +
    ggtitle("") +
    theme_light(base_size=20) +
    theme(axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          legend.direction = "horizontal", 
          legend.position = "bottom",
          legend.box = "horizontal") + 
    scale_colour_brewer(palette = palette) 
}


plot_k=plot_cluster(d_tsne_1_original, "cl_kmeans", "Accent")
plot_h=plot_cluster(d_tsne_1_original, "cl_hierarchical", "Set1")

## and finally: putting the plots side by side with gridExtra lib...
library(gridExtra)
grid.arrange(plot_k, plot_h,  ncol=2)

#----------------------------------------------------
#----------- Compute iHs for each gene/pop ---------
#----- INPUT: IMPUTED FILE or VCF FILE

dir='/home/fadel/DELGEME/Target/PopGenStats/iHs/Input/AMA1'
setwd(dir)
library(data.table)
simuatedFiles = list.files(path =".", pattern ='.txt', full.names =TRUE)
nsim = length(simuatedFiles)
out='/home/fadel/DELGEME/Target/PopGenStats/iHs/output/'
Nom='AMA1'
rdsFile =paste0(out,Nom,"_simmatout.xls")
lenFile =paste0(out,Nom,"_len.xls")
files = fread(simuatedFiles[1], header = F)
NbreLigne=nrow(files)
lenof = matrix(data= NA, nrow=NbreLigne, ncol= nsim)
simmatout=matrix(data= NA, nrow=NbreLigne, ncol= nsim)


library(rehh)
library(data.table)
library(vcfR) 

setwd('/home/fadel/DELGEME/Target/PopGenStats/TestOfSelection/PopGenome-Input/PfSEA/')
options(scipen = 999)

Output='/home/fadel/DELGEME/Target/PopGenStats/iHs/Output/'

simuatedFiles = list.files(path =".", pattern ='.vcf.gz$', full.names =TRUE)
nsim = length(simuatedFiles)

rdsFile =paste0(Output, "pfsea_iHS_pvalue.xls")
lenFile =paste0(Output, "pfsea_iHS_position.xls")

lenof = matrix(data= NA, nrow=5, ncol= nsim)
simmatout=matrix(data= NA, nrow=5, ncol= nsim)
for (i in 1:length(simuatedFiles)) 
{
  cat(paste0('File = ',gsub('_GT.vcf.gz', '', basename(simuatedFiles[i]))),'\n')
  Haplo <- data2haplohh(hap_file =simuatedFiles[i] , polarize_vcf = FALSE)
  res.scan =scan_hh(Haplo, limhaplo = 2, limehh = 0, limehhs = 0, scalegap = NA, maxgap = NA, discard_integration_at_border = F, threads = 1)
  ihs=ihh2ihs(res.scan,freqbin = 1, min_maf = 0.01)
  
  simmatout[1:length(ihs$ihs$POSITION),i]=round(ihs$ihs$LOGPVALUE,4)
  lenof[1:length(ihs$ihs$POSITION),i] = ihs$ihs$POSITION
  
}

colnames(simmatout) =c("Cameroon", "Congo", "CoteIvoire", "Ethiopia","Gabon","Gambia","Ghana","Guinea","Kenya","Madagascar","Malawi","Mali","Mauritanie","Nigeria","Senegal","Tanzania")
colnames(lenof) =c("Cameroon", "Congo", "CoteIvoire", "Ethiopia","Gabon","Gambia","Ghana","Guinea","Kenya","Madagascar","Malawi","Mali","Mauritanie","Nigeria","Senegal","Tanzania")

write.table(simmatout, rdsFile, sep = '\t', col.names = T, row.names = F, quote = F)
write.table(lenof, lenFile, sep = '\t', col.names = T, row.names = F, quote = F)

#--------------- COMPUTE IBD USING isoRelate ---------
#---- SEE THE R CODE isoRelate FROM THE DELGEME PROJECT FOLDER
#-------- INPUT FORMAT: PED AND MAP FILE


#-------------- HAPLOTYPE NETWORK -----------
#----- USING PEGAS AND ADE PACKAGES ---------
#------ INPUT FORMAT: Fasta file

#----- CREATE fasta FORMAT FOR EACH VACCINE CANDIDATE 
#---------------- Haplotype Function ----------------------
HaplotypeFormat <- function(genotype, first4columns, output)
{
  genotype= as.data.frame(genotype)
  
  hapFormat=data.frame(matrix(ncol = ncol(genotype), nrow = nrow(genotype)))
  
  for (i in 1:nrow(genotype)) 
  {
    k=1
    for (j in 1:ncol(genotype)) 
    {
      if(genotype[i,j]=="0") 
      {
        hapFormat[i,k]<- first4columns[i,3]
      }
      
      else if(genotype[i,j]=="1")                   
      {
        hapFormat[i,k]<- first4columns[i,4]
      }
      
      k=k+1
    }
    
  }
  
  hapFormat=as.data.frame(hapFormat)
  write.table(hapFormat, file =output, sep = '\t' , row.names = F, col.names =F, quote = FALSE)
}

#---------------- Fasta Format Function ----------------------

createFasta <- function(genotype, sample_names, outputFasta)
{
  vecteur<- vector(mode="character", length(sample_names))
  for (i in 1: length(sample_names))
  {
    vecteur[i] = paste(">",sample_names[i])
  }
  
  transposer= t(genotype)

  fichier=file(outputFasta,open="w")
  
  for (i in 1:nrow(transposer))
  {
    cat(vecteur[i],file=fichier,sep = '\n')
    sample_seqs <- as.character(transposer[i,])
    sample_seqs=paste0(sample_seqs,collapse ='')
    if(nchar(sample_seqs)>500)
    {
      
      cat(substring(sample_seqs, 1, 500), file=fichier, sep = '\n')
      RestOfLine <- nchar(substring(sample_seqs, 501, nchar(sample_seqs)))
      if(RestOfLine >500)
      {
        cat(substring(sample_seqs, 501, 1000), file=fichier, sep = '\n')
        RestOfLine2 <- nchar(substring(sample_seqs, 1001, nchar(sample_seqs)))
        if(RestOfLine2 >500)
        {
          cat(substring(sample_seqs, 1001, 1500), file=fichier, sep = '\n')
          RestOfLine3 <- nchar(substring(sample_seqs, 1501, nchar(sample_seqs)))
          if(RestOfLine3 >500)
          {
            cat(substring(sample_seqs, 1501, 2000), file=fichier, sep = '\n')
            cat(substring(sample_seqs, 2001, nchar(sample_seqs)), file=fichier, sep = '\n')
          }
          else
          {
            cat(substring(sample_seqs, 1501, nchar(sample_seqs)), file=fichier, sep = '\n')
          }
        }
        else
        {
          cat(substring(sample_seqs, 1001, nchar(sample_seqs)), file=fichier, sep = '\n')
        }
        
        
      }
      else
      {
        cat(substring(sample_seqs, 501, nchar(sample_seqs)), file=fichier, sep = '\n')
      }
      
    }
    else
      cat(sample_seqs, file=fichier, sep = '\n')
  }
  
  close(fichier)
}

#----------------------------------------------
files <- list.files(path="/home/fadel/DELGEME/Target/PopGenStats/Rtsne/Input", pattern='*_Imputed.txt', full.names = TRUE)

dir ="/home/fadel/DELGEME/Target/PopGenStats/Haplotype_Network/FastaFiles"
system(sprintf("mkdir %s", dir))

for(t in 1:length(files))
{
  cat(paste0('Files = ', files[t],'\n'))
  Filenames=fread(files[t], header = F)

  output =paste0(dir,'/',gsub('_Imputed.txt','.hap',basename(files[t])))
  first4columns<- subset(Filenames, select=c(1:4))
  genotype<- subset(Filenames, select=-c(1:4))
  HaplotypeFormat(genotype, first4columns, output)

}

Hapfiles <- list.files(dir, pattern='*.hap', full.names = TRUE)

sample_names=scan('/home/fadel/DELGEME/Target/SampleList.txt', what = 'character')

for(t in 1:length(Hapfiles))
{
  cat(paste0('File = ', basename(Hapfiles[t]),'\n'))
  Filenames=fread(Hapfiles[t], header = F)
  outputFasta =paste0(dir,'/',gsub('.hap','.fasta',basename(Hapfiles[t])))

  createFasta(Filenames, sample_names, outputFasta)
  
  print('**********************************************')

}

install.packages(c("igraph","network","sna","ndtv"))

library(data.table)
library(ape)
library(pegas)
library(haplotypes)
library(plyr)

files=list.files('/home/fadel/DELGEME/Target/PopGenStats/Haplotype_Network/FastaFiles', pattern = 'fasta', full.names = T)

Population=fread('/home/fadel/DELGEME/Target/MetaAfricanSamples.txt', header = F)

couleurs <- c('#253494', '#800026','#e611dc','#fccde5','#00441b','#fd8d3c','#00c957','#c3135d','#ffff33','#109dee',
              '#e31a1c','#11e6e8','#238b45', '#000000','#bdbdbd','#810f7c')

Out='/home/fadel/DELGEME/Target/PopGenStats/Haplotype_Network/Network-Plots/'

for (i in 1:length(files)) 
{
  cat(paste0('File= ', basename(files[i]),'\n'))
  
  d <- ape::read.dna(files[i], format='fasta')
  
  h <- pegas::haplotype(d)
  h <- sort(h, what = "label")
  net <-pegas::haploNet(h)
  
  #ind.hap<-with(stack(setNames(attr(h, 'index'), rownames(h))), table(hap=ind, individuals=rownames(d)[values]))
  
  ind.hap2<-with(stack(setNames(attr(h, 'index'), rownames(h))), table(hap=ind, individuals=Population$V2[values]))
  
  #write.table(t(ind.hap2), '/home/fadel/DELGEME/Target/PopGenStats/Haplotype_Network/rh5-1.txt', col.names = T, row.names = T, quote = F)
  
  mydata <- as.data.frame(ind.hap2)
  
  m <- mydata[mydata$Freq != 0,]
  #n <- mydata[mydata$Freq == 1,]
  
  locations <- as.character(Population$V2)
  
  new.hap <- table(m$hap, m$individuals) 
  
  FileName=gsub('.fasta','', basename(files[i]))
  
 # pdf(paste0(Out,FileName,'.pdf'), onefile = T, paper = 'USr')
  
  #-------------------------- PLOT ------------------
  par(bg="white")
 # par(mar=c(1,1,1,1))
  
  plot(net, size =3, col.link = "black", fast = TRUE, legend=F, label=F, pie = new.hap, 
       scale.ratio = 3.5, show.mutation = 0, bg=couleurs, fg="white")
  
  title(main = FileName,  col.main= "black", cex.main=2)  #paste0("Haplotype Network Plot\n of ",
  
  #plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
  #legend('topleft', unique(locations), col=couleurs, cex=1, pch=15, ncol=1,x.intersp = 0.65, bty = 'n')
  #legend(-40, 10, unique(Population$V2), col=couleurs, pch=20, cex=0.9, ncol = 2)
  #legend(-50, 15, unique(locations), col=couleurs, cex=1, pch=15, ncol=1,x.intersp = 0.65, bty = 'n')
  
  #dev.off()
 # png(file = paste0(Out,FileName,'.png'), width = 480, height = 480, units = "px", pointsize = 12, bg = "white", type = "cairo")
  #-------------------------- PLOT ------------------
 # par(bg="white")
 # par(mar=c(1,1,1,1))
  
 # plot(net, size =3, col.link = "black", fast = TRUE, legend=F, label=F, pie = ind.hap2, 
 #      scale.ratio = 2.00, show.mutation = 0, bg=couleurs, fg="white")
  
 # title(main = paste0("Haplotype Network Plot\n of ",FileName), line = -2.5, col.main= "black", cex.main=2)
  
 # dev.off()
  
}





