#!/bin/bash

#SBATCH -A mdiop
#SBATCH -J haplotype
#SBATCH -p bigmem
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1   
#SBATCH --mem=5000
#SBATCH --error="haplotype.error"
#SBATCH --output="haplotype.out" 
#SBATCH --mail-type=BEGIN,FAIL,END          
#SBATCH --mail-user=mdiop@mrc.gm


module load R

chmod +x HaplotypeNetwork.R
./HaplotypeNetwork.R