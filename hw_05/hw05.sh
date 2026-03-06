#!/bin/bash

module load R  

export DATADIR="/projectnb/bs859/data/tgen/annotated_imputed_vcfs"

#ls $DATADIR

##take a look at the example data.  The first 33 lines of the 
##VCF file are "information"
##the 34th line gives the IDs, and the 35th is the first variant.

#zcat $DATADIR/1000G_exome_chr20_example_softFiltered.calls.hg19_anno.vcf.gz |head -n 37 |cut -f1-10

## the VCF only stores genotype data (and annotations).  A separate file stores 
##phenotype and covariate data:

#head $DATADIR/example.ped
#head example.smmat.exonic.groups

##example code for doing SMMAT gene-based analyses using the imputed genotype
##data we used last session.  Send screen output to a log file:

#MAF < 0.05
Rscript --vanilla SMMAT_tgenCHR19.MAF5.R > SMMAT_tgenCHR19.MAF5.log

#MAF < 0.01
Rscript --vanilla SMMAT_tgenCHR19.MAF1.R > SMMAT_tgenCHR19.MAF1.log

#Using MAF 0.05 and flat weights
Rscript --vanilla SMMAT_tgenCHR19.FLAT.R > SMMAT_tgenCHR19.FLAT.log

#Making a scatterplot for the default vs flat weights 
Rscript --vanilla weight_smmat.R chr19.exonic05.csv chr19.exonic05.FLAT.csv


