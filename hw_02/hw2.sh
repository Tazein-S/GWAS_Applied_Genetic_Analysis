#set the DATADIR
DATADIR=/projectnb/bs859/data/tgen/cleaned

#let's check whats in this directory
ls -lat $DATADIR

##load R, plink, and eigensoft 
##smartpca is a program in the eigensoft module
module load R
module load plink/1.90b6.27
module load eigensoft  

##remove variants with minor allele frequency <0.02, and with genotype missingness >0.02. OUTPUT: TGEN_cleaned1
plink --bfile $DATADIR/TGEN_cleaned --maf 0.02 --geno 0.02 --make-bed --out TGEN_cleaned1 --allow-no-sex

##ld prune the variants using the parameters:  --indep-pairwise 10000kb  1 0.2. OUTPUT: TGEN_cleaned2
plink --bfile TGEN_cleaned1 --indep-pairwise 10000kb 1 0.2 --out TGEN_cleaned2 --allow-no-sex

##create a new data set that keeps only the pruned-in variants and has only variants from chromosomes 1-22
plink --bfile TGEN_cleaned1 --extract TGEN_cleaned2.prune.in --chr 1-22 --make-bed --out TGEN_cleaned3 --allow-no-sex

#look at the parameter file for smartpca
cat TGEN_cleaned3.par

# The output files will be "test.evec" (the eigenvectors, or as we call them, 
# the PCs), and in test.out, there will be a lot of useful information that 
# would otherwise be put on the screen while smartpca is running.
# smartpca tests the association of the PCs with case status (or whatever
# phenotype is in the input *.fam file in the 6th column)

#run smartpca  using the parameter file test.par
smartpca -p TGEN_cleaned3.par > TGEN_cleaned3_PCA.out

#Plot PCs by case status
Rscript --vanilla plotPCs.R  TGEN_cleaned3_PCA.evec 1 2 10 

###second filtering and PCA
##remove variants with minor allele frequency <0.02, and with genotype missingness >0.01. OUTPUT: TGEN_cleaned4
plink --bfile $DATADIR/TGEN_cleaned --maf 0.02 --geno 0.01 --make-bed --out TGEN_cleaned4 --allow-no-sex

##ld prune the variants using the parameters:  --indep-pairwise 10000kb  1 0.15. OUTPUT: TGEN_cleaned5
plink --bfile TGEN_cleaned4 --indep-pairwise 10000kb 1 0.15 --out TGEN_cleaned5 --allow-no-sex

##create a new data set that keeps only the pruned-in variants and has only variants from chromosomes 1-22
plink --bfile TGEN_cleaned4 --extract TGEN_cleaned5.prune.in --chr 1-22 --make-bed --out TGEN_cleaned6 --allow-no-sex

#look at the parameter file for smartpca
cat TGEN_cleaned6.par

#run smartpca  using the parameter file test.par
smartpca -p TGEN_cleaned6.par > TGEN_cleaned6_PCA.out

#Plot PCs by case status
Rscript --vanilla plotPCs.R  TGEN_cleaned6_PCA.evec 1 2 10 

## Getting the ANOVA stats
# For TGEN_cleaned3 - only p-values < 0.005
awk '/eigenvector_._Control_Case_/ && $2 < 0.005 {print $1, $2}' TGEN_cleaned3_PCA.out

# For TGEN_cleaned6 - only p-values < 0.005
awk '/eigenvector_._Control_Case_/ && $2 < 0.005 {print $1, $2}' TGEN_cleaned6_PCA.out

#For the outliers we found from the prev class, get their positions based on PCA values
echo "TGEN_cleaned3_PCA.evec" > outliers.txt
awk '($1=="WGACON_66:WGACON_66" || $1=="WGAAD_270:WGAAD_270") {print $0}' TGEN_cleaned3_PCA.evec >> outliers.txt
echo "" >> outliers.txt
echo "TGEN_cleaned6_PCA.evec" >> outliers.txt
awk '($1=="WGACON_66:WGACON_66" || $1=="WGAAD_270:WGAAD_270") {print $0}' TGEN_cleaned6_PCA.evec >> outliers.txt

#finding out where MAYO_8170 is on the PC1 vs PC2 plots
echo "TGEN_cleaned3_PCA.evec" > MAYO_8170.txt
awk '($1=="MAYO_8170:MAYO_8170") {print $0}' TGEN_cleaned3_PCA.evec >> MAYO_8170.txt
echo "" >> MAYO_8170.txt
echo "TGEN_cleaned6_PCA.evec" >> MAYO_8170.txt
awk '($1=="MAYO_8170:MAYO_8170") {print $0}' TGEN_cleaned3_PCA.evec >> MAYO_8170.txt

#creating the new plink formatted covariate file using TGEN_cleaned3_PCA.evec
#make a log file just in case any errors appear
Rscript --vanilla covariate.R > covariate.log
