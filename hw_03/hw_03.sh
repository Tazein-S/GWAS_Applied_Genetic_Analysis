#load in our packages
module load R 
module load plink/1.90b6.27

##make an alias for the directory with the data we are using today:
export DATADIR=/projectnb/bs859/data/tgen/cleaned/

#Use the file TGEN_pcs.txt to determine which PCs are associated with case status
plink --bfile $DATADIR/TGEN_cleaned --covar $DATADIR/TGEN_pcs.txt --logistic no-snp beta --ci .95 --out pc_case_status --allow-no-sex

#Which PCs are associated with case status at p<0.01?
awk 'NR==1||$8<0.01 {print $0}' pc_case_status.assoc.logistic > pc_case_status.sig.txt
head pc_case_status.sig.txt

#Perform a GWAS of Alzheimer Disease case status using log regression, adjusting for the PCs associated with case status at p<0.01 that you found in question 1
plink --bfile $DATADIR/TGEN_cleaned --covar $DATADIR/TGEN_pcs.txt --covar-name PC6,PC8 --logistic beta hide-covar --ci .95 --out logistadjPC --allow-no-sex

#How many SNPs in this GWAS have p-value < 0.0001?
awk 'NR==1||$12<0.0001{print $0}' logistadjPC.assoc.logistic > logistadjPC.sig.txt 

#Show the plink results for the most significant SNP.  
sort -k12n logistadjPC.sig.txt | head -2
awk 'NR==2 {print $2}' logistadjPC.sig.sorted.txt

#Show the 3 genotypes for the SNP  
awk '$2=="rs4693305" {print $0}' $DATADIR/TGEN_cleaned.bim

#Compute the odds ratio
plink --bfile $DATADIR/TGEN_cleaned --covar $DATADIR/TGEN_pcs.txt --covar-name PC6,PC8 --logistic hide-covar --out ORlogistadjPC --allow-no-sex

#Create a GRM from the TGEN data so that you can perform a GWAS using a logistic mixed model.  
#Use only chromosomes 1-22 (exclude X).  
#Don’t forget to prune before you make the GRM!  

#ld pruning
plink --bfile $DATADIR/TGEN_cleaned --indep-pairwise 10000kb 1 0.15 --allow-no-sex

#remove the variants that plink provided
plink --bfile $DATADIR/TGEN_cleaned --chr 1-22 --extract plink.prune.in --make-rel square --out grm --allow-no-sex

#What is the maximum relationship coefficient for MAYO_10139, other than with itself?
awk 'NR==1' grm.rel | tr '\t' '\n' | sort -gr | head -n 2

#run GMMAT -> creates N/A covariates + PC covar run
Rscript --vanilla GMMAT.R > GMMAT.log

#How many SNPs have p-value < 0.0001 in each of these two GWAS?
awk 'NR==1||$11<0.0001 {print $0}' test.glmm.score.nocov > glmm.nocov.sig.txt
wc glmm.nocov.sig.txt

awk 'NR==1||$11<0.0001 {print $0}' test.glmm.score.PC6PC8cov > glmm.PC6PC8.sig.txt
wc glmm.PC6PC8.sig.txt

#these awk commands replace the "PVAL" column header with "P", and the "POS" with "BP" . 
#This will make it match the P-value and base pair column headers from PLINK output.
##this will be useful when we use the R program qqplot.R to do a qqplot and manhattan plot 
##of the p-values
awk 'NR==1{$4="BP";$11="P"};{print $0}' test.glmm.score.nocov>glmm.score.nocov.txt
awk 'NR==1{$4="BP";$11="P"};{print $0}' test.glmm.score.PC6PC8cov>glmm.score.PC6PC8cov.txt

##Let's compare the effect estimates for one of the most associated variants:  rs11204005 
awk '(NR==1||($2=="rs4693305")){print $0}' logistadjPC.assoc.logistic
awk '(NR==1||($2=="rs4693305")){print $0}' glmm.score.nocov.txt
awk '(NR==1||($2=="rs4693305")){print $0}' glmm.score.PC6PC8cov.txt

### alternative tool for making prettier qq plots (using Umich qqplot code):
Rscript  --vanilla qq_umich_gc.R  logistadjPC.assoc.logistic log_PC6PC8 ADD 
Rscript --vanilla qq_umich_gc.R  glmm.score.nocov.txt glmm_nocov   
Rscript --vanilla qq_umich_gc.R  glmm.score.PC6PC8cov.txt glmmPC6PC8cov  

##Manhattan plots  for all of the GWAS models we tried
Rscript --vanilla gwaplot.R  logistadjPC.assoc.logistic "log with PC6 and PC8" logistadPC_manhattan
Rscript --vanilla gwaplot.R  glmm.score.nocov.txt "GLMM no cov analysis" GLMM_nocov_manhattan 
Rscript --vanilla gwaplot.R  glmm.score.PC6PC8cov.txt "GLMM PC6 and PC8 adjusted analysis" GLMM_PC68_manhattan 

