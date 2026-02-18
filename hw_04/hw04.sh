module load plink2/alpha3.7

#location of the imputed genotypes for TGEN chromosome 2:
export CHR2DIR=/projectnb/bs859/data/tgen/2023/chr2/

#location of the case/control files
export DATADIR=/projectnb/bs859/data/tgen/cleaned/

#How many imputed SNPs are in the chromosome 2 dosage vcf?
zcat $CHR2DIR/chr2.info.gz |wc

#How many imputed SNPs have both R2>=0.3 and MAF>=0.005?  
zcat $CHR2DIR/chr2.info.gz | awk 'NR>1 && $5 >= 0.005 && $7 >= 0.3' > filtered_chr2.txt
wc filtered_chr2.txt

##how many SNPs were in the original (pre-imputed) file
grep -v '^#' $CHR2DIR/chr2.vcf | wc -l

##Before we do an association anslysis, we will take a look
##at the R2 and MAF distribution of the SNPs and make sure that
##they are close to what we expect  
module load R
Rscript MAFinfo.R > MAFinfo.log

##Convert the vcf that comes from the imputation server to 
##plink2 "pgen" format, which is much faster to read in/out for analyes
##at the same time, filters the imputed variants so that variants with 
##Minor Allele Frequency (MAF)<0.005 are excluded, and variants with 
##imputation R2<0.3 are excluded.  VCF format does not include phenotypes,
##so this command also tells plink2 where to get the phenotype information:


##In plink2, you must specify the actual column number for the phenotype.
##in the fam file, the phenotype is the 6th column (fid iid mid fid sex pheno)

plink2  --double-id\
  --exclude-if-info "MAF<0.005" \
  --extract-if-info "R2>=0.3" \
  --vcf $CHR2DIR/chr2.dose.vcf.gz dosage=DS \
  --pheno /projectnb/bs859/data/tgen/cleaned/TGEN_cleaned.fam\
  --pheno-col-nums 6\
  --make-pgen \
  --out tgen_chr2_imputed

##In the last homework, we found that PCs 6 and 8 were associated with AD
##with p<0.01.  We will use these 2 PCs as covariates in our analyses.  
plink2 --pfile tgen_chr2_imputed \
--covar /projectnb/bs859/data/tgen/cleaned/TGEN_pcs.txt \
--covar-name PC6,PC8 --logistic hide-covar --out chr2_01pcs 

head chr2_01pcs.PHENO1.glm.logistic.hybrid

#making the QQ plots
module load R
Rscript  qqplot_pgen.R chr2_01pcs.PHENO1.glm.logistic.hybrid chr2.tgen ADD
Rscript  gwaplot_pgen.R chr2_01pcs.PHENO1.glm.logistic.hybrid "chr2.imputed" "chr2.imputed" 

##let's compare with the chr19 data that we used for the imputation:
##this is the same data as we had in the *.bed/*.bim/*.fam files we used last
##week.  I extracted chromosome 19 and did a little additional filtering, 
##and put it in the vcf format needed when uploading to the imputation server.
plink2  --double-id\
   --vcf $CHR19DIR/chr19.vcf \
  --pheno /projectnb/bs859/data/tgen/cleaned/TGEN_cleaned.fam\
  --pheno-col-nums 6\
   --covar /projectnb/bs859/data/tgen/cleaned/TGEN_pcs.txt\
   --covar-name PC6,PC8\
     --logistic hide-covar --out chr19_GENOTYPED 

##We can start by just looking at the most significant variants
sort -gk13 chr2_01pcs.PHENO1.glm.logistic.hybrid |head -n 5|cut -f 1-6,10-13

#finding the specific SNP
awk 'NR==1 || $2 == 127135234' chr2_01pcs.PHENO1.glm.logistic.hybrid

#finding the allele freq of T allele + R2 of the SNP
zcat $CHR2DIR/chr2.info.gz | awk 'NR==1 || $1 == "chr2:127135234:C:T"' | column -t

#was this SNP in the original vcf for chr 2?
grep -v '^#' $CHR2DIR/chr2.vcf | awk '$2 == 127892810' | column -t


