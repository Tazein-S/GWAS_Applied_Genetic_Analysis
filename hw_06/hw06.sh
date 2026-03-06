module load metal
module load R

DATADIR=/projectnb/bs859/data/bmi/uganda
ls $DATADIR

zcat $DATADIR/Uganda_BMI_GWAS_MAF01.txt.gz| head -n 3

#running METAL
metal metal.txt > uganda_metal.log

#look at the meta-analysis results a bit before running a qqplot
head -n 3 METAANALYSIS1.TBL

#why are we getting a warning for these variants?
zcat $DATADIR/Uganda_BMI_GWAS_MAF01.txt.gz | awk 'NR==1 || $1=="10:100003534:G:A"'

#make a QQ plot
Rscript --vanilla qqplot_meta.R

#thinning out the data for the Manhattan plot
awk 'NR>1 {split($1,a,":"); p=$10; if(p+0 < 0.01 || NR%100==0) print a[1]"\t"a[2]"\t"p}' METAANALYSIS1.TBL > thin.TBL

#make a Manhattan plot
Rscript --vanilla manhattan_meta.R

#Finding the variants that had GWAS sig
awk 'NR==1 || $10 < 5e-8' METAANALYSIS1.TBL

#finding the variants from AADM CHR7 that showed high sig
#First, SNP 7:141551191:C:T
awk 'NR==1 || $1 == "7:141551191:C:T"' METAANALYSIS1.TBL

#Then, SNP 7:141549317:A:G
awk 'NR==1 || $1 == "7:141549317:A:G"' METAANALYSIS1.TBL

#The effect, SE, P-Val, and effect Allele Freq for 7:141549317:A:G in each study
#Uganda 
zcat $DATADIR/Uganda_BMI_GWAS_MAF01.txt.gz | awk 'NR==1 || $1=="7:141549317:A:G"' | awk '{print $1, $6, $7, $8, $9, $10}'

#DCC
zcat $DATADIR/Uganda_BMI_GWAS_MAF01.txt.gz | awk 'NR==1 || $1=="7:141549317:A:G"' | awk '{print $1, $11, $12, $13, $14, $15}'

#DDS
zcat $DATADIR/Uganda_BMI_GWAS_MAF01.txt.gz | awk 'NR==1 || $1=="7:141549317:A:G"' | awk '{print $1, $16, $17, $18, $19, $20}'

#AADM
zcat $DATADIR/Uganda_BMI_GWAS_MAF01.txt.gz | awk 'NR==1 || $1=="7:141549317:A:G"' | awk '{print $1, $21, $22, $23, $24, $25}'










