module load gcta/1.94.1
module load R

##Chen et al. (2021)
##Fasting glucose (FG) is a blood biomarker used in the diagnosis of T2D
head -2 /projectnb/bs859/data/FG/MAGIC1000G_FG_EUR.tsv

## Reformat FG summary stats to GCTA format
## Columns in MAGIC1000G_FG_EUR.tsv (from readme):
## variant  chromosome  position  effect_allele  other_allele  effect_allele_frequency  n  beta  se  pvalue
awk 'NR==1{print "SNP A1 A2 freq b se p N"};NR>1{print $1,$4,$5,$6,$7,$8,$9,$10}' \
    /projectnb/bs859/data/FG/MAGIC1000G_FG_EUR.tsv > FG.txt

##Type 2 Diabetes (T2D) summary statistics
##1. Xue et al 2018
head /projectnb/bs859/data/T2D/T2D_Xue_et_al_2018.txt
# reformat Xue summary stats for GSMR
awk 'NR==1{print "SNP A1 A2 freq b se p N"};NR>1{print $3,$4,$5,$6,$7,$8,$9,$10}'  /projectnb/bs859/data/T2D/T2D_Xue_et_al_2018.txt > T2D.ss.txt
sort -t ' ' -k 1,1 -u T2D.ss.txt > T2D.ss.dedup.txt
wc T2D.ss.txt T2D.ss.dedup.txt
mv T2D.ss.dedup.txt T2D.ss.txt

##write the exposure and outcome file names to files for gsmr to read:
echo "FG FG.txt" > exposure.txt
echo "T2D T2D.ss.txt" > outcome.txt

###run gsmr to estimate causal effect of FG on T2D
gcta64 --gsmr-file exposure.txt outcome.txt \
    --bfile /projectnb/bs859/data/1000G/plinkformat/1000G_EUR \
    --gsmr-direction 0 \
    --out FG-T2D-gsmr \
    --effect-plot

##check log file for errors first!
more FG-T2D-gsmr.log

#creating the plot
RScript plot.R

##SNPs removed due to HEIDI outlier procedure (pleiotropic):
cat FG-T2D-gsmr.pleio_snps

##results:
cat FG-T2D-gsmr.gsmr

##note:  to run both forward and backward in the same run, you can
##use --gsmr-direction 2
gcta64 --gsmr-file exposure.txt outcome.txt \
    --bfile /projectnb/bs859/data/1000G/plinkformat/1000G_EUR \
    --gsmr-direction 2 \
    --out FG-T2D-gsmr-bidir \
    --effect-plot

cat FG-T2D-gsmr-bidir.gsmr


