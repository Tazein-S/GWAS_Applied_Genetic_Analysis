#set a shell variable to path to class data. Now we can use $DATADIR to refer to the longer directory name: 
export DATADIR=/projectnb/bs859/data/tgen  

#load plink software module:
module load plink/1.90b6.27 

#what's in the plink_tutorial directory we refer to as $DATADIR ?
ls -lat $DATADIR
   
# number of individuals
wc $DATADIR/TGEN.fam  

# number of variants 
wc $DATADIR/TGEN.bim 

##run some plink analyses. 
##compute the allele frequencies, and save them to a file freq1.frq
plink --bfile $DATADIR/TGEN --freq --allow-no-sex --out freq1   

## compute the genotype frequencies
plink --bfile $DATADIR/TGEN --freqx --allow-no-sex --out freqx  
    
#how many variants have no minor (A1) allele homozygotes
awk 'NR>1 && $5==0' freqx.frqx | wc -l

#missing sample/variant information
plink --bfile $DATADIR/TGEN --missing --allow-no-sex --out miss1 
   
#hwe test
plink --bfile $DATADIR/TGEN --hardy --allow-no-sex --out hwe1

## now filter in ONE command
plink --bfile $DATADIR/TGEN --maf 0.01 --geno 0.05 --mind 0.05 --hwe 1e-6 --allow-no-sex --make-bed --out wgas2

# first remove SNPs with MAF<0.01 and missing rate>0.05 and HWE p<0.000001 all in one step
plink --bfile $DATADIR/TGEN --maf 0.01 --geno 0.05 --hwe 1e-6 --allow-no-sex --make-bed --out wgas3_step1

# THEN you remove individuals with >5% missing genotypes.
plink --bfile wgas3_step1 --mind 0.05 --allow-no-sex --make-bed --out wgas3_step2

#which individuals (provide the IDs) were excluded in 1) but included in 2)?
module load R   ##need to load R module before using R
Rscript filter_diff.R>filter_diff.log

#ld prune
plink --bfile wgas2 --indep-pairwise 10000 kb  1  0.2 --out try1

#IBD estimation with the prune-in SNP subset:
plink --bfile wgas2 --chr 1-22 --extract try1.prune.in --genome --out ibd
wc ibd.genome
head ibd.genome

#how many pairs share 25% IBD
awk '$10>0.25{print $0}' ibd.genome > ibd.gt25.genome
wc ibd.gt25.genome

#how many unique individuals share more than at least 25% of their alleles with another person 
# sort -u makes it unique labels only
awk 'NR>1 {print $2"\n"$4}' ibd.gt25.genome | sort -u > ibd.gt25.unique.txt

#how many samples appear to be genetic duplicates?
awk '$10>0.90{print $0}' ibd.genome > ibd.duplicates.genome
wc ibd.duplicates.genome

#Summarize the non-genetic-duplicate relationships you see in the file that are 2nd degree relative or closer (pairs with PI_HAT<0.95 and PI_HAT>0.25)  
awk '$10>0.25 && $10<0.95 {print $0}' ibd.genome | sort > ibd.related.genome

#who has the most relations with other samples in our new set?
awk 'NR>1 {print $1}' ibd.gt25.genome | sort | uniq -c > most_relations.txt
sort -nr most_relations.txt

#remove the samples with IBD higher than 0.25
plink --bfile wgas2 --extract try1.prune.in --remove ibd.gt25.genome --allow-no-sex --make-bed --out wgas2_unrelated

#compute the heterozyogote deficit F statistics:
plink --bfile wgas2_unrelated --chr 1-22 --extract try1.prune.in --het --allow-no-sex --out Fstat

#running the basic stats on the Fstat.het 
Rscript Fstat.R>het_stats.log

#individuals with low F
awk 'NR>1 && $6 < -0.2 {print $0}' Fstat.het > low_F.txt

#individuals with high F
awk 'NR>1 && $6 > 0.1 {print $0}' Fstat.het > high_F.txt

