
DATA2='/projectnb/bs859/data/meta/downloads'

# Reformat for GCTA (if not already done)
zcat $DATA2/UKB.v2.albuminuria.n382500.sorted.tsv.gz | awk '{print $2,$5,$6,$7,$8,$9,$10,$11}' > togcta.txt

# Load GCTA
module load gcta

# Check available files
ls /projectnb/bs859/data/1000G/plinkformat/

# Conditional analysis on chr 12 using 1000G EUR LD reference
gcta64 --bfile /projectnb/bs859/data/1000G/plinkformat/1000G_EUR \
  --cojo-file togcta.txt --cojo-slct --chr 12 \
  --out chr12 > chr12.log

# Check results
# See the selected independent signal(s)
cat chr12.jma.cojo

# See the frequency mismatches (for concerns)
cat chr12.freq.badsnps

# Count bad snps
wc chr12.badsnps

#now run it on African samples
gcta64 --bfile /projectnb/bs859/data/1000G/plinkformat/1000G_AFR \
  --cojo-file togcta.txt --cojo-slct --chr 12 \
  --out chr12.AFR > chr12.AFR.log

# Check results when done
# See the selected independent signal(s)
cat chr12.AFR.jma.cojo

# See the frequency mismatches (for concerns)
wc chr12.AFR.freq.badsnps

# Count bad snps
wc chr12.AFR.badsnps

#getting the top variants in the region for the VEP
zcat $DATA2/UKB.v2.albuminuria.n382500.sorted.tsv.gz | \
  awk '($3==12 && $4>=69850008 && $4<=70008337 && $10<5e-7) {print $2}' \
  > toVEP_chr12.txt

wc -l toVEP_chr12.txt
cat toVEP_chr12.txt

