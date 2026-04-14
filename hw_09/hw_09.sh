export IGAP='/projectnb/bs859/data/igap'
export LDSCORES_DIR='/projectnb/bs859/data/ldscore_files'
export LIPIDS='/projectnb/bs859/data/lipids'

module load R
module load python2
module load ldsc
 
### LDL
munge_sumstats.py \
--sumstats $LIPIDS/jointGwasMc_LDL.txt.gz \
--snp rsid \
--a1 A1 \
--a2 A2 \
--signed-sumstats beta,0 \
--N-col N \
--merge-alleles $LDSCORES_DIR/w_hm3.snplist \
--out LDL

#grep a SNP that is empty in results in inital files

## Heritability of LDL (UKBB EUR LD scores)
## UKBB scores are one file per population
 
ldsc.py \
--h2 LDL.sumstats.gz \
--ref-ld $LDSCORES_DIR/UKBB.ALL.ldscore/UKBB.EUR.rsid \
--w-ld $LDSCORES_DIR/UKBB.ALL.ldscore/UKBB.EUR.rsid \
--out LDL_h2_UKBB


### HDL
munge_sumstats.py \
--sumstats $LIPIDS/jointGwasMc_HDL.txt.gz \
--snp rsid \
--a1 A1 \
--a2 A2 \
--signed-sumstats beta,0 \
--N-col N \
--merge-alleles $LDSCORES_DIR/w_hm3.snplist \
--out HDL

## Heritability of HDL (UKBB EUR LD scores)
 
ldsc.py \
--h2 HDL.sumstats.gz \
--ref-ld $LDSCORES_DIR/UKBB.ALL.ldscore/UKBB.EUR.rsid \
--w-ld $LDSCORES_DIR/UKBB.ALL.ldscore/UKBB.EUR.rsid \
--out HDL_h2_UKBB


## AD (Kunkle et al 2019, 21982 cases, 41944 controls)
munge_sumstats.py \
--sumstats $IGAP/Kunkle_etal_Stage1_results2019.txt.gz \
--snp MarkerName \
--N-cas 21982 \
--N-con 41944 \
--a1 Effect_allele \
--a2 Non_Effect_allele \
--signed-sumstats Beta,0 \
--merge-alleles $LDSCORES_DIR/w_hm3.snplist \
--out AD
 

## Genetic correlations among LDL, HDL, and AD (UKBB EUR)
## Gives LDL/HDL and LDL/AD correlations
ldsc.py \
--rg LDL.sumstats.gz,HDL.sumstats.gz,AD.sumstats.gz \
--ref-ld $LDSCORES_DIR/UKBB.ALL.ldscore/UKBB.EUR.rsid \
--w-ld $LDSCORES_DIR/UKBB.ALL.ldscore/UKBB.EUR.rsid \
--out LDL_HDL_rg_UKBB


## Gives HDL/AD correlation (missing from above since ldsc only computes
## pairwise correlations relative to phenotype 1)
ldsc.py \
--rg HDL.sumstats.gz,AD.sumstats.gz \
--ref-ld $LDSCORES_DIR/UKBB.ALL.ldscore/UKBB.EUR.rsid \
--w-ld $LDSCORES_DIR/UKBB.ALL.ldscore/UKBB.EUR.rsid \
--out HDL_AD_rg

 
## Heritability with --h2 alone for comparison
## LDL and HDL already done above (LDL_h2_UKBB, HDL_h2_UKBB)
## Run AD with --h2 to compare against the h2 printed inside the --rg output
 
ldsc.py \
--h2 AD.sumstats.gz \
--ref-ld $LDSCORES_DIR/UKBB.ALL.ldscore/UKBB.EUR.rsid \
--w-ld $LDSCORES_DIR/UKBB.ALL.ldscore/UKBB.EUR.rsid \
--out AD_h2_UKBB
 
 
## Cell type enrichment for AD (all 10 cell types)
## cell type group numbers (from $LDSCORES_DIR/1000G_Phase3_cell_type_groups/names):
## 1=Adrenal_Pancreas  2=Cardiovascular  3=CNS  4=Connective_Bone  5=GI
## 6=Hematopoietic     7=Kidney          8=Liver  9=Other  10=SkeletalMuscle
 
## 1: Adrenal_Pancreas
ldsc.py \
 --h2 AD.sumstats.gz \
--ref-ld-chr $LDSCORES_DIR/1000G_Phase3_cell_type_groups/cell_type_group.1.,$LDSCORES_DIR/1000G_EUR_Phase3_baseline/baseline. \
--w-ld-chr $LDSCORES_DIR/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
--overlap-annot \
--frqfile-chr $LDSCORES_DIR/1000G_Phase3_frq/1000G.EUR.QC. \
--pop-prev 0.10 \
--samp-prev 0.344 \
--print-coefficients \
--out AD_Adrenal_Pancreas

awk 'NR==2 {print $5, $6, $7}' AD_Adrenal_Pancreas.results
 
## 2: Cardiovascular
ldsc.py \
 --h2 AD.sumstats.gz \
--ref-ld-chr $LDSCORES_DIR/1000G_Phase3_cell_type_groups/cell_type_group.2.,$LDSCORES_DIR/1000G_EUR_Phase3_baseline/baseline. \
--w-ld-chr $LDSCORES_DIR/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
--overlap-annot \
--frqfile-chr $LDSCORES_DIR/1000G_Phase3_frq/1000G.EUR.QC. \
--pop-prev 0.10 \
--samp-prev 0.344 \
--print-coefficients \
--out AD_Cardiovascular

awk 'NR==2 {print $5, $6, $7}' AD_Cardiovascular.results
 
## 3: CNS (done in class)
ldsc.py \
 --h2 AD.sumstats.gz \
--ref-ld-chr $LDSCORES_DIR/1000G_Phase3_cell_type_groups/cell_type_group.3.,$LDSCORES_DIR/1000G_EUR_Phase3_baseline/baseline. \
--w-ld-chr $LDSCORES_DIR/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
--overlap-annot \
--frqfile-chr $LDSCORES_DIR/1000G_Phase3_frq/1000G.EUR.QC. \
--pop-prev 0.10 \
--samp-prev 0.344 \
--print-coefficients \
--out AD_CNS

awk 'NR==2 {print $5, $6, $7}' AD_CNS.results

 
## 4: Connective_Bone
ldsc.py \
 --h2 AD.sumstats.gz \
--ref-ld-chr $LDSCORES_DIR/1000G_Phase3_cell_type_groups/cell_type_group.4.,$LDSCORES_DIR/1000G_EUR_Phase3_baseline/baseline. \
--w-ld-chr $LDSCORES_DIR/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
--overlap-annot \
--frqfile-chr $LDSCORES_DIR/1000G_Phase3_frq/1000G.EUR.QC. \
--pop-prev 0.10 \
--samp-prev 0.344 \
--print-coefficients \
--out AD_Connective_Bone

awk 'NR==2 {print $5, $6, $7}' AD_Connective_Bone.results
 
## 5: GI
ldsc.py \
 --h2 AD.sumstats.gz \
--ref-ld-chr $LDSCORES_DIR/1000G_Phase3_cell_type_groups/cell_type_group.5.,$LDSCORES_DIR/1000G_EUR_Phase3_baseline/baseline. \
--w-ld-chr $LDSCORES_DIR/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
--overlap-annot \
--frqfile-chr $LDSCORES_DIR/1000G_Phase3_frq/1000G.EUR.QC. \
--pop-prev 0.10 \
--samp-prev 0.344 \
--print-coefficients \
--out AD_GI

awk 'NR==2 {print $5, $6, $7}' AD_GI.results

 
## 6: Hematopoietic
ldsc.py \
 --h2 AD.sumstats.gz \
--ref-ld-chr $LDSCORES_DIR/1000G_Phase3_cell_type_groups/cell_type_group.6.,$LDSCORES_DIR/1000G_EUR_Phase3_baseline/baseline. \
--w-ld-chr $LDSCORES_DIR/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
--overlap-annot \
--frqfile-chr $LDSCORES_DIR/1000G_Phase3_frq/1000G.EUR.QC. \
--pop-prev 0.10 \
--samp-prev 0.344 \
--print-coefficients \
--out AD_Hematopoietic

awk 'NR==2 {print $5, $6, $7}' AD_Hematopoietic.results
 
## 7: Kidney
ldsc.py \
 --h2 AD.sumstats.gz \
--ref-ld-chr $LDSCORES_DIR/1000G_Phase3_cell_type_groups/cell_type_group.7.,$LDSCORES_DIR/1000G_EUR_Phase3_baseline/baseline. \
--w-ld-chr $LDSCORES_DIR/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
--overlap-annot \
--frqfile-chr $LDSCORES_DIR/1000G_Phase3_frq/1000G.EUR.QC. \
--pop-prev 0.10 \
--samp-prev 0.344 \
--print-coefficients \
--out AD_Kidney

awk 'NR==2 {print $5, $6, $7}' AD_Kidney.results
 
## 8: Liver
ldsc.py \
 --h2 AD.sumstats.gz \
--ref-ld-chr $LDSCORES_DIR/1000G_Phase3_cell_type_groups/cell_type_group.8.,$LDSCORES_DIR/1000G_EUR_Phase3_baseline/baseline. \
--w-ld-chr $LDSCORES_DIR/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
--overlap-annot \
--frqfile-chr $LDSCORES_DIR/1000G_Phase3_frq/1000G.EUR.QC. \
--pop-prev 0.10 \
--samp-prev 0.344 \
--print-coefficients \
--out AD_Liver

awk 'NR==2 {print $5, $6, $7}' AD_Liver.results
 
## 9: Other
ldsc.py \
 --h2 AD.sumstats.gz \
--ref-ld-chr $LDSCORES_DIR/1000G_Phase3_cell_type_groups/cell_type_group.9.,$LDSCORES_DIR/1000G_EUR_Phase3_baseline/baseline. \
--w-ld-chr $LDSCORES_DIR/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
--overlap-annot \
--frqfile-chr $LDSCORES_DIR/1000G_Phase3_frq/1000G.EUR.QC. \
--pop-prev 0.10 \
--samp-prev 0.344 \
--print-coefficients \
--out AD_Other

awk 'NR==2 {print $5, $6, $7}' AD_Other.results
 
## 10: SkeletalMuscle
ldsc.py \
 --h2 AD.sumstats.gz \
--ref-ld-chr $LDSCORES_DIR/1000G_Phase3_cell_type_groups/cell_type_group.10.,$LDSCORES_DIR/1000G_EUR_Phase3_baseline/baseline. \
--w-ld-chr $LDSCORES_DIR/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
--overlap-annot \
--frqfile-chr $LDSCORES_DIR/1000G_Phase3_frq/1000G.EUR.QC. \
--pop-prev 0.10 \
--samp-prev 0.344 \
--print-coefficients \
--out AD_SkeletalMuscle

awk 'NR==2 {print $5, $6, $7}' AD_SkeletalMuscle.results
