module load R
module load prsice

Rscript $SCC_PRSICE_BIN/PRSice.R --dir . \
--prsice $SCC_PRSICE_BIN/PRSice \
--base /projectnb/bs859/data/hippo/Meta_Whole_hippocampus.txt.gz \
--target /projectnb/bs859/data/tgen/cleaned/TGEN_cleaned \
--stat Z \
--beta \
--snp SNP \
--A1 A1 \
--A2 A2 \
--pvalue P \
--binary-target T \
--cov-file /projectnb/bs859/data/tgen/cleaned/TGEN_pcs.txt \
--cov-col PC6,PC8 \
--clump-kb 500 \
--perm 1000 \
--seed 1443 \
--out hippo_base_TGEN_target
