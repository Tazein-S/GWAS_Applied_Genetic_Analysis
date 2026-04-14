module load R
module load prsice

# Q1 run 1: default (r2=0.1, kb=250)
Rscript $SCC_PRSICE_BIN/PRSice.R --dir . \
--prsice $SCC_PRSICE_BIN/PRSice \
--base /projectnb/bs859/data/igap/Kunkle_etal_Stage1_results2019.txt.gz \
--target /projectnb/bs859/data/tgen/cleaned/TGEN_cleaned \
--stat Beta \
--snp MarkerName \
--A1 Effect_allele \
--A2 Non_Effect_allele \
--pvalue Pvalue \
--binary-target T \
--cov-file /projectnb/bs859/data/tgen/cleaned/TGEN_pcs.txt \
--cov-col PC6,PC8 \
--clump-r2 0.10 \
--extract IGAP_base_TGEN_target.valid \
--out IGAP_base_TGEN_target

# Q1 run 2: r2=0.05
Rscript $SCC_PRSICE_BIN/PRSice.R --dir . \
--prsice $SCC_PRSICE_BIN/PRSice \
--base /projectnb/bs859/data/igap/Kunkle_etal_Stage1_results2019.txt.gz \
--target /projectnb/bs859/data/tgen/cleaned/TGEN_cleaned \
--stat Beta \
--snp MarkerName \
--A1 Effect_allele \
--A2 Non_Effect_allele \
--pvalue Pvalue \
--binary-target T \
--cov-file /projectnb/bs859/data/tgen/cleaned/TGEN_pcs.txt \
--cov-col PC6,PC8 \
--clump-r2 0.05 \
--extract IGAP_base_TGEN_target_r2_0.05.valid \
--out IGAP_base_TGEN_target_r2_0.05

# Q1 run 3: r2=0.15
Rscript $SCC_PRSICE_BIN/PRSice.R --dir . \
--prsice $SCC_PRSICE_BIN/PRSice \
--base /projectnb/bs859/data/igap/Kunkle_etal_Stage1_results2019.txt.gz \
--target /projectnb/bs859/data/tgen/cleaned/TGEN_cleaned \
--stat Beta \
--snp MarkerName \
--A1 Effect_allele \
--A2 Non_Effect_allele \
--pvalue Pvalue \
--binary-target T \
--cov-file /projectnb/bs859/data/tgen/cleaned/TGEN_pcs.txt \
--cov-col PC6,PC8 \
--clump-r2 0.15 \
--extract IGAP_base_TGEN_target_r2_0.15.valid \
--out IGAP_base_TGEN_target_r2_0.15

# Q1 run 4: r2=0.1, kb=500
Rscript $SCC_PRSICE_BIN/PRSice.R --dir . \
--prsice $SCC_PRSICE_BIN/PRSice \
--base /projectnb/bs859/data/igap/Kunkle_etal_Stage1_results2019.txt.gz \
--target /projectnb/bs859/data/tgen/cleaned/TGEN_cleaned \
--stat Beta \
--snp MarkerName \
--A1 Effect_allele \
--A2 Non_Effect_allele \
--pvalue Pvalue \
--binary-target T \
--cov-file /projectnb/bs859/data/tgen/cleaned/TGEN_pcs.txt \
--cov-col PC6,PC8 \
--clump-r2 0.10 \
--clump-kb 500 \
--extract IGAP_base_TGEN_target.valid \
--out IGAP_base_TGEN_target_r2_0.10_kb500
