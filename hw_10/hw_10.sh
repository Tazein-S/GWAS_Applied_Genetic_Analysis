export IGAPDIR='/projectnb/bs859/data/igap/'
export TGENDIR='/projectnb/bs859/data/tgen/cleaned/'

ls $IGAPDIR
zcat $IGAPDIR/Kunkle_etal_Stage1_results2019.txt.gz | head
ls $TGENDIR

source run_prsice.sh

cat IGAP_base_TGEN_target.summary 
cat IGAP_base_TGEN_target_r2_0.05.summary 
cat IGAP_base_TGEN_target_r2_0.15.summary
cat IGAP_base_TGEN_target_r2_0.10_kb500.summary

source run_prsiceHippo.sh
