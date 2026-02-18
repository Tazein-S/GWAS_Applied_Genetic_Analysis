

library(GMMAT)
pheno<-read.table("/projectnb/bs859/data/tgen/cleaned/TGEN_cleaned.fam",header=F)
colnames(pheno)<-c("FID","IID","fa","mo","sex","case")
pcs<-read.table("/projectnb/bs859/data/tgen/cleaned/TGEN_pcs.txt",header=T,as.is=T)

# merge the PC data with the fam file (pheno) data.

pheno1<-merge(pheno,pcs,by=c("FID","IID"),all.x=TRUE)

##Read in the GRM (genetic relationship matrix)
grm<-as.matrix(read.table("grm.rel",header=F))

# Read in the grm id file:
grm.ids<-read.table("grm.rel.id",header=F)
#apply the IDs to the two dimensions of the grm.  This is how
#gmmat will know which row and column belongs to each ID

dimnames(grm)[[1]]<-dimnames(grm)[[2]]<-grm.ids[,2]

## These two commands create the Null models (no SNPs) for the score tests.  We 
## are doing two different Null models -- model1.0 has no covariates.
## model2.0 has covariates PC6 and PC8

## note that the outcome here is "case-1" -- this is because case status in plink is coded 2=affected and 1=unaffected.
## GMMAT expects case status for a binary trait to be coded 1=affected and 0=unaffected.  Using "case-1" achievese the recoding.
#The first null model has only an intercept, no covariates:
model1.0<-glmmkin(case-1~1,data=pheno1,id="IID",kins=grm,family=binomial("logit")) #no covariate
model2.0<-glmmkin(case-1~PC6+PC8,data=pheno1,id="IID",kins=grm,family=binomial("logit")) # PC6 and PC8 covariates

## these two commands perform the score test for model 1 and model 2 for all of the SNPs in the
##  plink fileset specified by the path in geno.file:

geno.file <- "/projectnb/bs859/data/tgen/cleaned/TGEN_cleaned"
glmm.score(model1.0,infile=geno.file,outfile="test.glmm.score.nocov")
glmm.score(model2.0,infile=geno.file,outfile="test.glmm.score.PC6PC8cov")
