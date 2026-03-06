##SMMAT gene-based tests for TGEN chromosome 19, topmed imputed data
##where R2>0.3
## 
library(GMMAT)
pheno<-read.table("/projectnb/bs859/data/tgen/annotated_imputed_vcfs/tgen.psam",
header=T,comment.char="")
colnames(pheno)<-c("FID","IID","sex","case")
pcs<-read.table("/projectnb/bs859/data/tgen/cleaned/TGEN_pcs.txt",
	header=T,as.is=T)

##already run:  converts vcf to gds format
#SeqArray::seqVCF2GDS("anno1.chr19.vcf.gz", "chr19.gds")
###
# merge the PC data with the fam file (pheno) data.

pheno1<-merge(pheno,pcs,by=c("FID","IID"),all.x=TRUE)

##Read in the GRM (genetic relationship matrix) -- I saved the one from
##last week to the same directory as the annotated vcfs:
grm<-as.matrix(read.table("/projectnb/bs859/data/tgen/annotated_imputed_vcfs/grm.rel",header=F))

# Read in the grm id file:
grm.ids<-read.table("/projectnb/bs859/data/tgen/annotated_imputed_vcfs/grm.rel.id",header=F)
#apply the IDs to the two dimensions of the grm.  This is how
#gmmat will know which row and column belongs to each ID

dimnames(grm)[[1]]<-dimnames(grm)[[2]]<-grm.ids[,2]

## Null model (PC covariates, no SNPs, logistic model)

model.0<-glmmkin(case-1~PC6+PC8,data=pheno1,id="IID",kins=grm,family=binomial("logit")) 

print("##check Null model coefficients for significance:")
coef<-model.0$coef
sd<-sqrt(diag(model.0$cov))
coef.pval<-pchisq((coef/sd)^2,1,lower.tail=F)
coef.table<-cbind(coef,sd,coef.pval)
coef.table

## Performs the SMMAT efficient score test for all of the groups
## of SNPs in the groups file.  Limits to variants with MAF<=0.05, and uses
## the MAF weights proposed by Wu et al.
###see documentation to understand the different choices here:
###https://www.rdocumentation.org/packages/GMMAT/versions/1.5.0/topics/SMMAT


chr19.exonic05<-SMMAT(model.0, 
     "/projectnb/bs859/data/tgen/annotated_imputed_vcfs/chr19.gds", 
     group.file="chr19.exonic.smmat.groups", 
     group.file.se = " ",
     meta.file.prefix = NULL, MAF.range = c(1e-7, 0.05),
     MAF.weights.beta = c(1, 25), miss.cutoff = 1,
     missing.method = "impute2mean", method = "davies",
     tests =c("O", "E"), rho = c(0, 0.1^2, 0.2^2, 0.3^2, 0.4^2,
     0.5^2, 0.5, 1), use.minor.allele = FALSE,
     auto.flip = FALSE, Garbage.Collection = FALSE,
     is.dosage = TRUE, ncores = 1, verbose = TRUE)

##Create a subset from the results only the gene-based tests that 
##include more than one variant
chr19.exonic05.nvargt1<-subset(chr19.exonic05,n.variants>1)

##Create a subset from the results only the gene-based tests that only have one variant
chr19.exonic05.nvar1<-subset(chr19.exonic05,n.variants==1)

##Sort the results by p-value (smallest to largest)
chr19.exonic05.nvargt1<-chr19.exonic05.nvargt1[order(chr19.exonic05.nvargt1$E.pval),]

##Add a column called "cMAC" that is the cumulative minor allele count
##for each gene.  This is computed from the number of variants in the gene
##multiplied by the mean variant allele frequency multipleid
##by 2 times the number if individuals included in the test: 

chr19.exonic05.nvargt1$cMAC<-chr19.exonic05.nvargt1$n.variants*chr19.exonic05.nvargt1$freq.mean*length(model.0$id_include)*2

print("##the number of gene tests with one variant")
print(nrow(chr19.exonic05.nvar1))

print("##print a summary of the results:")
summary(chr19.exonic05.nvargt1)

print("##number of individuals in the analysis")
print(length(model.0$id_include))

print("#print the first 5 rows (smallest 5 gene p-values)")
chr19.exonic05.nvargt1[1:5,]

print("#print a table of the number of genes with cMAC>=10")
print("#This can be used to apply a Bonferroni correction, since we")
print("#don't count gene-based tests where cMAC<10.")

table(chr19.exonic05.nvargt1$cMAC>=10)

print("##the sig level Bonferroni corrected")
sig.level<-0.05/table(chr19.exonic05.nvargt1$cMAC>=10)["TRUE"]
print(sig.level)

print("#print the genes where the p-value is < bonferroni significance level")
subset(chr19.exonic05.nvargt1,E.pval<=sig.level)

##write all the results to a comma delimited text file
write.csv(chr19.exonic05.nvargt1,"chr19.exonic05.csv",row.names=F,quote=F)

