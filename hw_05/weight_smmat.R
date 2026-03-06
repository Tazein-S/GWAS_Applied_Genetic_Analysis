args <- commandArgs(trailingOnly = TRUE)
print(args)
myfile_default <- args[1]
myfile_flat    <- args[2]
dat_default <- read.csv(myfile_default, as.is=T, header=T, comment.char="")
dat_flat    <- read.csv(myfile_flat,    as.is=T, header=T, comment.char="")
jpeg("weights_comparison.jpeg")
plot(-log10(dat_default$E.pval), -log10(dat_flat$E.pval),
     xlab="-log10(P) Default Weights (1,25)",
     ylab="-log10(P) Flat Weights (1,1)",
     main="SMMAT-E: Flat vs Default MAF Weights")
abline(0,1)
dev.off()
