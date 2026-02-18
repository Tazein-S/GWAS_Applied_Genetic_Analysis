#read in the file, skipping the first line with the eigenvalues
data <- read.table("TGEN_cleaned3_PCA.evec", skip = 1, header = FALSE)

#split the first column by the colon (:) and the FID is the first part and the ID is the second part
fid <- sub(":.*", "", data[,1])  #before the colon
iid <- sub(".*:", "", data[,1])  #after the colon

#create a new dataframe (exclude the 12th column)
covar_data <- data.frame(fid, iid, data[, 2:11])

#save to a new file
write.table(covar_data, file = "TGEN_covariates.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
