het<-read.table("Fstat.het",header=T,as.is=T)  ### reading in output from the PLINK --het command  (I called the file "Fstat", and PLINK automatically puts the .het onto output of --het)

summary(het$F)  # standard base R summary statistics -- does not provide Standard Deviation
mean(het$F)   
sd(het$F)  #standard deviation

##plot the number of non missing genotypes (het$N.NM.) versus the F-statistic
jpeg("hetplot.jpeg")
plot(het$N.NM.,het$F,xlab="N Non-Missing",ylab="F")

dev.off()

summary(lm(het$N.NM.~het$F))   ## run a simple linear regression predicting the number of non-missing genotypes from the F statistic, and print the summary:

