infile<-gzfile("/projectnb/bs859/data/tgen/2023/chr2/chr2.info.gz")
info<-read.table(infile, header=T, as.is=T)
info[1:10,]
##MAF is minor allele frequency. Rsq is the imputation R2

##There are a lot of variants with MAF=0 (e.g., the TGEN participants
##all have the same homozygous genotype):
table(info$MAF==0)
##we will omit the MAF==0 variants, as they are not useful for association.
info<-subset(info, MAF>0)

dim(info)

##we are not omitting R2<0.3 or MAF<0.005 yet, as we first want to get an
##overall view of the imputation and how these parameters are related

## Produce summary of MAF and INFO 

summary(info$MAF)
summary(info$Rsq)

##Extract the SNP location from the snp name: 
info$BP<-sapply(info$SNP,function(x)as.numeric(strsplit(x,split=":")[[1]][2]))


info$binned.MAF<-cut(info$MAF,breaks=c(0,0.005,0.01,0.1,0.25,0.5),include.lowest=T)
info$binned.Rsq<-cut(info$Rsq,breaks=c(0,0.1,0.3,0.5,0.8,1),include.lowest=T)
info$binned.BP<-cut(info$BP,breaks=50)

table(info$binned.MAF,info$binned.Rsq,useNA="always")

par(mfrow=c(1,1))

jpeg("boxplot_Rsq_MAF.jpg")
boxplot(info$Rsq~info$binned.MAF,main="Imputed SNP distribution by Rsq and MAF",xlab="MAF")
dev.off()

table(info$binned.BP,info$binned.MAF)
jpeg("BPvsRsqbyMAF.jpg",quality=100,width=4000,height=3000)
par(mfrow=c(2,2))
for(i in unique(info$binned.MAF)[-c(1)]){
  x<-subset(info,binned.MAF==i)
  boxplot(x$Rsq~x$binned.BP,main=paste("MAF=",i,sep=""),xlab="BP",ylab="R2")}
dev.off()


