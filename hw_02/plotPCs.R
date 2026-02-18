##Rscript --vanilla plotPCs.R filename x  y NPC 
##Assumes smartpca file "filename" with Number of PCs+2 columns 
## x and y are PC (numbers)
## NPC is the number of PCs in the file 
##

args<-commandArgs(trailingOnly=TRUE)
print(args)
infile<-args[1]
x<-as.numeric(args[2])
y<-as.numeric(args[3])
N<-as.numeric(args[4])
#infile<-"test.evec"
#N<-10
#x<-1
#y<-2
yy<-read.table(infile,header=F,skip=1,col.names=c("ID",paste("PC",1:N,sep=""),"CASE"),as.is=T)

bitmap(paste(c(infile,".PC.",x,".",y,".jpeg"),collapse=""))
plot(yy[,x+1],yy[,y+1],col=as.numeric(as.factor(yy[,N+2])),xlab=paste("PC",x,sep=""),ylab=paste("PC",y,sep=""))
legend("bottom",legend=unique(as.factor(yy[,N+2])),col=as.numeric(unique(as.factor(yy[,N+2]))),horiz=TRUE,pch=1,bty="n")
dev.off()

