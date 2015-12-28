library(signal)
library(RColorBrewer)

loci=15000
mutRate<-3.1*10^(-9)
chrRegions<-diff(floor(c(0,0.5,4.5,12.5,16.5,17)*loci/17))
regRecRates<-c(0,7,1,7,0)/sum(c(0,7,1,7,0))
#chrRegRates<-rep(regRecRates,chrRegions)
#genome<-hypredGenome(num.chr=1,len.chr=c(1.0),num.snp=loci)
qtl.ids1<-c(floor(0.2*loci))
qtl.ids2<-c(floor(0.5*loci))
qtl.ids3<-c(floor(0.15*loci),floor(0.5*loci))
qtls<-list(qtl.ids1,qtl.ids2,qtl.ids3)

#varCol<-6

contFiles<-list.files(pattern="cont.*txt")
selectFiles<-list.files(pattern="select.*txt")

extractParams<-function(myFileList,fileExt=".txt") {
  "takes a list of file names and extracts parameters from filenames, removing the first 
  field"
  params<-strsplit(myFileList,fileExt)
  paramNames<-getFieldNames(params[[1]])[-1]
  params<-t(sapply(params,getFieldValues))
  colnames(params)<-paramNames
  return(data.frame(params,stringsAsFactors=FALSE))
}

getFieldNames<-function(myString){
  "extracts string part of a string of <str><num> pairs joined with _"
  x<-strsplit(myString,"[0-9_]")[[1]]
  x<-x[!sapply(x,function(x)  any(x==""))]
  return(x)
}

getFieldValues<-function(myString){
  "extracts numeric part of a string of <str><num> pairs joined with _"
  x<-strsplit(myString,"[a-zA-Z_]")[[1]]
  x<-x[!sapply(x,function(x)  any(x==""))]
  return(as.numeric(x))
}

getDistMax<-function(myData,qtlPosition,intervalStart,intervalStop){
  myMax<-apply(myData[,intervalStart:intervalStop],1,which.max)+intervalStart-1
  return(abs(myMax-qtlPosition))
}

findDist<-function(myData,qtlPositions) {
  dist2qtl<-data.frame()
  for (qtl in qtlPositions) {
    intervalStart=qtl-floor(dim(myData)[2]/15)
    intervalStop=qtl+floor(dim(myData)[2]/15)
    dist2qtl<-rbind(dist2qtl,getDistMax(myData,qtl,intervalStart,intervalStop))
  }
  row.names(dist2qtl)<-qtlPositions
  return(dist2qtl)
}

params<-extractParams(contFiles)
dataFiles<-data.frame(contFiles,selectFiles,params,stringsAsFactors=FALSE)


### ************* choose the columns to order dataFiles by (qtls and the one that varies)***********
sortFilesBy<-"G"
### *************
varNames<-c("q"="qtls","N"="popn","G"="gen","sl"="select","sq"="seqDepth","sm"="sims")
dataFiles<-dataFiles[order(dataFiles$q,dataFiles[,sortFilesBy]),]
if (!file.exists("./Plots")) {
  dir.create("./Plots")
}
today<-format(Sys.time(), "%Y%m%d")
colours<-c(brewer.pal(8,"Dark2"),brewer.pal(9,"Set1"))
pdf(file=paste0("./Plots/",today,"_SimPlots_",varNames[sortFilesBy],".pdf"), paper="a4", height=11, width=8)
par(mfrow=c(3,1))
maxDistFromQTLs<-data.frame()
for (d in 1:dim(dataFiles)[1]) {
  qtlPos<-qtls[dataFiles$q[d]+1][[1]]
  #plot single sim, raw freq
  contFreq<-as.matrix(read.table(dataFiles[d,1],header=FALSE))
  selectFreq<-as.matrix(read.table(dataFiles[d,2],header=FALSE))
  mainTitle<-paste0("Sim1 ",strsplit(dataFiles[d,1],"cont_|.txt")[[1]][2])
  plot(1:loci,contFreq[1,],type='p',xlab="SNP number",ylab="Allele1 frequency",
       main=mainTitle,ylim=c(0,1),xlim=c(0,loci),pch=16,cex=0.4,col="#00000077")
  points(1:loci, selectFreq[1,], type='p', pch=16, cex=0.4, col="#ff000077")
  lines(1:loci,sgolayfilt(contFreq[1,],p=3,n=1001),col="grey",lwd=4)
  lines(1:loci,sgolayfilt(selectFreq[1,],p=3,n=1001),col="dark red",lwd=4)
  mtext(text="qtl",at=qtlPos,cex=0.9)
  abline(v=qtlPos,lty=5)
  title(sub="q=qtl configuration, N=population size, G=generations, sl=%selected, sq=sequencing depth, sm=# of sims",cex=0.8)
  
  #plot 10 single sims, smoothed frequency difference, in different colours
  freqDiff<-selectFreq-contFreq
  smFreq<-t(apply(freqDiff,1,sgolayfilt,p=3,n=1001))
  myDist<-findDist(smFreq,qtlPos)
  mainTitle<-paste0("10 sims (smoothed)  ",strsplit(dataFiles[d,1],"cont_|.txt")[[1]][2])
  plot(1:loci,smFreq[1,],type='n',xlab="SNP number",ylab="Frequency difference", main=mainTitle,ylim=c(-1,1),xlim=c(0,loci),lwd=2)
  for (p in 1:dim(smFreq)[1]){
    lines(1:loci,smFreq[p,],type='l',lwd=2, col=colours[p])
  }
  mtext(text="qtl",at=qtlPos,cex=0.9)
  abline(v=qtlPos,lty=5)
  title(sub="q=qtl configuration, N=population size, G=generations, sl=%selected, sq=sequencing depth, sm=# of sims",cex=0.8)
  text(qtlPos,-1,labels=paste("Avr dist:", rowMeans(myDist)),cex=1)
  toAdd<-cbind(rep(dataFiles[d,sortFilesBy]),rep(dataFiles[d,"q"]),qtlPos,myDist)
  names(toAdd)<-c(varNames[sortFilesBy],"q","qtlPos", paste0("sim",1:dim(myDist)[2]))
  maxDistFromQTLs<-rbind(maxDistFromQTLs,toAdd)
    
  #plot mean and sd of all sims. no smoothing
  meanFreq<-colMeans(freqDiff)
#   seFreq<-apply(freqDiff,2,sd)/sqrt(dim(dataFiles)[1])
#   upper<-meanFreq+qt(0.975,df=dim(dataFiles)[1]-1)*seFreq
#   lower<-meanFreq-qt(0.975,df=dim(dataFiles)[1]-1)*seFreq
  sdFreq<-apply(freqDiff,2,sd)
  upper<-meanFreq+2*sdFreq
  lower<-meanFreq-2*sdFreq
  #upper<-apply(freqDiff,2,max)
  #lower<-apply(freqDiff,2,min)
  mainTitle<-paste0("10 sims (mean+-2*SD)  ",strsplit(dataFiles[d,1],"cont_|.txt")[[1]][2])
  plot(1:loci,meanFreq,type='n',xlab="SNP number",ylab="Frequency difference (mean+-2*SD)", main=mainTitle,ylim=c(-1,1),xlim=c(0,loci),lwd=1)
  polygon(c(1:loci,rev(1:loci)),c(upper,rev(lower)), col="#00000077",border=NA)
  lines(1:loci,meanFreq,lwd=1, col="blue")
  mtext(text="qtl",at=qtlPos,cex=0.9)
  abline(v=qtlPos,lty=5)
  title(sub="q=qtl configuration, N=population size, G=generations, sl=%selected, sq=sequencing depth, sm=# of sims",cex=0.8)
  text(qtlPos,-1,labels=paste("Avr dist:", rowMeans(myDist)),cex=1)
}
dev.off()

write.table(maxDistFromQTLs,file=paste0("./Plots/",today,"_maxDist_",varNames[sortFilesBy],".txt"),
            row.names=FALSE)

