options(scipen=TRUE)
varNames<-c("N"="popn","G"="gen","sl"="select","sq"="seqDepth","sm"="sims")
Xtext<-c("N"="Population size", "G"="Number of generations", "sl"="% selected", 
         "sq"="Sequencing depth", "sm"="Number of simulations")
var1<-"G"
expDate<-"20151209"
dist2qtl<-read.table(file=paste0("./Plots/",expDate,"_maxDist_",varNames[var1],".txt"),header=TRUE)
distCols<-grep("sim*",names(dist2qtl))
chrRegions<-diff(floor(c(0,0.5,4.5,12.5,16.5,17)*loci/17))
dist2qtl<-data.frame(dist2qtl[,1:3],"recRate"=dist2qtl$qtlPos<chrRegions[2],dist2qtl[4:13])
dist2qtl$recRate<-factor(dist2qtl$recRate)
levels(dist2qtl$recRate)<-c("low","high")
###
library(reshape)
dist2QTL<-melt(dist2qtl,id=names(dist2qtl)[1:4])
colours<-c("#0000ffdd","#22ff22dd")
depVar<-which(names(dist2QTL)==varNames[var1])

pdf(file=paste0("./Plots/",today,"_distMax.pdf"), paper="a4", height=11, width=8)
par(mfrow=c(3,2))
qtlToPlot=dist2QTL$q==0

boxplot(dist2QTL[qtlToPlot,"value"]~dist2QTL[qtlToPlot,depVar],ylab="Distance of peak max to QTL",
        xlab=Xtext[var1],outpch = NA,col="light grey")
title(main=paste0("QTL at ",dist2QTL[qtlToPlot,"qtlPos"][1], " (",
                  dist2QTL[qtlToPlot,"recRate"][1], " recombination rate)"))
stripchart(dist2QTL[qtlToPlot,"value"] ~ dist2QTL[qtlToPlot,depVar], vertical = TRUE, 
           method = "jitter", pch = 16, col = colours[1], add = TRUE)

qtlToPlot=dist2QTL$q==1

boxplot(dist2QTL[qtlToPlot,"value"]~dist2QTL[qtlToPlot,depVar],ylab="Distance of peak max to QTL",
        xlab=Xtext[var1],outpch = NA,col="light grey")
title(main=paste0("QTL at ",dist2QTL[qtlToPlot,"qtlPos"][1], " (",
                  dist2QTL[qtlToPlot,"recRate"][1], " recombination rate)"))
stripchart(dist2QTL[qtlToPlot,"value"] ~ dist2QTL[qtlToPlot,depVar], vertical = TRUE, 
           method = "jitter", pch = 16, col = colours[2], add = TRUE)


qtlToPlot=dist2QTL$q==2 & dist2QTL$qtlPos==2250

boxplot(dist2QTL[qtlToPlot,"value"]~dist2QTL[qtlToPlot,depVar],ylab="Distance of peak max to QTL",
        xlab=Xtext[var1],outpch = NA,col="light grey")
title(main=paste0("QTL at ",dist2QTL[qtlToPlot,"qtlPos"][1], " (",
                  dist2QTL[qtlToPlot,"recRate"][1], " recombination rate)"))
stripchart(dist2QTL[qtlToPlot,"value"] ~ dist2QTL[qtlToPlot,depVar], vertical = TRUE, 
           method = "jitter", pch = 16, col = colours[1], add = TRUE)

qtlToPlot=dist2QTL$q==2 & dist2QTL$qtlPos==7500

boxplot(dist2QTL[qtlToPlot,"value"]~dist2QTL[qtlToPlot,depVar],ylab="Distance of peak max to QTL",
        xlab=Xtext[var1],outpch = NA,col="light grey")
title(main=paste0("QTL at ",dist2QTL[qtlToPlot,"qtlPos"][1], " (",
                  dist2QTL[qtlToPlot,"recRate"][1], " recombination rate)"))
stripchart(dist2QTL[qtlToPlot,"value"] ~ dist2QTL[qtlToPlot,depVar], vertical = TRUE, 
           method = "jitter", pch = 16, col = colours[2], add = TRUE)


qtlToPlot=dist2QTL$recRate=="high"

boxplot(dist2QTL[qtlToPlot,"value"]~dist2QTL[qtlToPlot,depVar],ylab="Distance of peak max to QTL",
        xlab=Xtext[var1],outpch = NA,col="light grey")
title(main=paste0("QTL in ",dist2QTL[qtlToPlot,"recRate"][1], " recombination rate region"))
stripchart(dist2QTL[qtlToPlot,"value"] ~ dist2QTL[qtlToPlot,depVar], vertical = TRUE, 
           method = "jitter", pch = 16, col = colours[1], add = TRUE)

qtlToPlot=dist2QTL$recRate=="low"

boxplot(dist2QTL[qtlToPlot,"value"]~dist2QTL[qtlToPlot,depVar],ylab="Distance of peak max to QTL",
        xlab=Xtext[var1],outpch = NA,col="light grey")
title(main=paste0("QTL in ",dist2QTL[qtlToPlot,"recRate"][1], " recombination rate region"))
stripchart(dist2QTL[qtlToPlot,"value"] ~ dist2QTL[qtlToPlot,depVar], vertical = TRUE, 
           method = "jitter", pch = 16, col = colours[2], add = TRUE)  
dev.off()
