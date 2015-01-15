#!/usr/bin/env Rscript

### run ./coveragePlot.r par1

### par1 = 0 > simple
### par1 = 1 > moderate
### par1 = 2 > Hard

library("RColorBrewer")


start= 0.1
step= 0.15
end = 0.4

count = (end - start) / step + 1


data <- commandArgs(T)
print(data)
#dir <- "/cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/"
wdir <- "/Users/faezeh/Projects/compare_old_new_methylFlow_simulation/coverage/"
data[1]= 0
#wdir <- getwd();
print(wdir)
##### reading files ##################
if ( data[1] == "2"){
  print("Hard Setting Plot")
  dir_new <- paste(wdir,"/hard/new/",sep="");
  dir_old <- paste(wdir,"/hard/old/",sep="");
  
  
  coverageAvg_new <- read.table(paste(dir_new,"evalAvg.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
  
  mcfCoverage_new <- read.table(paste(dir_new,"mcf.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
  
  coverageAvg_old <- read.table(paste(dir_old,"evalAvg.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
  
  mcfCoverage_old <- read.table(paste(dir_old,"mcf.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
  
  dir <- paste(wdir, "/hard/fig/", sep="")
  
  
  #dir <- "/cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/hard-Auto/"
  #dir <- "/Users/faezeh/Desktop/project/methylFlow_analyses/cpg/hard-Auto/"
  
}

if ( data[1] == "1"){
  print("Moderate Setting Plot")
  dir_new <- paste(wdir,"/moderate/new/",sep="");
  dir_old <- paste(wdir,"/moderate/old/",sep="");
  
  
  coverageAvg_new <- read.table(paste(dir_new,"evalAvg.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
  
  mcfCoverage_new <- read.table(paste(dir_new,"mcf.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
  
  coverageAvg_old <- read.table(paste(dir_old,"evalAvg.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
  
  mcfCoverage_old <- read.table(paste(dir_old,"mcf.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
  dir <- paste(wdir, "/moderate/fig/", sep="")
  
  
  
}
if ( data[1] == "0"){
  print("Simple Setting Plot")
  dir_new <- paste(wdir,"/simple/new/",sep="");
  dir_old <- paste(wdir,"/simple/old/",sep="");
  
  
  coverageAvg_new <- read.table(paste(dir_new,"evalAvg.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
  
  mcfCoverage_new <- read.table(paste(dir_new,"mcf.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
  
  coverageAvg_old <- read.table(paste(dir_old,"evalAvg.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
  
  mcfCoverage_old <- read.table(paste(dir_old,"mcf.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
  
  dir <- paste(wdir, "/simple/fig/", sep="")
}

############## different plots for differnet CpG rate #################################################################


#### plot the abundance Error for different CpG rate

print("plot abundance Error vs Coverage rate")
pdf(paste(dir,"abdVCoverage.pdf",sep=""), width =3, height =2.5, pointsize=8)
par(mar= c(5,5,2,2))

# get the range for the x and y axis
xrange_new <- range(coverageAvg_new$var)
yrange_new <- range(coverageAvg_new$abdncError)
xrange_old <- range(coverageAvg_old$var)
yrange_old <- range(coverageAvg_old$abdncError)

xrange_min <- min(xrange_new[1], xrange_old[1])
xrange_max <- max(xrange_new[2], xrange_old[2])
yrange_min <- min(yrange_new[1], yrange_old[1])
yrange_max <- max(yrange_new[2], yrange_old[2])

#ntrees <- length(unique(coverageAvg$threshold))
ntrees <- count

# set up the plot
plot(0, 0,
     pch = "",
     ylim = c(yrange_min, yrange_max),
     xlim = c(xrange_min, xrange_max),
     xlab="Coverage",
     ylab="Abundance Error",
     cex.lab= 1.2,
     #     cex.axis = 0.5
)

ltys = seq(1:ntrees)
#colors <- rainbow(ntrees)
mypalette<-brewer.pal(2*ntrees,"Dark2")
pchs = seq(0,ntrees-1)
# add lines
for (i in 1:ntrees) {
  j = step * (i-1) + start
  print(j)
  
  sel <- which(coverageAvg_new$threshold == j)
  points(coverageAvg_new$var[sel[seq(1,length(sel),by=5)]],
         coverageAvg_new$abdncError[sel[seq(1,length(sel),by=5)]],
         col = "blue",
         pch = pchs[i],
         cex = 0.5,
         type = "b",
         lty = 2,
         lwd = 1.0
  )
  
  
  
  sel <- which(readLengthAvg_old$threshold == j)
  points(coverageAvg_old$var[sel[seq(1,length(sel),by=5)]],
         coverageAvg_old$abdncError[sel[seq(1,length(sel),by=5)]],
         col = "red",
         pch = pchs[i],
         cex = 0.5,
         type = "b",
         lty = 2,
         lwd = 1.0
  )
}

#legend("topright", legend= c("cpg-loss with thr = 0.1","cpg-loss with thr = 0.25","cpg-loss with thr = 0.4","region-loss with thr = 0.1","region-loss with thr = 0.25","region-loss with thr = 0.4"),
 #      title = "thresholds",
  #     pch = rep(pchs,2),
   #    col = c("blue", "blue","blue","red","red","red"),
    #   #  cex = 0.5,
     #  lwd= 2.5)


dev.off()
# cex scale the size
#pch = 16 is circle
#lx <- seq(xrange[1] + 30 ,0.4*(xrange[2]-xrange[1])+ xrange[1], length.out=ntrees) + 5
#ly <- rep(yrange[2] + 0.1, ntrees)
#points(lx, ly,
#col = colors[1:28],
#pch = 16, cex=1,
#lwd=2)

#txt <- unique(coverageAvg$threshold)
#sel1 <- c(1, (1:count))

#text(xrange[1], yrange[2]+.1, "Legend", cex=1.3, pos=4,
#lwd=2)
#text(lx[sel1], ly[sel1], pos = 1, offset = 2, txt[sel1], cex=1.3, srt=90,
#lwd=2)




#### plot the Methyl Call  Error for different Coverage rate

print("plot abundance Error vs Coverage rate")
pdf(paste(dir,"methylVCoverage.pdf",sep=""), width =3, height =2.5, pointsize=8)
par(mar= c(5,5,2,2))

# get the range for the x and y axis
xrange_new <- range(coverageAvg_new$var)
yrange_new <- range(coverageAvg_new$methylCallError)
xrange_old <- range(coverageAvg_old$var)
yrange_old <- range(coverageAvg_old$methylCallError)

xrange_min <- min(xrange_new[1], xrange_old[1])
xrange_max <- max(xrange_new[2], xrange_old[2])
yrange_min <- min(yrange_new[1], yrange_old[1])
yrange_max <- max(yrange_new[2], yrange_old[2])

#ntrees <- length(unique(coverageAvg$threshold))
ntrees <- count

# set up the plot
plot(0, 0,
     pch = "",
     ylim = c(yrange_min, yrange_max),
     xlim = c(xrange_min, xrange_max),
     xlab="Coverage",
     ylab="Methyl Call Error",
     cex.lab= 1.2
     #     cex.axis = 1.5
)

ltys = seq(1:ntrees)
#colors <- rainbow(ntrees)
mypalette<-brewer.pal(2*ntrees,"Dark2")
pchs = seq(0,ntrees-1)
# add lines
for (i in 1:ntrees) {
  j = step * (i-1) + start
  print(j)
  
  sel <- which(coverageAvg_new$threshold == j)
  points(coverageAvg_new$var[sel[seq(1,length(sel),by=5)]],
         coverageAvg_new$methylCallError[sel[seq(1,length(sel),by=5)]],
         col = "blue",
         pch = pchs[i],
         cex = 0.5,
         type = "b",
         lty = 2,
         lwd = 1.0
  )
  
  sel <- which(coverageAvg_old$threshold == j)
  points(coverageAvg_old$var[sel[seq(1,length(sel),by=5)]],
         coverageAvg_old$methylCallError[sel[seq(1,length(sel),by=5)]],
         col = "red",
         pch = pchs[i],
         cex = 0.5,
         type = "b",
         lty = 2,
         lwd = 1.0
  )
}
#legend("topright", legend=c("region-loss","cpg-loss"),
      # col=c("red", "blue"))
#legend("topright", legend = rep(seq(start,end, by = step),2),
#legend("topright", legend= c("cpg-loss with thr = 0.1","cpg-loss with thr = 0.25","cpg-loss with thr = 0.4","region-loss with thr = 0.1","region-loss with thr = 0.25","region-loss with thr = 0.4"),
 #title = "thresholds",
  #     pch = rep(pchs,2),
   #    col = c("blue", "blue","blue","red","red","red"),
    #   #  cex = 0.5,
     #  lwd= 2.5)

dev.off()

####### plot min cost flow error for differnet Coverage ########


print("plot min cost flow error vs Coverage")
pdf(paste(dir,"mcfCoverage.pdf",sep=""), width =3, height =2.5, pointsize=8)
par(mar= c(5,5,2,2))


agg_new = aggregate(mcfCoverage_new$minCostFlow, list(Coverage = mcfCoverage_new$var), FUN =  mean)
agg_old = aggregate(mcfCoverage_old$minCostFlow, list(Coverage = mcfCoverage_old$var), FUN =  mean)

xrange_new <- range(agg_new$Coverage)
yrange_new <- range(agg_new$x)
xrange_old <- range(agg_old$Coverage)
yrange_old <- range(agg_old$x)

xrange_min <- min(xrange_new[1], xrange_old[1])
xrange_max <- max(xrange_new[2], xrange_old[2])
yrange_min <- min(yrange_new[1], yrange_old[1])
yrange_max <- max(yrange_new[2], yrange_old[2])





plot(agg_new[seq(1,length(agg_new[,1]), by=2),1],agg_new[seq(1,length(agg_new[,1]), by=2),2], xlab="Coverage", 
     ylab="Min Cost Flow Error", 
     ylim=c(yrange_min, yrange_max), 
     xlim= c(xrange_min, xrange_max), 
     cex.lab = 1.2, 
     cex.axis = 1.0, 
     pch = 1, 
     cex = 0.7,
     lty = 1,
     type = "b",
     col= "blue")


points(agg_old[seq(1,length(agg_old[,1]), by=2),1],agg_old[seq(1,length(agg_old[,1]), by=2),2],  
cex.axis = 1.0,
cex = 0.7,
pch = 1, 
lty = 1,
col= "red")
#legend("topright", legend = c("CpG-loss", "region-loss"),
 #      pch = 19,
  #     col = c("blue","red"),
  #     #  cex = 0.5,
   #    lwd= 2.5)

dev.off()




