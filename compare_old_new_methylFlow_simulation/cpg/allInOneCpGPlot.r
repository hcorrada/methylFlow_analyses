#!/usr/bin/env Rscript

### run ./allInOneCpGPlot.r


library("RColorBrewer")


start= 0.1
step= 0.15
end = 0.4

count = (end - start) / step + 1



#dir <- "/cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/"
wdir <- "/Users/faezeh/Projects/methylFlow/exps/compare_old_new_methylFlow_simulation/cpg/"
#data[1]= 1
#wdir <- getwd();
print(wdir)
##### reading files ##################

  dir_hard <- paste(wdir,"/hard/new/",sep="");

  
  cpgAvg_hard <- read.table(paste(dir_hard,"evalAvg.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
  
  mcfCpG_hard <- read.table(paste(dir_hard,"mcf.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)




  dir_moderate <- paste(wdir,"/moderate/new/",sep="");
  
  cpgAvg_moderate <- read.table(paste(dir_moderate,"evalAvg.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
  
  mcfCpG_moderate <- read.table(paste(dir_moderate,"mcf.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
  

  

  dir_simple <- paste(wdir,"/simple/new/",sep="");

  cpgAvg_simple <- read.table(paste(dir_simple,"evalAvg.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
  
  mcfCpG_simple <- read.table(paste(dir_simple,"mcf.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)


#dir <- paste(wdir, "/allInOneFig/", sep="")
dir <- file.path(wdir, "allInOne")

############## different plots for differnet CpG rate #################################################################


#### plot the abundance Error for different CpG rate

print("plot abundance Error vs CpG rate")
pdf(paste(dir,"abdVCpG.pdf",sep=""), width =3, height =2.5, pointsize=8)
par(mar= c(5,5,2,2))

# get the range for the x and y axis
xrange_hard <- range(cpgAvg_hard$var)
xrange_moderate <- range(cpgAvg_moderate$var)
xrange_simple <- range(cpgAvg_simple$var)

yrange_hard <- range(cpgAvg_hard$abdncError)
yrange_moderate <- range(cpgAvg_moderate$abdncError)
yrange_simple <- range(cpgAvg_simple$abdncError)

xrange_min <- min(xrange_hard[1], xrange_moderate[1], xrange_simple[1])
xrange_max <- max(xrange_hard[2], xrange_moderate[2], xrange_simple[2])
yrange_min <- min(yrange_hard[1], yrange_moderate[1], yrange_simple[1])
yrange_max <- max(yrange_hard[2], yrange_moderate[2], yrange_simple[2])

#ntrees <- length(unique(cpgAvg$threshold))
ntrees <- count

# set up the plot
plot(0, 0,
     pch = "",
     ylim = c(yrange_min, yrange_max),
     xlim = c(xrange_min, xrange_max),
     xlab="CpG",
     ylab= "",
#    ylab="Abundance Error",
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
  
  sel <- which(cpgAvg_hard$threshold == j)
  points(cpgAvg_hard$var[sel[seq(1,length(sel),by=5)]],
         cpgAvg_hard$abdncError[sel[seq(1,length(sel),by=5)]],
         col = "red",
         pch = pchs[i],
         cex = 0.5,
         type = "b",
         lty = 3,
         lwd = 1.0
  )
  
  
  
  sel <- which(cpgAvg_moderate$threshold == j)
  points(cpgAvg_moderate$var[sel[seq(1,length(sel),by=5)]],
         cpgAvg_moderate$abdncError[sel[seq(1,length(sel),by=5)]],
         col = "green",
         pch = pchs[i],
         cex = 0.5,
         type = "b",
         lty = 3,
         lwd = 1.0
  )
  
  sel <- which(cpgAvg_simple$threshold == j)
  points(cpgAvg_simple$var[sel[seq(1,length(sel),by=5)]],
            cpgAvg_simple$abdncError[sel[seq(1,length(sel),by=5)]],
            col = "blue",
            pch = pchs[i],
            cex = 0.5,
            type = "b",
            lty = 3,
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

#txt <- unique(cpgAvg$threshold)
#sel1 <- c(1, (1:count))

#text(xrange[1], yrange[2]+.1, "Legend", cex=1.3, pos=4,
#lwd=2)
#text(lx[sel1], ly[sel1], pos = 1, offset = 2, txt[sel1], cex=1.3, srt=90,
#lwd=2)




#### plot the Methyl Call  Error for different CpG rate

print("plot abundance Error vs CpG rate")
pdf(paste(dir,"methylVCpG.pdf",sep=""), width =3, height =2.5, pointsize=8)
par(mar= c(5,5,2,2))

# get the range for the x and y axis
xrange_hard <- range(cpgAvg_hard$var)
xrange_moderate <- range(cpgAvg_moderate$var)
xrange_simple <- range(cpgAvg_simple$var)

yrange_hard <- range(cpgAvg_hard$methylCallError)
yrange_moderate <- range(cpgAvg_moderate$methylCallError)
yrange_simple <- range(cpgAvg_simple$methylCallError)

xrange_min <- min(xrange_hard[1], xrange_moderate[1], xrange_simple[1])
xrange_max <- max(xrange_hard[2], xrange_moderate[2], xrange_simple[2])
yrange_min <- min(yrange_hard[1], yrange_moderate[1], yrange_simple[1])
yrange_max <- max(yrange_hard[2], yrange_moderate[2], yrange_simple[2])

#ntrees <- length(unique(cpgAvg$threshold))
ntrees <- count

# set up the plot
plot(0, 0,
     pch = "",
     ylim = c(yrange_min, yrange_max),
     xlim = c(xrange_min, xrange_max),
     xlab="CpG",
     ylab= "",
#  ylab="Methyl Call Error",
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
  
  sel <- which(cpgAvg_hard$threshold == j)
  points(cpgAvg_hard$var[sel[seq(1,length(sel),by=5)]],
         cpgAvg_hard$methylCallError[sel[seq(1,length(sel),by=5)]],
         col = "red",
         pch = pchs[i],
         cex = 0.5,
         type = "b",
         lty = 3,
         lwd = 1.0
  )
  
  sel <- which(cpgAvg_moderate$threshold == j)
  points(cpgAvg_moderate$var[sel[seq(1,length(sel),by=5)]],
         cpgAvg_moderate$methylCallError[sel[seq(1,length(sel),by=5)]],
         col = "green",
         pch = pchs[i],
         cex = 0.5,
         type = "b",
         lty = 3,
         lwd = 1.0
  )
  
  sel <- which(cpgAvg_simple$threshold == j)
  points(cpgAvg_simple$var[sel[seq(1,length(sel),by=5)]],
        cpgAvg_simple$methylCallError[sel[seq(1,length(sel),by=5)]],
        col = "blue",
        pch = pchs[i],
        cex = 0.5,
        type = "b",
        lty = 3,
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

####### plot min cost flow error for differnet CpG ########


print("plot min cost flow error vs CpG")
pdf(paste(dir,"mcfCpG.pdf",sep=""), width =3, height =2.5, pointsize=8)
par(mar= c(5,5,2,2))


agg_hard = aggregate(mcfCpG_hard$minCostFlow, list(CpG = mcfCpG_hard$var), FUN =  mean)
agg_moderate = aggregate(mcfCpG_moderate$minCostFlow, list(CpG = mcfCpG_moderate$var), FUN =  mean)
agg_simple = aggregate(mcfCpG_simple$minCostFlow, list(CpG = mcfCpG_simple$var), FUN =  mean)


xrange_hard <- range(agg_hard$CpG)
xrange_moderate <- range(agg_moderate$CpG)
xrange_simple <- range(agg_simple$CpG)

yrange_hard <- range(agg_hard$x)
yrange_moderate <- range(agg_moderate$x)
yrange_simple <- range(agg_simple$x)

xrange_min <- min(xrange_hard[1], xrange_moderate[1], xrange_simple[1])
xrange_max <- max(xrange_hard[2], xrange_moderate[2], xrange_simple[2])
yrange_min <- min(yrange_hard[1], yrange_moderate[1], yrange_simple[1])
yrange_max <- max(yrange_hard[2], yrange_moderate[2], yrange_simple[2])





plot(agg_hard[seq(1,length(agg_hard[,1]), by=2),1],agg_hard[seq(1,length(agg_hard[,1]), by=2),2], xlab="CpG",
#ylab="Min Cost Flow Error",
     ylab= "",
     ylim=c(yrange_min, yrange_max), 
     xlim= c(xrange_min, xrange_max), 
     cex.lab = 1.2, 
     cex.axis = 1.0, 
     pch = 1, 
     cex = 0.5,
     lty = 1,
     type = "b",
     col= "red")


lines(agg_moderate[seq(1,length(agg_moderate[,1]), by=2),1],agg_moderate[seq(1,length(agg_moderate[,1]), by=2),2],
cex.axis = 1.0,
cex = 0.5,
pch = 5,
lty = 1,
type = "b",
col= "green")

lines(agg_simple[seq(1,length(agg_simple[,1]), by=2),1],agg_simple[seq(1,length(agg_simple[,1]), by=2),2],
cex.axis = 1.0,
cex = 0.5,
pch = 2,
lty = 1,
type = "b",
col= "blue")


#legend("topright", legend = c("CpG-loss", "region-loss"),
 #      pch = 19,
  #     col = c("blue","red"),
  #     #  cex = 0.5,
   #    lwd= 2.5)

dev.off()




