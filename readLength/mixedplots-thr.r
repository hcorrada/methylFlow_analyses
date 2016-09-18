#!/usr/bin/env Rscript

### run ./ReadLengthPlot.r


library("RColorBrewer")
#library(ggplot2)

start= 0.4
step= 0.2
end = 0.8

count = (end - start) / step + 1



#dir <- "/cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/"
wdir <- "/Users/faezeh/Projects/methylFlow/exps/readLength/"
#data[1]= 1
#wdir <- getwd();
print(wdir)
##### reading files ##################

  dir_hard <- paste(wdir,"hard-Auto/",sep="");

  
  readLengthAvg_hard <- read.table(paste(dir_hard,"evalAvg.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
  
  mcfReadLength_hard <- read.table(paste(dir_hard,"mcf.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)




  dir_moderate <- paste(wdir,"moderate-Auto/",sep="");
  
  readLengthAvg_moderate <- read.table(paste(dir_moderate,"evalAvg.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
  
  mcfReadLength_moderate <- read.table(paste(dir_moderate,"mcf.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
  

  

  dir_simple <- paste(wdir,"simple-Auto/",sep="");

  readLengthAvg_simple <- read.table(paste(dir_simple,"evalAvg.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
  
  mcfReadLength_simple <- read.table(paste(dir_simple,"mcf.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)


#dir <- paste(wdir, "/allInOneFig/", sep="")
dir <- file.path(wdir, "mixed_thr_")

####### different plots for readLength ##############################################
sel_s  <- with(readLengthAvg_simple,which(threshold %in% c("0.4","0.6","0.8")))
sel_m  <- with(readLengthAvg_moderate,which(threshold %in% c("0.4","0.6","0.8")))
sel_h  <- with(readLengthAvg_hard,which(threshold %in% c("0.4","0.6","0.8")))


xrange_hard <- range(readLengthAvg_hard$methylCallError[sel_h])
xrange_moderate <- range(readLengthAvg_moderate$methylCallError[sel_m])
xrange_simple <- range(readLengthAvg_simple$methylCallError[sel_s])

yrange_hard <- range(readLengthAvg_hard$abdncError[sel_h])
yrange_moderate <- range(readLengthAvg_moderate$abdncError[sel_m])
yrange_simple <- range(readLengthAvg_simple$abdncError[sel_s])


#### plot the abundance Error vs methylcall for simple

print("plot abundance Error vs methylCall error vs different thresholds for simple run of readLength")
pdf(paste(dir,"abdVmethylVthrVreadLenght_Simple.pdf",sep=""), width =3, height =2.5, pointsize=8)
par(mar= c(5,5,2,2))

# get the range for the x and y axis


#xrange_min <- min(xrange_hard[1], xrange_moderate[1], xrange_simple[1])
#xrange_max <- max(xrange_hard[2], xrange_moderate[2], xrange_simple[2])
#yrange_min <- min(yrange_hard[1], yrange_moderate[1], yrange_simple[1])
#yrange_max <- max(yrange_hard[2], yrange_moderate[2], yrange_simple[2])

#ntrees <- length(unique(readLengthAvg$threshold))
ntrees <- count

# set up the plot
plot(0, 0,
     pch = "",
#yaxt='n',
     ylim = yrange_simple,
     xlim = xrange_simple,
     xlab="MethylCall Error",
ylab="Abundance Error",
# ylab="Abundance Error",
     cex.lab= 1.2,
     #     cex.axis = 0.5
)

ltys = seq(1:ntrees)
colors <-brewer.pal(2*ntrees,"Dark2")
pchs = seq(0,ntrees-1)
# add lines
for (i in 1:ntrees) {
  j = step * (i-1) + start
  print(j)
  
  sel <- which(abs(readLengthAvg_simple$threshold - j) < 0.0001)
  

#loess_fit <- loess(readLengthAvg_simple$abdncError[sel] ~ readLengthAvg_simple$methylCallError[sel], readLengthAvg_simple)
#lines(readLengthAvg_simple$methylCallError[sel], predict(loess_fit), col =  colors[i])

  points(readLengthAvg_simple$methylCallError[sel[seq(1,length(sel),by=2)]],
         readLengthAvg_simple$abdncError[sel[seq(1,length(sel),by=2)]],
         col = colors[i],
         pch = pchs[i],
         cex = 0.5
         #     type = "b",
         #    lty = 3,
         #lwd = 1.0
  )
}
  dev.off()
  
  
  print("plot abundance Error vs methylCall error vs different thresholds for moderate run of readLength")
  pdf(paste(dir,"abdVmethylVthrVreadLenght_Moderate.pdf",sep=""), width =3, height =2.5, pointsize=8)
  par(mar= c(5,5,2,2))
  
  # get the range for the x and y axis
  
  
  #xrange_min <- min(xrange_hard[1], xrange_moderate[1], xrange_simple[1])
  #xrange_max <- max(xrange_hard[2], xrange_moderate[2], xrange_simple[2])
  #yrange_min <- min(yrange_hard[1], yrange_moderate[1], yrange_simple[1])
  #yrange_max <- max(yrange_hard[2], yrange_moderate[2], yrange_simple[2])
  
  #ntrees <- length(unique(readLengthAvg$threshold))
  ntrees <- count
  
  # set up the plot
  plot(0, 0,
  pch = "",
  #yaxt='n',
  ylim = yrange_moderate,
  xlim = xrange_moderate,
  xlab="MethylCall Error",
  ylab="Abundance Error",
  # ylab="Abundance Error",
  cex.lab= 1.2,
  #     cex.axis = 0.5
  )
  
  ltys = seq(1:ntrees)
  colors <-brewer.pal(2*ntrees,"Dark2")
  pchs = seq(0,ntrees-1)
  # add lines
  for (i in 1:ntrees) {
      j = step * (i-1) + start
      print(j)
      
  
  sel <- which(abs(readLengthAvg_moderate$threshold - j) < 0.0001)
  
  #loess_fit <- loess(readLengthAvg_moderate$abdncError[sel] ~ readLengthAvg_moderate$methylCallError[sel], readLengthAvg_moderate)
  #lines(readLengthAvg_moderate$methylCallError[sel], predict(loess_fit), col =  colors[i])

  points(readLengthAvg_moderate$methylCallError[sel[seq(1,length(sel),by=2)]],
         readLengthAvg_moderate$abdncError[sel[seq(1,length(sel),by=2)]],
         col = colors[i],
         pch = pchs[i],
         cex = 0.5
         # type = "b",
         #   lty = 3,
         #lwd = 1.0
  )
  }
  
  dev.off()
  
  
  print("plot abundance Error vs methylCall error vs different thresholds for hard run of readLength")
  pdf(paste(dir,"abdVmethylVthrVreadLenght_Hard.pdf",sep=""), width =3, height =2.5, pointsize=8)
  par(mar= c(5,5,2,2))
  
  # get the range for the x and y axis
  
  
  #xrange_min <- min(xrange_hard[1], xrange_moderate[1], xrange_simple[1])
  #xrange_max <- max(xrange_hard[2], xrange_moderate[2], xrange_simple[2])
  #yrange_min <- min(yrange_hard[1], yrange_moderate[1], yrange_simple[1])
  #yrange_max <- max(yrange_hard[2], yrange_moderate[2], yrange_simple[2])
  
  #ntrees <- length(unique(readLengthAvg$threshold))
  ntrees <- count
  
  # set up the plot
  plot(0, 0,
  pch = "",
  #yaxt='n',
  ylim = yrange_hard,
  xlim = xrange_hard,
  xlab="MethylCall Error",
  ylab="Abundance Error",
  # ylab="Abundance Error",
  cex.lab= 1.2,
  #     cex.axis = 0.5
  )
  
  ltys = seq(1:ntrees)
#colors <- rainbow(ntrees)
  colors <-brewer.pal(2*ntrees,"Dark2")
  pchs = seq(0,ntrees-1)
  # add lines
  for (i in 1:ntrees) {
      j = step * (i-1) + start
      print(j)
      
  
  sel <- which(abs(readLengthAvg_hard$threshold - j) < 0.0001)
  
  #loess_fit <- loess(readLengthAvg_hard$abdncError[sel] ~ readLengthAvg_hard$methylCallError[sel], readLengthAvg_hard)
  #lines(readLengthAvg_hard$methylCallError[sel], predict(loess_fit), col =  colors[i])
 
  points(readLengthAvg_hard$methylCallError[sel[seq(1,length(sel),by=2)]],
            readLengthAvg_hard$abdncError[sel[seq(1,length(sel),by=2)]],
            col = colors[i],
            pch = pchs[i],
            cex = 0.5
            #     type = "b",
            #        lty = 3,
            #    lwd = 1.0
  )
  }
legend("topright", legend= c("thr = 0.4","thr = 0.6","thr = 0.8"),
#, "Moderate, thr = 0.2","Moderate, thr = 0.4","Moderate, thr = 0.6", "Moderate, thr = 0.8","Hard, thr = 0.2","Hard, thr = 0.4","Hard, thr = 0.6","Hard, thr = 0.8"
#title = "thresholds",
pch = rep(pchs,3),
col = colors[1:ntrees],
cex = 0.7,
pt.cex = 0.7,
lty = 3,
lwd= 1.0)

  
  dev.off()


