#!/usr/bin/env Rscript

### run ./ReadLengthPlot.r


library("RColorBrewer")
#library(ggplot2)

start= 0.4
step= 0.2
end = 0.8

count = (end - start) / step + 1



#dir <- "/cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/"
wdir <- "/Users/faezeh/Projects/methylFlow/exps/compelet/"
#data[1]= 1
#wdir <- getwd();
print(wdir)
##### reading files ##################

  dir_hard <- paste(wdir,"hard-Auto/",sep="");

  
  avg_hard <- read.table(paste(dir_hard,"evalAvg.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
  
  mcf_hard <- read.table(paste(dir_hard,"mcf.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)




  dir_moderate <- paste(wdir,"moderate-Auto/",sep="");
  
  avg_moderate <- read.table(paste(dir_moderate,"evalAvg.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
  
  mcf_moderate <- read.table(paste(dir_moderate,"mcf.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
  

  

  dir_simple <- paste(wdir,"simple-Auto/",sep="");

  avg_simple <- read.table(paste(dir_simple,"evalAvg.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
  
  mcf_simple <- read.table(paste(dir_simple,"mcf.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)


#dir <- paste(wdir, "/allInOneFig/", sep="")
dir <- file.path(wdir, "allInOne")

####### different plots for readLength ##############################################


xrange_hard <- range(avg_hard$methylCallError)
xrange_moderate <- range(avg_moderate$methylCallError)
xrange_simple <- range(avg_simple$methylCallError)

yrange_hard <- range(avg_hard$abdncError)
yrange_moderate <- range(avg_moderate$abdncError)
yrange_simple <- range(avg_simple$abdncError)


#### plot the abundance Error vs methylcall for simple

print("plot abundance Error vs methylCall error vs different thresholds for simple run of readLength")
pdf(paste(dir,"abdVmethylVthrVreadLenght_Simple.pdf",sep=""), width =3, height =2.5, pointsize=8)
par(mar= c(5,5,2,2))

# get the range for the x and y axis


xrange_min <- min(xrange_hard[1], xrange_moderate[1], xrange_simple[1])
xrange_max <- max(xrange_hard[2], xrange_moderate[2], xrange_simple[2])
yrange_min <- min(yrange_hard[1], yrange_moderate[1], yrange_simple[1])
yrange_max <- max(yrange_hard[2], yrange_moderate[2], yrange_simple[2])

#ntrees <- length(unique(readLengthAvg$threshold))
ntrees <- count

# set up the plot
plot(0, 0,
     pch = "",
#yaxt='n',
ylim = c(yrange_min - 0.1, yrange_max + 0.1),
xlim = c(xrange_min - 0.05, xrange_max + 0.1),
xlab="Average MethylCall Error",
ylab="Average Abundance Error",
# ylab="Abundance Error",
     cex.lab= 1.2,
     #     cex.axis = 0.5
)


# add lines

#loess_fit <- loess(avg_simple$abdncError[sel] ~ avg_simple$methylCallError[sel], avg_simple)
#lines(avg_simple$methylCallError[sel], predict(loess_fit), col =  colors[i])

  points(avg_simple$methylCallError,
         avg_simple$abdncError,
         col = "dark green",
         #  pch = pchs[i],
         cex = 0.5,
         type = "b",
         lty = 3,
         lwd = 1.0
  )
  
  lines(avg_moderate$methylCallError,
  avg_moderate$abdncError,
  col = "blue",
  # pch = pchs[i],
  cex = 0.5,
  type = "b",
  lty = 3,
  lwd = 1.0
  )


lines(avg_hard$methylCallError,
avg_hard$abdncError,
col = "red",
#pch = pchs[i],
cex = 0.5,
type = "b",
lty = 3,
lwd = 1.0
)

legend("topright", legend = c("simple", "moderate", "hard"),
#pch = pchs,
col = c("green", "blue", "red"),
#  cex = 0.5,
cex = 0.5,
lty = 3,
pch = 1,
lwd = 1.0)


  dev.off()

