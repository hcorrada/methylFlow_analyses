
#!/usr/bin/env Rscript

### run ./cpgPlot.r par1 par2

### par1 = 0 > Auto-lambda
### par1 = 1 > Non-Auto - lambda is hard coded

### par2 = 0 > simple
### par2 = 1 > moderate
### par2 = 2 > Hard


start= 0.1
step= 0.05
end = 0.4

count = (end - start) / step + 1


data <- commandArgs(T)
print(data)
#dir <- "/cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/"
dir <- "/Users/faezeh/Desktop/project/methylFlow_analyses/cpg/"
data[1]= 0
data[2] = 0
##### reading files ##################
if (data[1] == "0"){
  if ( data[2] == "2"){
    print("Hard Setting Plot")
    CpGAvg <- read.table(paste(dir,"hard-Auto/evalAvg.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
    
    mcfCpG <- read.table(paste(dir,"hard-Auto/mcf.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
    
    #dir <- "/cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/hard-Auto/"
    dir <- "/Users/faezeh/Desktop/project/methylFlow_analyses/cpg/hard-Auto/"
    
  }
  
  if ( data[2] == "1"){
    print("Moderate Setting Plot")
    CpGAvg <- read.table(paste(dir,"moderate-Auto/evalAvg.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
    
    mcfCpG <- read.table(paste(dir,"moderate-Auto/mcf.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
    
    #dir <- "/cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/moderate-Auto/"
    dir <- "/Users/faezeh/Desktop/project/methylFlow_analyses/cpg/moderate-Auto/"
    
    
  }
  if ( data[2] == "0"){
    print("Simple Setting Plot")
    CpGAvg <- read.table(paste(dir,"simple-Auto/evalAvg.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
    
    mcfCpG <- read.table(paste(dir,"simple-Auto/mcf.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
    
    #dir <- "/cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/simple-Auto/"
    dir <- "/Users/faezeh/Desktop/project/methylFlow_analyses/cpg/simple-Auto/"
    
    
  }
}

if (data[1] == "1"){
  if ( data[2] == "2"){
    print("Hard Setting Plot")
    CpGAvg <- read.table(paste(dir,"hard/evalAvg.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
    
    mcfCpG <- read.table(paste(dir,"hard/mcf.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
    
    #dir <- "/cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/hard/"
    dir <- "/Users/faezeh/Desktop/project/methylFlow_analyses/cpg/hard/"
  }
  
  if ( data[2] == "1"){
    print("Moderate Setting Plot")
    CpGAvg <- read.table(paste(dir,"moderate/evalAvg.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
    
    mcfCpG <- read.table(paste(dir,"moderate/mcf.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
    
    #dir <- "/cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/moderate/"
    dir <- "/Users/faezeh/Desktop/project/methylFlow_analyses/cpg/moderate/"
    
  }
  if ( data[2] == "0"){
    print("Simple Setting Plot")
    CpGAvg <- read.table(paste(dir,"simple/evalAvg.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
    
    mcfCpG <- read.table(paste(dir,"simple/mcf.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
    
    #dir <- "/cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/simple/"
    dir <- "/Users/faezeh/Desktop/project/methylFlow_analyses/cpg/simple/"
    
  }
}


############## different plots for differnet CpG rate #################################################################


#### plot the abundance Error for different CpG rate

print("plot abundance Error vs CpG rate")
pdf(paste(dir,"abdVCpG.pdf",sep=""))
par(mar= c(5,5,2,2))

# get the range for the x and y axis
xrange <- range(CpGAvg$var)
yrange <- range(CpGAvg$abdncError)
ntrees <- length(unique(CpGAvg$threshold))
# set up the plot
plot(0, 0,
     pch = "",
     ylim = c(yrange[1], yrange[2] + .1),
     xlim = xrange,
     xlab="CpG",
     ylab="Abundance Error",
     cex.lab= 2,
     cex.axis = 1.5)

ltys = seq(1:ntrees)
#colors <- rainbow(ntrees)
mypalette<-brewer.pal(ntrees,"Dark2")
pchs = seq(0,6)
# add lines
for (i in 1:ntrees) {
  j = step * (i-1) + start
  print(j)
  
  sel <- which(CpGAvg$threshold == j)
  points(CpGAvg$var[sel],
         CpGAvg$abdncError[sel],
         col = mypalette[i],
         pch = pchs[i],
         cex = 0.5,
         type = "b",
         lwd = 2.5)
}

legend("topright", legend = seq(start,end, by = step),
       pch = pchs,
       col = mypalette,
       #  cex = 0.5,
       lwd= 2.5)

dev.off()
# cex scale the size
#pch = 16 is circle
#lx <- seq(xrange[1] + 30 ,0.4*(xrange[2]-xrange[1])+ xrange[1], length.out=ntrees) + 5
#ly <- rep(yrange[2] + 0.1, ntrees)
#points(lx, ly,
#col = colors[1:28],
#pch = 16, cex=1,
#lwd=2)

#txt <- unique(CpGAvg$threshold)
#sel1 <- c(1, (1:count))

#text(xrange[1], yrange[2]+.1, "Legend", cex=1.3, pos=4,
#lwd=2)
#text(lx[sel1], ly[sel1], pos = 1, offset = 2, txt[sel1], cex=1.3, srt=90,
#lwd=2)




#### plot the Methyl Call  Error for different CpG rate

print("plot methyl call Error vs CpG rate ")
pdf(paste(dir, "methylVCpG.pdf",sep=""))

par(mar= c(5,5,2,2))

# get the range for the x and y axis
xrange <- range(CpGAvg$var)
yrange <- range(CpGAvg$methylCallError)
ntrees <- length(unique(CpGAvg$threshold))
# set up the plot
plot(0, 0,
     pch = "",
     ylim = c(yrange[1], yrange[2] + .1),
     xlim = xrange,
     xlab="CpG",
     ylab="Methyl Call Error" ,
     cex.lab= 2,
     cex.axis = 1.5)


# add lines
ltys = seq(1:ntrees)
#colors <- rainbow(ntrees)
mypalette<-brewer.pal(ntrees,"Dark2")
pchs = seq(0,6)
for (i in 1:ntrees) {
  j = step * (i-1) + start
  print(j)
  
  sel <- which(CpGAvg$threshold == j)
  points(CpGAvg$var[sel],
         CpGAvg$methylCallError[sel],
         col = mypalette[i],
         pch = pchs[i],
         cex = 0.5,
         type = "b",
         lwd = 2.5)
}

legend("topright", legend = seq(start,end, by = step),
       pch = pchs,
       col = mypalette,
       #  cex = 0.5,
       lwd= 2.5)

dev.off()

# cex scale the size
#pch = 16 is circle
#lx <- seq(xrange[1] + 30 ,0.4*(xrange[2]-xrange[1])+ xrange[1], length.out=ntrees) + 5
#ly <- rep(yrange[2] + .1, ntrees)
#points(lx, ly,
#col = colors[1:28],
#pch = 16, cex=1,
#lwd=2)

#txt <- unique(CpGAvg$threshold)
#sel1 <- c(1, (1:count))

#text(xrange[1], yrange[2]+.1, "Legend", cex=1.3, pos=4,
#lwd=2)
#text(lx[sel1], ly[sel1] ,pos = 1, offset = 2, txt[sel1], cex=1.3, srt=90,
#lwd=2)


#### plot #TP for different CpG rate

print("plot TP vs CpG rate ")
pdf(paste(dir, "TPVCpG.pdf",sep=""))


# get the range for the x and y axis
xrange <- range(CpGAvg$var)
yrange <- range(CpGAvg$TP)
ntrees <- length(unique(CpGAvg$threshold))
# set up the plot
plot(0, 0,
     pch = "",
     ylim = c(yrange[1], yrange[2] + .1),
     xlim = xrange,
     xlab="CpG",
     ylab="TP",
     cex.lab= 2,
     cex.axis = 1.5)

ltys = seq(1:ntrees)
#colors <- rainbow(ntrees)
mypalette<-brewer.pal(ntrees,"Dark2")
pchs = seq(0,6)

# add lines
for (i in 1:ntrees) {
  j = step * (i-1) + start
  print(j)
  
  sel <- which(CpGAvg$threshold == j)
  lines(CpGAvg$var[sel],
        CpGAvg$TP[sel],
        col = mypalette[i],
        pch = pchs[i],
        cex = 0.5,
        type = "b",
        lwd = 2.5)
}

legend("bottomright", legend = seq(start,end, by = step),
       pch = pchs,
       col = mypalette,
       #  cex = 0.5,
       lwd= 2.5)

dev.off()
# cex scale the size
#pch = 16 is circle
#lx <- seq(xrange[1] + 25 ,0.3*(xrange[2]-xrange[1])+ xrange[1], length.out=ntrees) + 5
#ly <- rep(yrange[2] + .1, ntrees)
#points(lx, ly,
#col = colors[1:28],
#pch = 16, cex=1,
#lwd=2)

#txt <- unique(CpGAvg$threshold)
#sel1 <- c(1, (1:count))

#text(xrange[1], yrange[2]+.1, "Legend", cex=1, pos=4,
#lwd=2)
#text(lx[sel1], ly[sel1] - 0.02, txt[sel1], cex=1, srt=90,
#lwd=2)





######## plot sensitivity  for different CpG  ###############
#### TP / (TP + FN)


print("plot sensitivity vs CpG")
pdf(paste(dir,"sensitivityVCpG.pdf",sep=""))

sensitivity = CpGAvg$TP/(CpGAvg$TP + CpGAvg$FN)

# get the range for the x and y axis
xrange <- range(CpGAvg$var)
yrange <- range(sensitivity)

yy  <- ifelse(is.na(yrange[2]) , 0, yrange[2])
yrange[2] <- yy

yy  <- ifelse(is.na(yrange[1]) , 0, yrange[1])
yrange[1] <- yy

ntrees <- length(unique(CpGAvg$threshold))
# set up the plot
plot(0, 0,
     pch = "",
     ylim = c(yrange[1], yrange[2] + .1),
     xlim = xrange,
     xlab="CpG",
     ylab="sensitivity" ,
     cex.lab= 2,
     cex.axis = 1.5)


ltys = seq(1:ntrees)
#colors <- rainbow(ntrees)
mypalette<-brewer.pal(ntrees,"Dark2")
pchs = seq(0,6)

# add lines
for (i in 1:ntrees) {
  j = step * (i-1) + start
  print(j)
  
  sel <- which(CpGAvg$threshold == j)
  lines(CpGAvg$var[sel],
        sensitivity[sel],
        col = mypalette[i],
        pch = pchs[i],
        cex = 0.5,
        type = "b",
        lwd = 2.5)
}

legend("bottomright", legend = seq(start,end, by = step),
       pch = pchs,
       col = mypalette,
       #  cex = 0.5,
       lwd= 2.5)

dev.off()

# cex scale the size
#pch = 16 is circle
#lx <- seq(xrange[1] + 25 ,0.3*(xrange[2]-xrange[1])+ xrange[1], length.out=ntrees) + 5
#ly <- rep(yrange[2] + .1, ntrees)
#points(lx, ly,
#col = colors[1:28],
#pch = 16, cex=1,
#lwd=2)

#txt <- unique(CpGAvg$threshold)
#sel1 <- c(1, (1:count))

#text(xrange[1], yrange[2]+.1, "Legend", cex=1, pos=4 ,
#lwd=2)
#text(lx[sel1], ly[sel1] - 0.02, txt[sel1], cex=1, srt=90,
#lwd=2)




######## plot Precision(Positive Predictive Rate)  for different CpG ###############
##### TP / (TP + FP )


print("plot precision vs CpG")
pdf(paste(dir,"precisionVCpG.pdf",sep=""))

precision = rep(0, length(CpGAvg$var));

sel <- which(CpGAvg$TP + CpGAvg$FP != 0)
precision[sel] = CpGAvg$TP[sel]/(CpGAvg$TP[sel] + CpGAvg$FP[sel])


# get the range for the x and y axis
xrange <- range(CpGAvg$var)
yrange <- range(precision)
ntrees <- length(unique(CpGAvg$threshold))
# set up the plot
plot(0, 0,
     pch = "",
     ylim = c(yrange[1], yrange[2] + .1),
     xlim = xrange,
     xlab="CpG",
     ylab="precision",
     cex.lab= 2,
     cex.axis = 1.5 )

ltys = seq(1:ntrees)
#colors <- rainbow(ntrees)
mypalette<-brewer.pal(ntrees,"Dark2")
pchs = seq(0,6)

# add lines
for (i in 1:ntrees) {
  j = step * (i-1) + start
  print(j)
  
  sel <- which(CpGAvg$threshold == j)
  lines(CpGAvg$var[sel],
        precision[sel],
        col = mypalette[i],
        pch = pchs[i],
        cex = 0.5,
        type = "b",
        lwd = 2.5)
}

legend("topright", legend = seq(start,end, by = step),
       pch = pchs,
       col = mypalette,
       #  cex = 0.5,
       lwd= 2.5)

dev.off()

# cex scale the size
#pch = 16 is circle
#lx <- seq(xrange[1] + 25 ,0.3*(xrange[2]-xrange[1])+ xrange[1], length.out=ntrees) + 5
#ly <- rep(yrange[2] + .1, ntrees)
#points(lx, ly,
#col = colors[1:28],
#pch = 16, cex=1,
#lwd=2)

#txt <- unique(CpGAvg$threshold)
#sel1 <- c(1, (1:count))

#text(xrange[1], yrange[2]+.1, "Legend", cex=1, pos=4,
#lwd=2)
#text(lx[sel1], ly[sel1] - 0.02, txt[sel1], cex=1, srt=90,
#lwd=2)


#dev.off()


######## plot FDR false discovery rate  for CpG ###############


print("plot false discovery rate vs CpG ")
pdf(paste(dir,"FDRVCpG.pdf",sep=""))

FDR = CpGAvg$FP/(CpGAvg$TP + CpGAvg$FP)

# get the range for the x and y axis
xrange <- range(CpGAvg$var)
yrange <- range(FDR)
ntrees <- length(unique(CpGAvg$threshold))
# set up the plot
plot(0, 0,
     pch = "",
     ylim = c(yrange[1], yrange[2] + .1),
     xlim = xrange,
     xlab="CpG",
     ylab="FDR" ,
     cex.lab= 2,
     cex.axis = 1.5)

ltys = seq(1:ntrees)
#colors <- rainbow(ntrees)
mypalette<-brewer.pal(ntrees,"Dark2")
pchs = seq(0,6)

# add lines
for (i in 1:ntrees) {
  j = step * (i-1) + start
  print(j)
  
  sel <- which(CpGAvg$threshold == j)
  lines(CpGAvg$var[sel],
        FDR[sel],
        col = mypalette[i],
        pch = pchs[i],
        cex = 0.5,
        type = "b",
        lwd = 2.5)
}

legend("topright", legend = seq(start,end, by = step),
       pch = pchs,
       col = mypalette,
       #  cex = 0.5,
       lwd= 2.5)

dev.off()
# cex scale the size
#pch = 16 is circle
#lx <- seq(xrange[1] + 25 ,0.3*(xrange[2]-xrange[1])+ xrange[1], length.out=ntrees) + 5
#ly <- rep(yrange[2] + .1, ntrees)
#points(lx, ly,
#col = colors[1:28],
#pch = 16, cex=1,
#lwd=2)

#txt <- unique(CpGAvg$threshold)
#sel1 <- c(1, (1:count))

#text(xrange[1], yrange[2]+.1, "Legend", cex=1, pos=4,
#lwd=2)
#text(lx[sel1], ly[sel1] - 0.02, txt[sel1], cex=1, srt=90,
#lwd=2)


#dev.off()


####### plot min cost flow error for differnet CpG ########


print("plot min cost flow error vs CpG")
pdf(paste(dir,"mcfCpG.pdf",sep=""))
par(mar= c(5,5,2,2))

agg = aggregate(mcfCpG$minCostFlow, list(CpG = mcfCpG$var), FUN =  mean)

plot(agg, ylab="Min Cost Flow Error", cex.lab = 2, cex.axis = 1.5, pch = 19, col= "blue")

dev.off()



###### plot FP  for CpG ###########



print("plot FP vs CpG ")
pdf(paste(dir,"FPVCpG.pdf",sep=""))

FP = (CpGAvg$FP)

# get the range for the x and y axis
xrange <- range(CpGAvg$var)
yrange <- range(FP)
ntrees <- length(unique(CpGAvg$threshold))
# set up the plot
plot(0, 0,
     pch = "",
     ylim = c(yrange[1], yrange[2] + 10),
     xlim = xrange,
     xlab="CpG",
     ylab="FP" ,
     cex.lab= 1.5,
     cex.axis = 1.5)

ltys = seq(1:ntrees)
#colors <- rainbow(ntrees)
mypalette<-brewer.pal(ntrees,"Dark2")
pchs = seq(0,6)

# add lines
for (i in 1:ntrees) {
  j = step * (i-1) + start
  print(j)
  
  sel <- which(CpGAvg$threshold == j)
  lines(CpGAvg$var[sel],
        FP[sel],
        col = mypalette[i],
        pch = pchs[i],
        cex = 0.5,
        type = "b",
        lwd = 2.5)
}

legend("topright", legend = seq(start,end, by = step),
       pch = pchs,
       col = mypalette,
       #  cex = 0.5,
       lwd= 2.5)

dev.off()
# cex scale the size
#pch = 16 is circle
#lx <- seq(xrange[1] + 25 ,0.3*(xrange[2]-xrange[1])+ xrange[1], length.out=ntrees) + 5
#ly <- rep(yrange[2] + 0.15*(yrange[2]-yrange[1]), ntrees)
#points(lx, ly,
#col = colors[1:28],
#pch = 16, cex=1,
#lwd=2)

#txt <- unique(CpGAvg$threshold)
#sel1 <- c(1, (1:count))

#text(xrange[1], yrange[2]+.1, "Legend", cex=1, pos=4,
#lwd=2)
#text(lx[sel1], ly[sel1] - 0.01*(yrange[2]-yrange[1]), txt[sel1], cex=1, srt=90,
#lwd=2)


#dev.off()



