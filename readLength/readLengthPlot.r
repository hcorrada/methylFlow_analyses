#!/usr/bin/env Rscript

### ./readLengthPlot.r par1 par2

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
dir <- "/cbcb/project-scratch/fdorri/Code/methylFlow/testing/readLength/"




##### reading files ##################
if (data[1] == "0"){
if ( data[2] == "2"){
    print("Hard Setting Plot")
    readLengthAvg <- read.table(paste(dir,"hard-Auto/evalAvg.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
    
    mcfReadLength <- read.table(paste(dir,"hard-Auto/mcf.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
    
    dir <- "/cbcb/project-scratch/fdorri/Code/methylFlow/testing/readLength/hard-Auto/"
    
}

if ( data[2] == "1"){
    print("Moderate Setting Plot")
    readLengthAvg <- read.table(paste(dir,"moderate-Auto/evalAvg.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
    
    mcfReadLength <- read.table(paste(dir,"moderate-Auto/mcf.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
    
    dir <- "/cbcb/project-scratch/fdorri/Code/methylFlow/testing/readLength/moderate-Auto/"
    
    
}
if ( data[2] == "0"){
    print("Simple Setting Plot")
    readLengthAvg <- read.table(paste(dir,"simple-Auto/evalAvg.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
    
    mcfReadLength <- read.table(paste(dir,"simple-Auto/mcf.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
    
    dir <- "/cbcb/project-scratch/fdorri/Code/methylFlow/testing/readLength/simple-Auto/"
    
    
}
}

if (data[1] == "1"){
    if ( data[2] == "2"){
        print("Hard Setting Plot")
        readLengthAvg <- read.table(paste(dir,"hard/evalAvg.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
        
        mcfReadLength <- read.table(paste(dir,"hard/mcf.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
        
        dir <- "/cbcb/project-scratch/fdorri/Code/methylFlow/testing/readLength/hard/"
        
    }
    
    if ( data[2] == "1"){
        print("Moderate Setting Plot")
        readLengthAvg <- read.table(paste(dir,"moderate/evalAvg.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
        
        mcfReadLength <- read.table(paste(dir,"moderate/mcf.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
        
        dir <- "/cbcb/project-scratch/fdorri/Code/methylFlow/testing/readLength/moderate/"
        
        
    }
    if ( data[2] == "0"){
        print("Simple Setting Plot")
        readLengthAvg <- read.table(paste(dir,"simple/evalAvg.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
        
        mcfReadLength <- read.table(paste(dir,"simple/mcf.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
        
        dir <- "/cbcb/project-scratch/fdorri/Code/methylFlow/testing/readLength/simple/"
        
        
    }
}






############## different plots for differnet read Length #################################################################


#### plot the abundance Error for different read Length

print("plot abundance Error vs read Length")
pdf(paste(dir,"abdVreadLength.pdf",sep=""))
par(mar= c(5,5,2,2))


# get the range for the x and y axis
xrange <- range(readLengthAvg$var)
yrange <- range(readLengthAvg$abdncError)
ntrees <- length(unique(readLengthAvg$threshold))
# set up the plot
plot(0, 0,
pch = "",
ylim = c(yrange[1], yrange[2] + .1),
xlim = xrange,
xlab="readLength",
ylab="Abundance Error",
cex.lab= 1.5,
cex.axis = 1.5)

colors <- rainbow(ntrees)
# add lines
for (i in 1:ntrees) {
    j = step * (i-1) + start
    print(j)
    
    sel <- which(readLengthAvg$threshold == j)
    lines(readLengthAvg$var[sel],
    readLengthAvg$abdncError[sel],
    col = colors[i],
    lwd=2)
}
# cex scale the size
#pch = 16 is circle
lx <- seq(xrange[1] + 0.15*(xrange[2]-xrange[1]) ,0.4*(xrange[2]-xrange[1])+ xrange[1], length.out=ntrees) + 5
ly <- rep(yrange[2] + .1, ntrees)
points(lx, ly,
col = colors[1:28],
pch = 16, cex=1,
lwd=2)

txt <- unique(readLengthAvg$threshold)
sel1 <- c(1, (1:count))

text(xrange[1], yrange[2]+.1, "Legend", cex=1.3, pos=4,
lwd=2)
text(lx[sel1], ly[sel1],pos = 1, offset = 2, txt[sel1], cex=1.3, srt=90,
lwd=2)


dev.off()


#### plot the Methyl Call  Error for different read Length

print("plot methyl call Error vs readLength ")
pdf(paste(dir,"methylVreadLength.pdf",sep=""))

par(mar= c(5,5,2,2))

# get the range for the x and y axis
xrange <- range(readLengthAvg$var)
yrange <- range(readLengthAvg$methylCallError)
ntrees <- length(unique(readLengthAvg$threshold))
# set up the plot
par(mar=c(5, 5, 1, 0.5))

plot(0, 0,
pch = "",
ylim = c(yrange[1], yrange[2] + .1),
xlim = xrange,
xlab="readLength",
ylab="Methyl Call Error" ,
cex.lab= 2,
cex.axis = 1.5)


colors <- rainbow(ntrees)
# add lines
for (i in 1:ntrees) {
    j = step * (i-1) + start
    print(j)
    
    sel <- which(readLengthAvg$threshold == j)
    lines(readLengthAvg$var[sel],
    readLengthAvg$methylCallError[sel],
    col = colors[i],
    lwd=2)
}
# cex scale the size
#pch = 16 is circle
lx <- seq(xrange[1] + 0.15*(xrange[2]-xrange[1]) ,0.4*(xrange[2]-xrange[1])+ xrange[1], length.out=ntrees) + 5
ly <- rep(yrange[2] + .1, ntrees)
points(lx, ly,
col = colors[1:28],
pch = 16, cex=1,
lwd=2)

txt <- unique(readLengthAvg$threshold)
sel1 <- c(1, (1:count))

text(xrange[1], yrange[2]+.1, "Legend", cex=1.3, pos=4,
lwd=2)
text(lx[sel1], ly[sel1],pos = 1, offset = 2, txt[sel1], cex=1.3, srt=90,
lwd=2)


dev.off()

#### plot #TP for different read Length

print("plot TP vs readLength ")
pdf(paste(dir,"TPVreadLength.pdf",sep=""))


# get the range for the x and y axis
xrange <- range(readLengthAvg$var)
yrange <- range(readLengthAvg$TP)
ntrees <- length(unique(readLengthAvg$threshold))
# set up the plot
plot(0, 0,
pch = "",
ylim = c(yrange[1], yrange[2] + .1),
xlim = xrange,
xlab="readLength",
ylab="TP" ,
cex.lab= 1.5,
cex.axis = 1.5 )

colors <- rainbow(ntrees)
# add lines
for (i in 1:ntrees) {
    j = step * (i-1) + start
    print(j)
    
    sel <- which(readLengthAvg$threshold == j)
    lines(readLengthAvg$var[sel],
    readLengthAvg$TP[sel],
    col = colors[i],
    lwd=2)
}
# cex scale the size
#pch = 16 is circle
lx <- seq(xrange[1] + 0.15*(xrange[2]-xrange[1]) ,0.4*(xrange[2]-xrange[1])+ xrange[1], length.out=ntrees) + 5
ly <- rep(yrange[2] + .1, ntrees)
points(lx, ly,
col = colors[1:28],
pch = 16, cex=1,
lwd=2)

txt <- unique(readLengthAvg$threshold)
sel1 <- c(1, (1:count))

text(xrange[1], yrange[2]+.1, "Legend", cex=1, pos=4,
lwd=2)
text(lx[sel1], ly[sel1] - 0.02, txt[sel1], cex=1, srt=90,
lwd=2)


dev.off()




######## plot sensitivity  for different readLength  ###############
#### TP / (TP + FN)


print("plot sensitivity vs readLength")
pdf(paste(dir,"sensitivityVreadLength.pdf",sep=""))

sensitivity = readLengthAvg$TP/(readLengthAvg$TP + readLengthAvg$FN)

# get the range for the x and y axis
xrange <- range(readLengthAvg$var)
yrange <- range(sensitivity)

yy  <- ifelse(is.na(yrange[2]) , 0, yrange[2])
yrange[2] <- yy

yy  <- ifelse(is.na(yrange[1]) , 0, yrange[1])
yrange[1] <- yy

print(yrange[1])
print(yrange[2])

ntrees <- length(unique(readLengthAvg$threshold))
# set up the plot
plot(0, 0,
pch = "",
ylim = c(yrange[1], yrange[2] + .1),
xlim = xrange,
xlab="readLength",
ylab="sensitivity" ,
cex.lab= 1.5,
cex.axis = 1.5 )

colors <- rainbow(ntrees)
# add lines
for (i in 1:ntrees) {
    j = step * (i-1) + start
    print(j)
    
    sel <- which(readLengthAvg$threshold == j)
    lines(readLengthAvg$var[sel],
    sensitivity[sel],
    col = colors[i],
    lwd=2)
}
# cex scale the size
#pch = 16 is circle
lx <- seq(xrange[1] + 0.15*(xrange[2]-xrange[1]) ,0.4*(xrange[2]-xrange[1])+ xrange[1], length.out=ntrees) + 5
ly <- rep(yrange[2] + .1, ntrees)
points(lx, ly,
col = colors[1:28],
pch = 16, cex=1,
lwd=2)

txt <- unique(readLengthAvg$threshold)
sel1 <- c(1, (1:count))

text(xrange[1], yrange[2]+.1, "Legend", cex=1, pos=4,
lwd=2)
text(lx[sel1], ly[sel1] - 0.02, txt[sel1], cex=1, srt=90,
lwd=2)


dev.off()



######## plot Precision(Positive Predictive Rate)  for different readLength ###############
##### TP / (TP + FP )


print("plot precision vs readLength")
pdf(paste(dir,"precisionVreadLength.pdf",sep=""))

precision = rep(0, length(readLengthAvg$var));

sel <- which(readLengthAvg$TP + readLengthAvg$FP != 0)
precision[sel] = readLengthAvg$TP[sel]/(readLengthAvg$TP[sel] + readLengthAvg$FP[sel])


# get the range for the x and y axis
xrange <- range(readLengthAvg$var)
yrange <- range(precision)
ntrees <- length(unique(readLengthAvg$threshold))
# set up the plot
plot(0, 0,
pch = "",
ylim = c(yrange[1], yrange[2] + .1),
xlim = xrange,
xlab="readLength",
ylab="precision",
cex.lab= 1.5,
cex.axis = 1.5  )

colors <- rainbow(ntrees)
# add lines
for (i in 1:ntrees) {
    j = step * (i-1) + start
    print(j)
    
    sel <- which(readLengthAvg$threshold == j)
    lines(readLengthAvg$var[sel],
    precision[sel],
    col = colors[i],
    lwd=2)
}
# cex scale the size
#pch = 16 is circle
lx <- seq(xrange[1] + 0.15*(xrange[2]-xrange[1]) ,0.4*(xrange[2]-xrange[1])+ xrange[1], length.out=ntrees) + 5
ly <- rep(yrange[2] + .1, ntrees)
points(lx, ly,
col = colors[1:28],
pch = 16, cex=1,
lwd=2)

txt <- unique(readLengthAvg$threshold)
sel1 <- c(1, (1:count))

text(xrange[1], yrange[2]+.1, "Legend", cex=1, pos=4,
lwd=2)
text(lx[sel1], ly[sel1] - 0.02, txt[sel1], cex=1, srt=90,
lwd=2)


dev.off()


######## plot FDR false discovery rate  for readLength ###############


print("plot false discovery rate vs readLength ")
pdf(paste(dir,"FDRVreadLength.pdf",sep=""))
FDR = rep(0, length(readLengthAvg$var));

sel <- which((readLengthAvg$TP + readLengthAvg$FP) != 0)
FDR[sel] = readLengthAvg$FP[sel]/(readLengthAvg$TP[sel] + readLengthAvg$FP[sel])

# get the range for the x and y axis
xrange <- range(readLengthAvg$var)
yrange <- range(FDR)
ntrees <- length(unique(readLengthAvg$threshold))
# set up the plot
plot(0, 0,
pch = "",
ylim = c(yrange[1], yrange[2] + .1),
xlim = xrange,
xlab="readLength",
ylab="FDR" ,
cex.lab= 1.5,
cex.axis = 1.5 )

colors <- rainbow(ntrees)
# add lines
for (i in 1:ntrees) {
    j = step * (i-1) + start
    print(j)
    
    sel <- which(readLengthAvg$threshold == j)
    lines(readLengthAvg$var[sel],
    FDR[sel],
    col = colors[i],
    lwd=2)
}
# cex scale the size
#pch = 16 is circle
lx <- seq(xrange[1] + 0.15*(xrange[2]-xrange[1]) ,0.4*(xrange[2]-xrange[1])+ xrange[1], length.out=ntrees) + 5
ly <- rep(yrange[2] + .1, ntrees)
points(lx, ly,
col = colors[1:28],
pch = 16, cex=1,
lwd=2)

txt <- unique(readLengthAvg$threshold)
sel1 <- c(1, (1:count))

text(xrange[1], yrange[2]+.1, "Legend", cex=1, pos=4,
lwd=2)
text(lx[sel1], ly[sel1] - 0.02, txt[sel1], cex=1, srt=90,
lwd=2)


dev.off()


####### plot min cost flow error for differnet read length ########


print("plot min cost flow error vs ReadLength")
pdf(paste(dir,"mcfReadLength.pdf",sep=""))
par(mar= c(5,5,2,2))


agg = aggregate(mcfReadLength$minCostFlow, list(readLength = mcfReadLength$var), FUN =  mean)

plot(agg, ylab= " Min Cost Flow Error", cex.lab = 2, cex.axis = 1.5, col= "blue", pch = 19)

dev.off()

######## plot FP   for readLength ###############




print("plot FP vs Read Length ")
pdf(paste(dir,"FPVREadLength.pdf",sep=""))

FP = (readLengthAvg$FP)

# get the range for the x and y axis
xrange <- range(readLengthAvg$var)
yrange <- range(FP)
ntrees <- length(unique(readLengthAvg$threshold))
# set up the plot
plot(0, 0,
pch = "",
ylim = c(yrange[1], yrange[2] + 10),
xlim = xrange,
xlab="Read Length",
ylab="FP" ,
cex.lab= 1.5,
cex.axis = 1.5 )

colors <- rainbow(ntrees)
# add lines
for (i in 1:ntrees) {
    j = step * (i-1) + start
    print(j)
    
    sel <- which(readLengthAvg$threshold == j)
    lines(readLengthAvg$var[sel],
    FP[sel],
    col = colors[i],
    lwd=2)
}
# cex scale the size
#pch = 16 is circle
lx <- seq(xrange[1] + 0.15*(xrange[2]-xrange[1]) ,0.4*(xrange[2]-xrange[1])+ xrange[1], length.out=ntrees) + 5
ly <- rep(yrange[2] + .1, ntrees)
points(lx, ly,
col = colors[1:28],
pch = 16, cex=1,
lwd=2)

txt <- unique(readLengthAvg$threshold)
sel1 <- c(1, (1:count))

text(xrange[1], yrange[2]+.1, "Legend", cex=1, pos=4,
lwd=2)
text(lx[sel1], ly[sel1] - 0.02, txt[sel1], cex=1, srt=90,
lwd=2)


dev.off()




