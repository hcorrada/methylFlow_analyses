#!/usr/bin/env Rscript

### run ./coveragePlot.r par1 par2

### par1 = 0 > Auto-lambda
### par1 = 1 > Non-Auto - lambda is hard coded

### par2 = 0 > simple
### par2 = 1 > moderate
### par2 = 2 > Hard

library("RColorBrewer")

start= 0.2
step= 0.2
end = 0.8

count = (end - start) / step + 1


data <- commandArgs(T)
print(data)
#dir <- "/Users/faezeh/Desktop/project/methylFlow_analyses/cpg/"

#dir <- "/cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/"
#dir <- "/Users/faezeh/Desktop/project/methylFlow_analyses/coverage/"
wdir <- getwd();

#data[1]= 0
#data[2] = 2
##### reading files ##################
if (data[1] == "0"){
    if ( data[2] == "2"){
        dir <- paste(wdir,"/hard-Auto/",sep="");
        print(dir);
        print("Hard Setting Plot")
        coverageAvg <- read.table(paste(dir,"evalAvg.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
        
        mcfCoverage <- read.table(paste(dir,"mcf.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
        
    
    }
    
    if ( data[2] == "1"){
        
        print("Moderate Setting Plot")
        dir <- paste(wdir,"/moderate-Auto/",sep="");
        print(dir);

        coverageAvg <- read.table(paste(dir,"evalAvg.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
        
        mcfCoverage <- read.table(paste(dir,"mcf.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
        
    }
    if ( data[2] == "0"){
        
        print("Simple Setting Plot")
        dir <- paste(wdir,"/simple-Auto/",sep="");
        print(dir);
        
        coverageAvg <- read.table(paste(dir,"evalAvg.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
        
        mcfCoverage <- read.table(paste(dir,"mcf.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
        
        
        
    }
}

if (data[1] == "1"){
if ( data[2] == "2"){
    
    print("Hard Setting Plot")
    dir <- paste(wdir,"/hard/",sep="");
    print(dir);
    coverageAvg <- read.table(paste(dir,"evalAvg.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
    
    mcfCoverage <- read.table(paste(dir,"mcf.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
    

}

if ( data[2] == "1"){
    
    print("Moderate Setting Plot")
    dir <- paste(wdir,"/moderate/",sep="");
    print(dir);
    coverageAvg <- read.table(paste(dir,"evalAvg.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
    
    mcfCoverage <- read.table(paste(dir,"mcf.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
    
}
if ( data[2] == "0"){
    
    print("Simple Setting Plot")
    dir <- paste(wdir,"/simple/",sep="");
    print(dir);
    coverageAvg <- read.table(paste(dir,"evalAvg.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
    
    mcfCoverage <- read.table(paste(dir,"mcf.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
    

}
}


############## different plots for differnet coverage rate #################################################################


#### plot the abundance Error for different coverage rate

print("plot abundance Error vs coverage rate")
pdf(paste(dir,"abdVcoverage.pdf",sep=""))
par(mar= c(5,5,2,2))

# get the range for the x and y axis
xrange <- range(coverageAvg$var)
yrange <- range(coverageAvg$abdncError)
ntrees <- length(unique(coverageAvg$threshold))
# set up the plot
plot(0, 0,
pch = "",
ylim = c(yrange[1], yrange[2] + .1),
xlim = xrange,
xlab="coverage",
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
    
    sel <- which(abs(coverageAvg$threshold - j) < 0.0001)
    points(coverageAvg$var[sel],
    coverageAvg$abdncError[sel],
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

#txt <- unique(coverageAvg$threshold)
#sel1 <- c(1, (1:count))

#text(xrange[1], yrange[2]+.1, "Legend", cex=1.3, pos=4,
#lwd=2)
#text(lx[sel1], ly[sel1], pos = 1, offset = 2, txt[sel1], cex=1.3, srt=90,
#lwd=2)




#### plot the Methyl Call  Error for different coverage rate

print("plot methyl call Error vs coverage rate ")
pdf(paste(dir, "methylVcoverage.pdf",sep=""))

par(mar= c(5,5,2,2))

# get the range for the x and y axis
xrange <- range(coverageAvg$var)
yrange <- range(coverageAvg$methylCallError)
ntrees <- length(unique(coverageAvg$threshold))
# set up the plot
plot(0, 0,
pch = "",
ylim = c(yrange[1], yrange[2] + .1),
xlim = xrange,
xlab="coverage",
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
    
    sel <- which(abs(coverageAvg$threshold - j) < 0.0001)
    points(coverageAvg$var[sel],
    coverageAvg$methylCallError[sel],
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

#txt <- unique(coverageAvg$threshold)
#sel1 <- c(1, (1:count))

#text(xrange[1], yrange[2]+.1, "Legend", cex=1.3, pos=4,
#lwd=2)
#text(lx[sel1], ly[sel1] ,pos = 1, offset = 2, txt[sel1], cex=1.3, srt=90,
#lwd=2)


#### plot #TP for different coverage rate

print("plot TP vs coverage rate ")
pdf(paste(dir, "TPVcoverage.pdf",sep=""))


# get the range for the x and y axis
xrange <- range(coverageAvg$var)
yrange <- range(coverageAvg$TP)
ntrees <- length(unique(coverageAvg$threshold))
# set up the plot
plot(0, 0,
pch = "",
ylim = c(yrange[1], yrange[2] + .1),
xlim = xrange,
xlab="coverage",
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
    
    sel <- which(abs(coverageAvg$threshold - j) < 0.0001)
    lines(coverageAvg$var[sel],
    coverageAvg$TP[sel],
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

#txt <- unique(coverageAvg$threshold)
#sel1 <- c(1, (1:count))

#text(xrange[1], yrange[2]+.1, "Legend", cex=1, pos=4,
#lwd=2)
#text(lx[sel1], ly[sel1] - 0.02, txt[sel1], cex=1, srt=90,
#lwd=2)





######## plot sensitivity  for different coverage  ###############
#### TP / (TP + FN)


print("plot sensitivity vs coverage")
pdf(paste(dir,"sensitivityVcoverage.pdf",sep=""))

sensitivity = coverageAvg$TP/(coverageAvg$TP + coverageAvg$FN)

# get the range for the x and y axis
xrange <- range(coverageAvg$var)
yrange <- range(sensitivity)

yy  <- ifelse(is.na(yrange[2]) , 0, yrange[2])
yrange[2] <- yy

yy  <- ifelse(is.na(yrange[1]) , 0, yrange[1])
yrange[1] <- yy

ntrees <- length(unique(coverageAvg$threshold))
# set up the plot
plot(0, 0,
pch = "",
ylim = c(yrange[1], yrange[2] + .1),
xlim = xrange,
xlab="coverage",
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
    
    sel <- which(abs(coverageAvg$threshold - j) < 0.0001)
    lines(coverageAvg$var[sel],
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

#txt <- unique(coverageAvg$threshold)
#sel1 <- c(1, (1:count))

#text(xrange[1], yrange[2]+.1, "Legend", cex=1, pos=4 ,
#lwd=2)
#text(lx[sel1], ly[sel1] - 0.02, txt[sel1], cex=1, srt=90,
#lwd=2)




######## plot Precision(Positive Predictive Rate)  for different Coverage ###############
##### TP / (TP + FP )


print("plot precision vs coverage")
pdf(paste(dir,"precisionVcoverage.pdf",sep=""))

precision = rep(0, length(coverageAvg$var));

sel <- which(coverageAvg$TP + coverageAvg$FP != 0)
precision[sel] = coverageAvg$TP[sel]/(coverageAvg$TP[sel] + coverageAvg$FP[sel])


# get the range for the x and y axis
xrange <- range(coverageAvg$var)
yrange <- range(precision)
ntrees <- length(unique(coverageAvg$threshold))
# set up the plot
plot(0, 0,
pch = "",
ylim = c(yrange[1], yrange[2] + .1),
xlim = xrange,
xlab="coverage",
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
    
    sel <- which(abs(coverageAvg$threshold - j) < 0.0001)
    lines(coverageAvg$var[sel],
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

#txt <- unique(coverageAvg$threshold)
#sel1 <- c(1, (1:count))

#text(xrange[1], yrange[2]+.1, "Legend", cex=1, pos=4,
#lwd=2)
#text(lx[sel1], ly[sel1] - 0.02, txt[sel1], cex=1, srt=90,
#lwd=2)


#dev.off()


######## plot FDR false discovery rate  for Coverage ###############


print("plot false discovery rate vs coverage ")
pdf(paste(dir,"FDRVCoverage.pdf",sep=""))

FDR = coverageAvg$FP/(coverageAvg$TP + coverageAvg$FP)

# get the range for the x and y axis
xrange <- range(coverageAvg$var)
yrange <- range(FDR)
ntrees <- length(unique(coverageAvg$threshold))
# set up the plot
plot(0, 0,
pch = "",
ylim = c(yrange[1], yrange[2] + .1),
xlim = xrange,
xlab="coverage",
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
    
    sel <- which(abs(coverageAvg$threshold - j) < 0.0001)
    lines(coverageAvg$var[sel],
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

#txt <- unique(coverageAvg$threshold)
#sel1 <- c(1, (1:count))

#text(xrange[1], yrange[2]+.1, "Legend", cex=1, pos=4,
#lwd=2)
#text(lx[sel1], ly[sel1] - 0.02, txt[sel1], cex=1, srt=90,
#lwd=2)


#dev.off()


####### plot min cost flow error for differnet coverage ########


print("plot min cost flow error vs Coverage")
pdf(paste(dir,"mcfCoverage.pdf",sep=""))
par(mar= c(5,5,2,2))

agg = aggregate(mcfCoverage$minCostFlow, list(coverage = mcfCoverage$var), FUN =  mean)

plot(agg, ylab="Min Cost Flow Error", ylim=c(0,1), cex.lab = 2, cex.axis = 1.5, pch = 19, col= "blue")

dev.off()



###### plot FP  for coverage ###########



print("plot FP vs coverage ")
pdf(paste(dir,"FPVCovergae.pdf",sep=""))

FP = (coverageAvg$FP)

# get the range for the x and y axis
xrange <- range(coverageAvg$var)
yrange <- range(FP)
ntrees <- length(unique(coverageAvg$threshold))
# set up the plot
plot(0, 0,
pch = "",
ylim = c(yrange[1], yrange[2] + 10),
xlim = xrange,
xlab="Coverage",
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
    
    sel <- which(abs(coverageAvg$threshold - j) < 0.0001)
    lines(coverageAvg$var[sel],
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

#txt <- unique(coverageAvg$threshold)
#sel1 <- c(1, (1:count))

#text(xrange[1], yrange[2]+.1, "Legend", cex=1, pos=4,
#lwd=2)
#text(lx[sel1], ly[sel1] - 0.01*(yrange[2]-yrange[1]), txt[sel1], cex=1, srt=90,
#lwd=2)


#dev.off()

