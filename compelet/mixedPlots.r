#!/usr/bin/env Rscript

### run ./ReadLengthPlot.r


library("RColorBrewer")
#library(ggplot2)

data <- commandArgs(T)
#data[1]== 0 >>>> for different threshold
#data[1]== 1 >>>> for different lambda
#data[1]== 0 >>>> for different noise

#data[1]=1

#dir <- "/cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/"
wdir <- "/Users/faezeh/Projects/methylFlow/exps/compelet/"
#data[1]= 1
#wdir <- getwd();
print(wdir)
##### reading files ##################
if (data[1] == "0"){
    
    start= 0.4
    step= 0.2
    end = 0.8
    
    count = (end - start) / step + 1


  dir_hard <- paste(wdir,"hard-Auto/",sep="");

  
  avg_hard <- read.table(paste(dir_hard,"evalAvg.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
  
  mcf_hard <- read.table(paste(dir_hard,"mcf.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)




  dir_moderate <- paste(wdir,"moderate-Auto/",sep="");
  
  avg_moderate <- read.table(paste(dir_moderate,"evalAvg.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
  
  mcf_moderate <- read.table(paste(dir_moderate,"mcf.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
  

  

  dir_simple <- paste(wdir,"simple-Auto/",sep="");

  avg_simple <- read.table(paste(dir_simple,"evalAvg.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
  
  mcf_simple <- read.table(paste(dir_simple,"mcf.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
  
  
  dir <- file.path(wdir, "allInOne")


xrange_hard <- range(avg_hard$methylCallError)
xrange_moderate <- range(avg_moderate$methylCallError)
xrange_simple <- range(avg_simple$methylCallError)

yrange_hard <- range(avg_hard$abdncError)
yrange_moderate <- range(avg_moderate$abdncError)
yrange_simple <- range(avg_simple$abdncError)


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
col = c(" dark green", "blue", "red"),
#  cex = 0.5,
cex = 0.5,
lty = 3,
pch = 1,
lwd = 1.0)


dev.off()



}

####################### lambda ##################

if (data[1] == "1"){
    start= 0.4
    step= 0.2
    end = 0.8
    

    count = (end - start) / step + 1
    


dir_hard <- paste(wdir,"hard-lambda/",sep="");


avg_hard <- read.table(paste(dir_hard,"evalAvg.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)

mcf_hard <- read.table(paste(dir_hard,"mcf.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)




dir_moderate <- paste(wdir,"moderate-lambda/",sep="");

avg_moderate <- read.table(paste(dir_moderate,"evalAvg.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)

mcf_moderate <- read.table(paste(dir_moderate,"mcf.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)




dir_simple <- paste(wdir,"simple-lambda/",sep="");

avg_simple <- read.table(paste(dir_simple,"evalAvg.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)

mcf_simple <- read.table(paste(dir_simple,"mcf.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)



dir <- file.path(wdir, "lambda_")

print("plot min cost flow error vs Coverage simple")
pdf(paste(dir,"mcf_simple.pdf",sep=""))
par(mar= c(5,5,2,2))

agg = aggregate(mcf_simple$minCostFlow, list(lambda = mcf_simple$var), FUN =  mean)

plot(agg, log="x", ylab="Min Cost Flow Error", cex.lab = 2, cex.axis = 1.5, pch = 19, col= "blue")

dev.off()



print("plot min cost flow error vs Coverage moderare")
pdf(paste(dir,"mcf_moderate.pdf",sep=""))
par(mar= c(5,5,2,2))

agg = aggregate(mcf_moderate$minCostFlow, list(lambda = mcf_moderate$var), FUN =  mean)

plot(agg, log="x", ylab="Min Cost Flow Error", cex.lab = 2, cex.axis = 1.5, pch = 19, col= "blue")

dev.off()


print("plot min cost flow error vs Coverage hard")
pdf(paste(dir,"mcf_hard.pdf",sep=""))
par(mar= c(5,5,2,2))

agg = aggregate(mcf_hard$minCostFlow, list(lambda = mcf_hard$var), FUN =  mean)

plot(agg, log="x", ylab="Min Cost Flow Error", cex.lab = 2, cex.axis = 1.5, pch = 19, col= "blue")

dev.off()




print("plot min cost flow error vs Coverage")
pdf(paste(dir,"mcf_all.pdf",sep=""), width =3, height =2.5, pointsize=8)
par(mar= c(5,5,2,2))


agg_hard = aggregate(mcf_hard$minCostFlow, list(lambda = mcf_hard$var), FUN =  mean)
agg_moderate = aggregate(mcf_moderate$minCostFlow, list(lambda = mcf_moderate$var), FUN =  mean)
agg_simple = aggregate(mcf_simple$minCostFlow, list(lambda = mcf_simple$var), FUN =  mean)


xrange_hard <- range(agg_hard$lambda)
xrange_moderate <- range(agg_moderate$lambda)
xrange_simple <- range(agg_simple$lambda)

yrange_hard <- range(agg_hard$x)
yrange_moderate <- range(agg_moderate$x)
yrange_simple <- range(agg_simple$x)

xrange_min <- min(xrange_hard[1], xrange_moderate[1], xrange_simple[1])
xrange_max <- max(xrange_hard[2], xrange_moderate[2], xrange_simple[2])
yrange_min <- min(yrange_hard[1], yrange_moderate[1], yrange_simple[1])
yrange_max <- max(yrange_hard[2], yrange_moderate[2], yrange_simple[2])

yrange_min <- 0;
yrange_max <- 1;




plot(agg_hard[seq(1,length(agg_hard[,1]), by=2),1],agg_hard[seq(1,length(agg_hard[,1]), by=2),2], log="x", xlab="Lambda",
ylab="Min Cost Flow Error",
#ylim=c(yrange_min, yrange_max),
ylim=c(0,0.8),
xlim= c(xrange_min, xrange_max),
xaxt="n",
cex.lab = 1.2,
cex.axis = 1.0,
pch = 1,
cex = 0.5,
lty = 1,
type = "b",
col= "red")

ticks <- seq(-2, 2, by=1)
labels <- sapply(ticks, function(i) as.expression(bquote(10^ .(i))))
axis(1, at=c(0.01, 0.1, 1, 10, 100), labels=labels)

lines(agg_moderate[seq(1,length(agg_moderate[,1]), by=2),1],agg_moderate[seq(1,length(agg_moderate[,1]), by=2),2],
cex.axis = 1.0,
cex = 0.5,
pch = 5,
lty = 1,
type = "b",
col= "black")

lines(agg_simple[seq(1,length(agg_simple[,1]), by=2),1],agg_simple[seq(1,length(agg_simple[,1]), by=2),2],
cex.axis = 1.0,
cex = 0.5,
pch = 2,
lty = 1,
type = "b",
col= "blue")


legend("topright", legend = c("Hard", "Moderate", "Simple"),
pch = c(1,5,2),
col = c("red","black", "blue"),
cex = 0.85,
pt.cex = 0.5,
lty = 1,
lwd= 1.0)


dev.off()






}


if (data[1] == "2"){
    start= 0.4
    step= 0.2
    end = 0.8
    

    count = (end - start) / step + 1
    
    
    
    dir_hard <- paste(wdir,"hard-noise/",sep="");
    
    
    avg_hard <- read.table(paste(dir_hard,"evalAvg.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
    
    mcf_hard <- read.table(paste(dir_hard,"mcf.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
    
    
    
    
    dir_moderate <- paste(wdir,"moderate-noise/",sep="");
    
    avg_moderate <- read.table(paste(dir_moderate,"evalAvg.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
    
    mcf_moderate <- read.table(paste(dir_moderate,"mcf.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
    
    
    
    
    dir_simple <- paste(wdir,"simple-noise/",sep="");
    
    avg_simple <- read.table(paste(dir_simple,"evalAvg.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
    
    mcf_simple <- read.table(paste(dir_simple,"mcf.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
    
    
    
    dir <- file.path(wdir, "noise_")
    
    
    
    print("plot min cost flow error vs Coverage simple")
    pdf(paste(dir,"mcf_simple.pdf",sep=""))
    par(mar= c(5,5,2,2))
    
    agg = aggregate(mcf_simple$minCostFlow, list(Noise = mcf_simple$var), FUN =  mean)
    
    plot(agg, ylab="Min Cost Flow Error", cex.lab = 2, cex.axis = 1.5, pch = 19, col= "blue")
    
    dev.off()
    
    
    
    print("plot min cost flow error vs Coverage moderare")
    pdf(paste(dir,"mcf_moderate.pdf",sep=""))
    par(mar= c(5,5,2,2))
    
    agg = aggregate(mcf_moderate$minCostFlow, list(Noise = mcf_moderate$var), FUN =  mean)
    
    plot(agg, ylab="Min Cost Flow Error", cex.lab = 2, cex.axis = 1.5, pch = 19, col= "blue")
    
    dev.off()
    
    
    print("plot min cost flow error vs Coverage hard")
    pdf(paste(dir,"mcf_hard.pdf",sep=""))
    par(mar= c(5,5,2,2))
    
    agg = aggregate(mcf_hard$minCostFlow, list(Noise = mcf_hard$var), FUN =  mean)
    
    plot(agg, ylab="Min Cost Flow Error", cex.lab = 2, cex.axis = 1.5, pch = 19, col= "blue")
    
    dev.off()
    
    
    
    
    print("plot min cost flow error vs Coverage")
    pdf(paste(dir,"mcf_all.pdf",sep=""), width =3, height =2.5, pointsize=8)
    par(mar= c(5,5,2,2))
    
    
    agg_hard = aggregate(mcf_hard$minCostFlow, list(Noise = mcf_hard$var), FUN =  mean)
    agg_moderate = aggregate(mcf_moderate$minCostFlow, list(Noise = mcf_moderate$var), FUN =  mean)
    agg_simple = aggregate(mcf_simple$minCostFlow, list(Noise = mcf_simple$var), FUN =  mean)
    
    
    xrange_hard <- range(agg_hard$Noise)
    xrange_moderate <- range(agg_moderate$Noise)
    xrange_simple <- range(agg_simple$Noise)
    
    yrange_hard <- range(agg_hard$x)
    yrange_moderate <- range(agg_moderate$x)
    yrange_simple <- range(agg_simple$x)
    
    xrange_min <- min(xrange_hard[1], xrange_moderate[1], xrange_simple[1])
    xrange_max <- max(xrange_hard[2], xrange_moderate[2], xrange_simple[2])
    yrange_min <- min(yrange_hard[1], yrange_moderate[1], yrange_simple[1])
    yrange_max <- max(yrange_hard[2], yrange_moderate[2], yrange_simple[2])
    
    #yrange_min <- 0;
    #yrange_max <- 1;
    
    
    
    
    plot(agg_hard[seq(1,length(agg_hard[,1]), by=2),1],agg_hard[seq(1,length(agg_hard[,1]), by=2),2], xlab="Noise",
    ylab="Min Cost Flow Error",
    #ylim=c(yrange_min, yrange_max),
    ylim=c(0,0.8),
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
    col= "black")
    
    lines(agg_simple[seq(1,length(agg_simple[,1]), by=2),1],agg_simple[seq(1,length(agg_simple[,1]), by=2),2],
    cex.axis = 1.0,
    cex = 0.5,
    pch = 2,
    lty = 1,
    type = "b",
    col= "blue")
    
    
    legend("topright", legend = c("Hard", "Moderate", "Simple"),
    pch = c(1,5,2),
    col = c("red","black", "blue"),
    cex = 0.85,
    pt.cex = 0.5,
    lty = 1,
    lwd= 1.0)
    
    
    dev.off()
    
    

    
}

