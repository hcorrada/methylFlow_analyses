#!/usr/bin/env Rscript
### run ./analysis.r

library("devtools")
path = "/Users/faezeh/Projects/methylFlow/src/methylFlow/methylFlowr"
#path = "../../src/methylFlow/methylFlowr"
load_all(path)

#install_github("rafalib","ririzarr")
library("rafalib")

datadir_region <- "/Users/faezeh/Projects/methylFlow/exps/bscapture/region_output"
pd <- data.frame(subject=rep(4:6,2), status=c("T", "T", "T","N", "N", "N"))
chr <- "10"

pd$dirname=sprintf("CAP_%s_%d", pd$status, pd$subject)


objs_region <- lapply(seq(len=nrow(pd)), function(i) {
  curdir <- file.path(datadir_region, pd$dirname[i],chr)
  print(curdir)
  read.methylflow.dir(curdir, pd$dirname[i], has.header=TRUE)
})

names(objs_region) <- c("T1", "T2", "T3", "N1", "N2", "N3")



objs_region <- lapply(objs_region, processMethylpats)


filteredObjs_region <- lapply(objs_region, mfFilterBy, minComponentCoverage=100, minComponentWidth=86)
names(filteredObjs_region) <- names(objs_region)

figdir <- "/Users/faezeh/Projects/methylFlow/exps/bscapture/figs_region/"
figdir_region <- file.path(figdir, chr)

if (!file.exists(figdir_region)) dir.create(figdir_region)


pdf(file.path(figdir_region,"fragment_length.pdf"),height=8,width=8)
mypar(1,1)

widths <- lapply(filteredObjs_region, function(obj) width(components(obj)))
ncpgs <- lapply(filteredObjs_region, function(obj) ncpgs(obj, level="pattern", summary="max"))

names(ncpgs) <- names(widths) <- gsub("CAP_", "", names(widths))
boxplot(widths, main="reconstructed fragment size")
boxplot(ncpgs, main="number of cpgs in reconstructed fragments", cex.lab =0.2)
dev.off()






####
filteredObjs_region <- lapply(filteredObjs_region, mfFilterBy, minNumberOfPatterns=1)

npats_region <- sapply(filteredObjs_region, npatterns, by.component=TRUE)
compCoverage_region <- sapply(filteredObjs_region, counts, level="component", kind="raw")


pdf(file.path(figdir_region,"patt_by_coverage.pdf"), height=4, width=6)
mypar(2,3)

for (i in seq(along=filteredObjs_region)) {
  plot(compCoverage_region[[i]], npats_region[[i]],main=names(filteredObjs_region)[i],xlab="coverage",ylab="num. patterns")
}
dev.off()



rawCpgs_region <- lapply(filteredObjs_region, makeCpgGR, kind="raw")
estimatedCpgs_region <- lapply(filteredObjs_region, makeCpgGR, kind="estimated")

betaGR_region <- lapply(seq(along=rawCpgs_region), function(i) {
  rawGR <- rawCpgs_region[[i]]
  estGR <- estimatedCpgs_region[[i]]
  #    readGR <- readCpgs[[i]]
  
  olaps <- findOverlaps(rawGR, estGR)
  gr <- rawGR[queryHits(olaps),]
  out <- gr
  mcols(out) <- NULL
  
  out$rawBeta <- gr$Beta
  out$rawCov <- gr$Cov
  
  gr <- estGR[subjectHits(olaps),]
  out$estimatedBeta <- gr$Beta
  
  #    olaps <- findOverlaps(readGR, out)
  #   out$readBeta <- readGR$Beta[queryHits(olaps)]
  out
})

names(betaGR_region) <- names(filteredObjs_region)
pdf(file.path(figdir_region,"meth_percentage.pdf"), height=7, width=9)
mypar(2,3)

for (i in seq(along=betaGR_region)) {
  gr <- betaGR_region[[i]]
  plot(gr$rawBeta,gr$estimatedBeta,
       bty='l',
       main=names(filteredObjs_region)[i],
       ylab="pattern methyl Percentage",
       xlab="region methyl Precentage",
       cex = 0.5,
       cex.lab=.9)
}
dev.off()



err_region = c()
for (i in seq(along=betaGR_region)) {
  gr <- betaGR_region[[i]]
  err_region[i] <- sqrt(sum(gr$rawBeta - gr$estimatedBeta)^2/length(betaGR_region[[i]]))
}

pd$error_region <- err_region



RRSE <- function(y_true, y_pred){
  RRSE <- sqrt(sum((y_true-y_pred)^2)/sum((y_true-mean(y_true))^2))
  return(RRSE)
}


RMSE <- function(y_true, y_pred){
  RMSE <- sqrt(mean((y_true-y_pred)^2))
  return(RMSE)
}

for (i in seq(along=betaGR_region)) {
  gr <- betaGR_region[[i]]
  
  print(i)
  #logit_raw <- pmin(19,pmax(1,qlogis(gr$rawBeta, location = 10, scale = 1, lower.tail = TRUE, log.p = FALSE)))
  
  #logit_estimated <- pmin(19,pmax(1,qlogis(gr$estimatedBeta, location = 10, scale = 1, lower.tail = TRUE, log.p = FALSE)))
  
  
  logit_raw <- pmin(10,pmax(-10,qlogis(gr$rawBeta, location = 0, scale = 1, lower.tail = TRUE, log.p = FALSE)))
  
  logit_estimated <- pmin(10,pmax(-10,qlogis(gr$estimatedBeta, location = 0, scale = 1, lower.tail = TRUE, log.p = FALSE)))
  
  RRSE(gr$rawBeta, gr$estimatedBeta)
  #print("RRSE = ")
  #print("RRSE")
  
  RRSE(logit_raw, logit_estimated)
  #print("logit RRSE = ")
  #print("RRSE")
  
  RMSE(gr$rawBeta, gr$estimatedBeta)
  #print("RMSE = ")
  #print("RMSE")
  
  corr = cor(gr$rawBeta, gr$estimatedBeta)
  print("corr = ")
  print(corr)
  
  cor_logit = cor(logit_raw, logit_estimated)
  print("logit corr = ")
  print(cor_logit)
}




allBeta <- betaGR_region[[1]]
#only take first sample
tmp <- allBeta$estimatedBeta

mcols(allBeta) <- NULL
mcols(allBeta)[[names(filteredObjs_region)[1]]] <- tmp

for (i in seq(2,length(betaGR_region))) {
  gr <- betaGR_region[[i]]
  olaps <- findOverlaps(gr, allBeta)
  
  sname <- names(filteredObjs_region)[i]
  mcols(allBeta)[[sname]] <- NA
  mcols(allBeta)[[sname]][subjectHits(olaps)] <- gr$estimatedBeta[queryHits(olaps)]
}


library(bumphunter)
clusters <- clusterMaker(seqnames(allBeta), start(allBeta), maxGap=1000)

#install.packages("doParallel", repos="http://R-Forge.R-project.org")
library(doParallel)
registerDoParallel(cores=2)
betamat <- as.matrix(mcols(allBeta))

bumps <- bumphunter(betamat, model.matrix(~1+status,data=pd),
                    pos=start(allBeta), cluster=clusters, B=0, smooth=FALSE, cutoff=.2)

plotRegion <- GRanges(chr, IRanges(start=bumps$tab$start[1], bumps$tab$end[1]))
#plotRegion <- plotRegion * -8

#plotRegion <- GRanges("13", IRanges(start=111127500,end=111128600))

regionBeta <- subsetByOverlaps(allBeta, plotRegion)


betamat <- as.matrix(mcols(regionBeta))
Index <- splitit(pd$status)
betamns <- sapply(Index, function(ind) rowMeans(betamat[,ind]))

fits <- lapply(1:2, function(i) loess(betamns[,i]~start(regionBeta)))

pdf(file.path(figdir_region,"countKey.pdf"), width=1,height=3)
mypar(1,1)
p <- seq(0,1,len=11)
keyCols <- colorRampPalette(brewer.pal(5,"Blues"))(11)
#keyCols <- rgb(grayrgb[1],grayrgb[2],grayrgb[3],alpha=255*p,max=255)
plot(1,1,type="n",xlim=c(0,1),ylim=c(0,length(keyCols)),bty="n",
     xaxt="n",yaxt="n",xlab="",ylab="")
rect(xleft=0,xright=1,
     ybottom=seq(along=keyCols)-1,
     ytop=seq(along=keyCols),col=keyCols)
axis(side=2,at=seq(along=keyCols)-.5,p,las=2)
dev.off()


load_all(path)
pdf(file.path(figdir_region, "region_sum_v4.pdf"),width=8,height=6)

layoutmat <- matrix(c(10,rep(1,6),11,12,13,rep(2,5),rep(3,5),rep(4,5),rep(5,5),rep(6,5),rep(7,5), rep(8,3),rep(9,5)),nc=1)
#layoutmat <- matrix(c(10,8,rep(1,4),9,rep(2,5),rep(3,3),rep(4,5),rep(5,3),rep(6,5),rep(7,3)),nc=2)

layout(layoutmat)
par(mar=c(0, 3, 0, 1.1))
par(oma=c(0,0,0,0))

xlim <- c(start(plotRegion), end(plotRegion)+70)

plot(0, type="n", xlim=xlim, ylim=c(0,1), ylab="",xaxt="n",yaxt="n")
lat <- pretty(xlim,n=5)
axis(side=1, at=lat, cex.axis=.8, mgp=c(1,.3,0))
axis(side=2, at=c(0, 0.5, 1),las=1,cex.axis=.9, mgp=c(1,.4,0))
title(xlab="genomic position (chr13)", ylab="DNA methylation", cex.lab=0.7)

for (i in seq(along=filteredObjs_region)) {
  points(start(regionBeta), mcols(regionBeta)[[names(filteredObjs_region[i])]], col=ifelse(pd$status[i]=="N",1,2), cex=.3, pch=19)
}

for (i in 1:2) {
  curve(predict(fits[[i]], x), from=min(start(regionBeta)), to=max(start(regionBeta)), col=i,add=TRUE,lwd=1.3)
}
legend("bottomright", c("Normal","Tumor"), pch=22,pt.bg=1:2)

indexes <- c(1,4,2,5,3,6)
#indexes <- c(3)
for (i in indexes) {
  # plotPatterns(filteredObjs_region[[i]], plotRegion,
  #              xaxt="n", yaxt="n", ylab="", xlab="",
  #             bty="n")
  #    mtext(side=3,paste(gsub("CAP_","",names(filteredObjs)[i])), cex=.7)
  keep <- components(filteredObjs_region[[i]]) %over% plotRegion
  
  npats <- npatterns(filteredObjs_region[[i]], by.component=TRUE)[keep]
  cids <- components(filteredObjs_region[[i]])$cid[keep]
  
  npatterns <- max(npats)
  plot(0, type="n", xlim=c(start(plotRegion) , end(plotRegion)+70), ylim=c(0,npatterns+1),  xaxt="n", yaxt="n", ylab="", xlab="",
       bty="n")
  orig_npatterns <- npatterns
  
  for (j in seq(along=cids)) {
    cid <- cids[j]
    patterns <- patterns(filteredObjs_region[[i]])[patterns(filteredObjs_region[[i]])$cid == cid,]
    maxabd = max(patterns$abundance)
    patterns <- patterns[patterns$abundance > 0.1*maxabd]
    o <- order(patterns$abundance, decreasing = FALSE)
    patterns <- patterns[o,]
    npatterns <- min(length(patterns), 5)
    m <- length(patterns)
    n <- length(patterns) - npatterns + 1
    patterns_plot <- patterns[n:m,]
    #patterns <- patterns[n-npatterns+1:n,]
    scale <- orig_npatterns/npatterns
    #npatterns <- length(patterns)
    #npatterns <- 4
  #  if(length(grep("N",names(filteredObjs_region[i])))> 0 )
    scale2 <- npatterns / 5
    shift <- (5 - npatterns)/2 - 0.2/scale2
    segments(min(start(patterns))-1, (shift+seq(along=patterns))*scale*scale2,
             max(end(patterns))+2, (shift+seq(along=patterns))*scale*scale2)
    
    #        text(end(patterns), seq(along=patterns),
    #            labels=sprintf("%.3f", patterns$abundance), pos=4,cex=.8)
    
   # scaledAbundance <- patterns$abundance / max(patterns$abundance)
    scaledAbundance <- patterns_plot$abundance / sum(patterns$abundance)
    quantizedAbundance <- cut(scaledAbundance, breaks=seq(0,1,len=11))
    
    #cols <- colorRampPalette(brewer.pal(5,"Blues"))(10)
    cols <- colorRampPalette(c("#EFF3FF", "#2163A7", "#08519C"))(10)
    #cols <- colorRampPalette(c("green","red"))(10)
    col <- cols[quantizedAbundance]
    
    
    rect(max(end(patterns))+10,(shift+(seq(along=patterns)-0.5))*scale*scale2,max(end(patterns)) + 40,(shift+(seq(along=patterns)+0.5))*scale*scale2,
         col=col)
   # text(max(end(patterns))+45, (npatterns), label="abundance of ", cex=.7,srt=90)
    text(max(end(patterns))+45, (shift+3*(npatterns)/4)*scale*scale2, label= names(filteredObjs_region[i]) ,cex=.7,srt=90)
    #        points(end(patterns)+10, seq(along=patterns), pch=22, bg=col)
   #text(max(end(patterns))+60, seq(along=patterns), label=as.character(quantizedAbundance))
    locs <- rep(start(patterns), patterns$ncpgs) +
      unlist(patterns$locs)
    ylocs <- rep(seq(along=patterns), patterns$ncpgs)
    methyl <- unlist(patterns_plot$meth)
    col  <- ifelse(methyl=="M", "black", "white")
    points(locs,(shift+ylocs)*scale*scale2,pch=21,bg=col)
  
  }
  
 
}

p <- seq(0,1,len=11)
#keyCols <- colorRampPalette(brewer.pal(5,"Blues"))(10)
keyCols <- colorRampPalette(c("#EFF3FF", "#2163A7", "#08519C"))(10)
#keyCols <- colorRampPalette(c("green","red"))(10)
#keyCols <- rgb(grayrgb[1],grayrgb[2],grayrgb[3],alpha=255*p,max=255)
plot(1,1,type="n",xlim=c(-0.5,length(keyCols)),ylim=c(0,1),bty="n",
     xaxt="n",yaxt="n",xlab="",ylab="")
#rect(xleft=0,xright=1,
#    ybottom=seq(along=keyCols)-1,
#   ytop=seq(along=keyCols),col=keyCols)
rect(xleft=seq(along=keyCols)/2+1,xright=seq(along=keyCols)/2+1.5,
     ybottom=0.1,
     ytop=0.6,col=keyCols)


axis(side=1,at=p*5+1.5,p,las=2, cex=0.8)

dev.off()

