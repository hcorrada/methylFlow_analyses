#!/usr/bin/env Rscript
### run ./sam.r par1

### par1 = input directory: readLength/coverage/cpg
### par2 = type: simple/moderate/hard


library("devtools")
path = "~/Projects/methylFlow/src/methylFlow/methylFlowr"
load_all(path)

#data <- commandArgs(T)

dir <- "/Users/faezeh/Projects/methylFLow/exps/compare_old_new_methylFlow_simulation/"
#install_github("rafalib","ririzarr")
library("rafalib")

resultdir <- file.path(paste(dir, "comp_meth_pct",sep=""))
print(resultdir)

pd <- data.frame(type=c("cpg_loss","cpg_loss","cpg_loss", "region_loss", "region_loss", "region_loss"), status=c("Simple", "Moderate", "Hard", "Simple", "Moderate", "Hard"))
pd$dirname=paste(pd$status,"/", pd$type, sep="")


objs <- lapply(seq(len=nrow(pd)), function(i) {
  curdir <- file.path(resultdir, pd$dirname[i])
  read.methylflow.dir(curdir, pd$dirname[i],has.header=TRUE)
})

names(objs) <- pd$dirname
objs <- lapply(objs, processMethylpats)

filteredObjs <- lapply(objs, mfFilterBy, minComponentCoverage=1, minComponentWidth=10)
names(filteredObjs) <- names(objs)


figdir <- file.path(resultdir, "figs")
if (!file.exists(figdir)) dir.create(figdir)

filteredObjs <- lapply(filteredObjs, mfFilterBy, minNumberOfPatterns=1)

rawCpgs <- lapply(filteredObjs, makeCpgGR, kind="raw")
estimatedCpgs <- lapply(filteredObjs, makeCpgGR, kind="estimated")


betaGR <- lapply(seq(along=rawCpgs), function(i) {
  rawGR <- rawCpgs[[i]]
  estGR <- estimatedCpgs[[i]]
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

pdf(file.path(figdir,"meth_pct_1.pdf"), height=7, width=9)
mypar(2,3)

for (i in seq(along=betaGR)) {
  gr <- betaGR[[i]]
  plot(gr$rawBeta,gr$estimatedBeta,
       bty='l',
       main=names(filteredObjs)[i],
       ylab="pattern methyl Percentage",
       xlab="region methyl Precentage",
       cex = 0.5,
       cex.lab=2)
}


dev.off()

pdf(file.path(figdir,"meth_pct_2.pdf"), height=3, width=9)
mypar(1,3)

for (i in c(1,2,3)) {
  gr_cpg <- betaGR[[i]]
  plot(gr_cpg$rawBeta,gr_cpg$estimatedBeta,
       bty='l',
       main=pd$status[i],
       ylab="pattern methyl Percentage",
       xlab="region methyl Precentage",
       col="blue",
       cex = 0.5,
       cex.lab=1.5)
  
  gr_region <- betaGR[[i+3]]
  points(gr_region$rawBeta,gr_region$estimatedBeta,
       bty='l',
       col="red")
}

legend("bottomright", c("cpg_loss","region_loss"), col=c("blue", "red"), pch=1)


dev.off()


pdf(file.path(figdir,"meth_pct_3.pdf"), height=3, width=9)
mypar(1,3)

for (i in c(4,5,6)) {
    gr <- betaGR[[i]]
    plot(gr$rawBeta,gr$estimatedBeta,
    bty='l',
    main=names(filteredObjs)[i],
    ylab="pattern methyl Percentage",
    xlab="region methyl Precentage",
    cex = 0.5,
    cex.lab=2)
}


dev.off()



err <- c()
for (i in seq(along=betaGR)) {
  gr <- betaGR[[i]]
  err[i] <- sqrt(sum(gr$rawBeta - gr$estimatedBeta)^2/length(betaGR[[i]]))
}

pd$error <- err
error <- matrix(err,nrow=3)
colnames(error) <- c("cpg_loss", "region_loss")
rownames(error) <-  c("Simple","Moderate","Hard")
error <- as.table(error)
error

pdf(file.path(figdir,"meth_sqr_err_loss.pdf"))
barplot(error,legend=T, beside=T)
dev.off()

#################################

objsNew <- lapply(seq(len=1), function(i) {
  curDirNew <- file.path(resultdir, newDir)
  read.methylflow.dir(curDirNew, has.header=TRUE)
})

objsNew <- lapply(objsNew, processMethylpats)
filteredObjsNew <- lapply(objsNew, mfFilterBy, minComponentCoverage=1, minComponentWidth=10)
names(filteredObjsNew) <- names(objsNew)

###############
objsOld <- lapply(seq(len=1), function(i) {
  curDirOld <- file.path(resultdir, oldDir)
  print(curDirOld)
  read.methylflow.dir(curDirOld, has.header=TRUE)
})

#names(objsOld) <- file.path(dir,datadir)

objsOld <- lapply(objsOld, processMethylpats)
filteredObjsOld <- lapply(objsOld, mfFilterBy, minComponentCoverage=1, minComponentWidth=10)
names(filteredObjsOld) <- names(objsOld)

#################################

figdir <- file.path(resultdir,"figs2")

if (!file.exists(figdir)) dir.create(figdir)

pdf(file.path(figdir,"fragment_length_new.pdf"),height=9,width=9)
#mypar(1,1)

widthsNew <- lapply(filteredObjsNew, function(obj) width(components(obj)))
ncpgsNew <- lapply(filteredObjsNew, function(obj) ncpgs(obj, level="pattern", summary="max"))

#names(ncpgs) <- names(widths) <- gsub("CAP_", "", names(widths))
boxplot(widthsNew, main="reconstructed fragment size")
boxplot(ncpgsNew, main="number of cpgs in reconstructed fragments")
dev.off()


figdir <- file.path(resultdir,"figs2")

if (!file.exists(figdir)) dir.create(figdir)

pdf(file.path(figdir,"fragment_length_old.pdf"),height=9,width=9)
#mypar(1,1)

widthsOld <- lapply(filteredObjsOld, function(obj) width(components(obj)))
ncpgsOld <- lapply(filteredObjsOld, function(obj) ncpgs(obj, level="pattern", summary="max"))

#names(ncpgs) <- names(widths) <- gsub("CAP_", "", names(widths))
boxplot(widthsOld, main="reconstructed fragment size")
boxplot(ncpgsOld, main="number of cpgs in reconstructed fragments")
dev.off()


  
  #for (i in seq(along=betaGR)) {
  #    gr <- betaGR[[i]]
  #    plot(gr$rawBeta,gr$readBeta,
  #         bty='l',
  #         main=names(filteredObjs)[i],
  #         ylab="read methyl Percentage",
  #         xlab="region methyl Precentage",
  #         cex = 0.2,
  #         cex.lab=.9)
  #}
#



