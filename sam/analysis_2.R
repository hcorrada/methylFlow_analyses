#!/usr/bin/env Rscript
### run ./analysis.r

library("devtools")
path = "/Users/faezeh/Projects/methylFlow/src/methylFlow/methylFlowr"
#path = "../../src/methylFlow/methylFlowr"
load_all(path)

#install_github("rafalib","ririzarr")
library("rafalib")

datadir_cpg <- "/Users/faezeh/Projects/methylFlow/exps/bscapture/cpg_output"
datadir_region <- "/Users/faezeh/Projects/methylFlow/exps/bscapture/region_output"
pd <- data.frame(subject=rep(4:6,2), status=c("T", "T", "T","N", "N", "N"))
j <- "13"

pd$dirname=sprintf("CAP_%s_%d", pd$status, pd$subject)

objs_cpg <- lapply(seq(len=nrow(pd)), function(i) {
    curdir <- file.path(datadir_cpg, pd$dirname[i],j)
    print(curdir)
    read.methylflow.dir(curdir, pd$dirname[i], has.header=TRUE)
})

objs_region <- lapply(seq(len=nrow(pd)), function(i) {
    curdir <- file.path(datadir_region, pd$dirname[i],j)
    print(curdir)
    read.methylflow.dir(curdir, pd$dirname[i], has.header=TRUE)
})

names(objs_cpg) <- pd$dirname
names(objs_region) <- pd$dirname


objs_cpg <- lapply(objs_cpg, processMethylpats)
objs_region <- lapply(objs_region, processMethylpats)

filteredObjs_cpg <- lapply(objs_cpg, mfFilterBy, minComponentCoverage=100, minComponentWidth=86)
names(filteredObjs_cpg) <- names(objs_cpg)


filteredObjs_region <- lapply(objs_region, mfFilterBy, minComponentCoverage=100, minComponentWidth=86)
names(filteredObjs_region) <- names(objs_region)

figdir <- "/Users/faezeh/Projects/methylFlow/exps/bscapture/"
figdir_region <- file.path(figdir, "figs_region")
figdir_cpg <- file.path(figdir, "figs_cpg")
if (!file.exists(figdir_region)) dir.create(figdir_region)
if (!file.exists(figdir_cpg)) dir.create(figdir_cpg)


pdf(file.path(figdir_region,"fragment_length.pdf"),height=8,width=8)
mypar(1,1)

widths <- lapply(filteredObjs_region, function(obj) width(components(obj)))
ncpgs <- lapply(filteredObjs_region, function(obj) ncpgs(obj, level="pattern", summary="max"))

names(ncpgs) <- names(widths) <- gsub("CAP_", "", names(widths))
boxplot(widths, main="reconstructed fragment size")
boxplot(ncpgs, main="number of cpgs in reconstructed fragments", cex.lab =2)
dev.off()




pdf(file.path(figdir_cpg,"fragment_length.pdf"),height=8,width=8)
mypar(1,1)

widths <- lapply(filteredObjs_cpg, function(obj) width(components(obj)))
ncpgs <- lapply(filteredObjs_cpg, function(obj) ncpgs(obj, level="pattern", summary="max"))

names(ncpgs) <- names(widths) <- gsub("CAP_", "", names(widths))
boxplot(widths, main="reconstructed fragment size")
boxplot(ncpgs, main="number of cpgs in reconstructed fragments", cex.lab =2)
dev.off()


####
filteredObjs_region <- lapply(filteredObjs_region, mfFilterBy, minNumberOfPatterns=1)
filteredObjs_cpg <- lapply(filteredObjs_cpg, mfFilterBy, minNumberOfPatterns=1)


npats_cpg <- sapply(filteredObjs_cpg, npatterns, by.component=TRUE)
compCoverage_cpg <- sapply(filteredObjs_cpg, counts, level="component", kind="raw")

npats_region <- sapply(filteredObjs_region, npatterns, by.component=TRUE)
compCoverage_region <- sapply(filteredObjs_region, counts, level="component", kind="raw")


pdf(file.path(figdir_cpg,"patt_by_coverage.pdf"), height=4, width=6)
mypar(2,3)
  
for (i in seq(along=filteredObjs_cpg)) {
    plot(compCoverage_cpg[[i]], npats_cpg[[i]],main=names(filteredObjs_cpg)[i],xlab="coverage",ylab="num. patterns")
}
dev.off()


pdf(file.path(figdir_region,"patt_by_coverage.pdf"), height=4, width=6)
mypar(2,3)

for (i in seq(along=filteredObjs_region)) {
    plot(compCoverage_region[[i]], npats_region[[i]],main=names(filteredObjs_region)[i],xlab="coverage",ylab="num. patterns")
}
dev.off()



rawCpgs_cpg <- lapply(filteredObjs_cpg, makeCpgGR, kind="raw")
estimatedCpgs_cpg <- lapply(filteredObjs_cpg, makeCpgGR, kind="estimated")


rawCpgs_region <- lapply(filteredObjs_region, makeCpgGR, kind="raw")
estimatedCpgs_region <- lapply(filteredObjs_region, makeCpgGR, kind="estimated")



#readGRs <- lapply(seq(len=nrow(pd)), function(i) {
#    curdir <- file.path("reads", pd$dirname[i])
#    tab <- read.delim(file.path(curdir, "13.methylation.withsub.tsv"),
#                      header=FALSE,stringsAsFactors=FALSE)
#    names(tab) <- c("id", "start", "width", "cstrand", "methylpat", "subst")
#    tab$end <- tab$start + tab$width - 1
#    tab$chr <- "13"
#    tab$width <- NULL
#    gr <- tab2gr(tab)
#    gr
#})
#
#.parseMethylpats <- function(x) {
#    tmp <- strsplit(x, ",")
#    ncpgs <- ifelse(x=="*", 0, sapply(tmp,length))
#    tmp2 <- lapply(seq(along=x), function(i) {
#        if (ncpgs[i] == 0) return(NULL)
#        strsplit(tmp[[i]], ":")
#    })
#    
#    locs <- lapply(tmp2, function(y) {
#      if (is.null(y)) return(0)
#      as.integer(sapply(y,"[",1))
#    })
#
#    meth <- lapply(tmp2, function(y) {
#        if (is.null(y)) return("")
#        sapply(y,"[",2)
#    })
#    list(ncpgs=ncpgs,locs=locs,meth=meth)
#}
#
#readGRs <- lapply(readGRs, function(x) {
#    tmp <- .parseMethylpats(x$methylpat)
#    x$ncpgs <- tmp$ncpgs
#    x$locs <- tmp$locs
#    x$meth <- tmp$meth
#    x
#})
#
#readCpgs <- lapply(readGRs, function(gr) {
#    keep <- gr$ncpgs > 0
#    gr <- gr[keep,]
#    locs <- rep(start(gr), gr$ncpgs) + unlist(gr$locs)
#
#    Cov <- rep(1, length(locs))
#    Meth <- 1*(unlist(gr$meth) == "M")
#
#    Cov <- tapply(Cov, locs, sum)
#    Meth <- tapply(Meth, locs, sum)
#    Loc <- tapply(locs, locs, function(x) x[1])
#    chr <- tapply(rep(as.character(seqnames(gr)), gr$ncpgs), locs, function(x) x[1])
#    newGR <- GRanges(chr, IRanges(start=Loc, width=1), Cov=Cov, Meth=Meth)
#    newGR$Beta <- Meth / Cov
#    names(newGR) <- NULL
#    newGR
#})
#
#readCpgs <- lapply(readCpgs, renameSeqlevels, value=structure("chr13", names="13"))
#
betaGR_cpg <- lapply(seq(along=rawCpgs_cpg), function(i) {
    rawGR <- rawCpgs_cpg[[i]]
    estGR <- estimatedCpgs_cpg[[i]]
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

pdf(file.path(figdir_cpg,"meth_percentage_1.pdf"), height=7, width=9)
mypar(2,3)

for (i in seq(along=betaGR_cpg)) {
    gr <- betaGR_cpg[[i]]
    plot(gr$rawBeta,gr$estimatedBeta,
         bty='l',
         main=names(filteredObjs_cpg)[i],
         ylab="pattern methyl Percentage",
         xlab="region methyl Precentage",
         cex = 1,
         cex.lab=2)
}
dev.off()

pdf(file.path(figdir_region,"meth_percentage_1.pdf"), height=7, width=9)
mypar(2,3)

for (i in seq(along=betaGR_region)) {
    gr <- betaGR_region[[i]]
    plot(gr$rawBeta,gr$estimatedBeta,
    bty='l',
    main=names(filteredObjs_region)[i],
    ylab="pattern methyl Percentage",
    xlab="region methyl Precentage",
    cex = 1,
    cex.lab=2)
}
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



pdf(file.path(figdir_cpg,"meth_pct_2.pdf"), height=9, width=10)
mypar(2,3)

#for (i in c(1,2,3)) {
for (i in seq(along=betaGR_cpg)) {
    gr_cpg <- betaGR_cpg[[i]]
    plot(gr_cpg$rawBeta,gr_cpg$estimatedBeta,
    bty='l',
    main=pd$dirname[i],
    ylab="pattern methyl Percentage",
    xlab="region methyl Precentage",
    col="blue",
    cex = 1,
    cex.lab=3)
    
    gr_region <- betaGR_region[[i]]
    points(gr_region$rawBeta,gr_region$estimatedBeta,
    bty='l',
    col="red")
}

legend("bottomright", c("cpg_loss","region_loss"), col=c("blue", "red"), pch=1)


dev.off()

err_cpg = c()
for (i in seq(along=betaGR_cpg)) {
    gr <- betaGR_cpg[[i]]
    err_cpg[i] <- sqrt(sum(gr$rawBeta - gr$estimatedBeta)^2/length(betaGR_cpg[[i]]))
}

err_region = c()
for (i in seq(along=betaGR_region)) {
    gr <- betaGR_region[[i]]
    err_region[i] <- sqrt(sum(gr$rawBeta - gr$estimatedBeta)^2/length(betaGR_region[[i]]))
}

pd$error_cpg <- err_cpg
pd$error_region <- err_region

error <- rbind(err_cpg, err_region)

colnames(error)= pd$dirname


pdf(file.path(figdir_cpg,"meth_sqr_err_loss.pdf"))
barplot(error,legend=T, beside=T)
dev.off()





allBeta <- betaGR_region[[1]]
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
clusters <- clusterMaker(seqnames(allBeta), start(allBeta), maxGap=100)

#install.packages("doParallel", repos="http://R-Forge.R-project.org")
library(doParallel)
registerDoParallel(cores=2)
betamat <- as.matrix(mcols(allBeta))

bumps <- bumphunter(betamat, model.matrix(~1+status,data=pd),
                    pos=start(allBeta), cluster=clusters, B=0, smooth=FALSE, cutoff=.2)



plotRegion <- GRanges("13", IRanges(start=bumps$tab$start[1], bumps$tab$end[1]))
keep <- components(filteredObjs_region[[6]]) %over% plotRegion
plotRegion <- plotRegion * -8

#plotRegion <- GRanges("13", IRanges(start=111128000,end=111128200))

regionBeta <- subsetByOverlaps(allBeta, plotRegion)


betamat <- as.matrix(mcols(regionBeta))
Index <- split(pd$status, pd$status)
betamns <- sapply(Index, function(ind) rowMeans(betamat[,ind]))

fits <- lapply(1:2, function(i) loess(betamns[,i]~start(regionBeta)))




load_all(path)
pdf(file.path(figdir_region, "region.pdf"),width=8,height=6)

layoutmat <- matrix(c(8,rep(1,4),9,rep(2,3),rep(3,3),rep(4,3),rep(5,3),rep(6,3),rep(7,3)),nc=1)
layout(layoutmat)
par(mar=c(1, 5, 0, 1))
par(oma=c(0,0,0,0))
    
xlim <- c(start(plotRegion), end(plotRegion)+40)

plot(0, type="n", xlim=xlim, ylim=c(0,1), ylab="",xaxt="n",yaxt="n")
lat <- pretty(xlim-50,n=4)
axis(side=1, at=lat, cex.axis=1, mgp=c(1,.3,0))
axis(side=2, at=c(0, 0.5, 1),las=1,cex.axis=1, mgp=c(1,.4,0))
title(xlab="genomic position (chr13)", ylab="%methyl", cex.lab=1.2)

for (i in seq(along=filteredObjs_region)) {
    points(start(regionBeta), mcols(regionBeta)[[names(filteredObjs_region[i])]], col=ifelse(pd$status[i]=="N",1,2), cex=.5, pch=19)
}

for (i in 1:2) {
    curve(predict(fits[[i]], x), from=min(start(regionBeta)), to=max(start(regionBeta)), col=i,add=TRUE,lwd=1.5)
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
  plot(0, type="n", xlim=c(start(plotRegion) , end(plotRegion)+40), ylim=c(0,npatterns+1),  xaxt="n", yaxt="n", ylab="", xlab="",
       bty="n")
  
  for (j in seq(along=cids)) {
    cid <- cids[j]
    patterns <- patterns(filteredObjs_region[[i]])[patterns(filteredObjs_region[[i]])$cid == cid,]
    
    o <- order(patterns$abundance)
    patterns <- patterns[o,]
    
    npatterns <- length(patterns)
    
    segments(min(start(patterns)), seq(along=patterns),
             max(end(patterns))-10, seq(along=patterns))
    
    
    
    #        text(end(patterns), seq(along=patterns),
    #            labels=sprintf("%.3f", patterns$abundance), pos=4,cex=.8)
    
    scaledAbundance <- patterns$abundance / max(patterns$abundance)
    quantizedAbundance <- cut(scaledAbundance, breaks=seq(0,1,len=11))
    
    cols <- colorRampPalette(brewer.pal(5,"Blues"))(11)
    col <- cols[quantizedAbundance]
    
    rect(max(end(patterns))-3,seq(along=patterns)-.5,max(end(patterns)) + 12,seq(along=patterns)+.5,
         col=col)
    text(max(end(patterns))+15, (npatterns+1)/2, label="abundance",cex=1,srt=90)
    #        points(end(patterns)+10, seq(along=patterns), pch=22, bg=col)
    text(max(end(patterns))+17, (npatterns+1)/2, label= names(filteredObjs_region[i]) ,cex=1,srt=90)
    
    locs <- rep(start(patterns), patterns$ncpgs) +
      unlist(patterns$locs)
    ylocs <- rep(seq(along=patterns), patterns$ncpgs)
    methyl <- unlist(patterns$meth)
    col  <- ifelse(methyl=="M", "black", "white")
    points(locs,ylocs,pch=21,bg=col)
  }
  
  

}
dev.off()


load_all(path)
entrStats <- lapply(filteredObjs_region, getEntropyStats)
jointGR <- disjoin(unlist(do.call("GRangesList", entrStats)))

makeVec <- function(gr, jointGR, slot) {
    olaps <- findOverlaps(jointGR, gr)
    out <- rep(NA,length(jointGR))
    out[queryHits(olaps)] <- mcols(gr)[[slot]][subjectHits(olaps)]
    out
}
entrMat <- sapply(entrStats, makeVec, jointGR, "entr")
nentrMat <- sapply(entrStats, makeVec, jointGR, "normEntr")
giniMat <- sapply(entrStats, makeVec, jointGR, "gini")
maxEntrMat <- sapply(entrStats, makeVec, jointGR, "maxEntr")

colData=DataFrame(pd)
rownames(colData)=pd$dirname
statsSE <- SummarizedExperiment(rowData=jointGR,
                                colData=colData,
                                assays=SimpleList(entropy=entrMat,
                                    normEntropy=nentrMat,
                                    maxEntropy=maxEntrMat,
                                    gini=giniMat))

pdf(file.path(figdir_region, "entropy.pdf"),width=9,height=3)
mypar(1,3)

for (i in 1:3) {
    mat = assay(statsSE,"normEntropy")[,c(i+3,i)]
    plot(mat,main="normalized entropy" )
    abline(0,1)
}

for (i in 1:3) {
    mat = assay(statsSE,"normEntropy")[,c(i+3,i)]
    d = mat[,2]-mat[,1]
    hist(d,nc=50,xlab=sprintf("%s-%s", colnames(mat)[2], colnames(mat)[1]),main="normalized entropy")
}

for (i in 1:3) {
    mat = assay(statsSE,"maxEntropy")[,c(i+3,i)]
    plot(mat, main="log num patterns")
    abline(0,1)
}

for (i in 1:3) {
    mat = assay(statsSE,"maxEntropy")[,c(i+3,i)]
    d = mat[,2]-mat[,1]
    hist(d,nc=50,xlab=sprintf("%s-%s", colnames(mat)[2], colnames(mat)[1]),main="log num patterns")
}
dev.off()

jointGR <- disjoin(unlist(do.call(GRangesList, betaGR_region)))
betaMat <- sapply(betaGR_region, makeVec, jointGR, "estimatedBeta")
normEntrMat <- sapply(entrStats, makeVec, jointGR, "normEntr")
entrMat <- sapply(entrStats, makeVec, jointGR, "entr")
cpgSE <- SummarizedExperiment(rowData=jointGR,
                              colData=colData,
                              assays=SimpleList(beta=betaMat,
                                  entropy=entrMat,
                                  normEntropy=normEntrMat))

pdf(file.path(figdir_region,"methEntropy.pdf"),width=9,height=3)
mypar(1,3)

for (i in 1:3) {
    mat = assay(cpgSE, "entropy")[,c(i+3,i)]
    plot(mat,main="entropy")
}

for (i in 1:3) {
    mat = assay(cpgSE,"entropy")[,c(i+3,i)]
    d = mat[,2]-mat[,1]
    xlab="tumor - normal entropy"
    main=paste("sample", gsub("CAP_[N|T]_","",colnames(mat)[1]))
    hist(d,nc=50,xlab=xlab,main=main)
}

for (i in 1:3) {
    mat = assay(cpgSE, "beta")[,c(i+3,i)]
    plot(mat,main="methylation")
}

for (i in 1:3) {
    mat = assay(cpgSE,"beta")[,c(i+3,i)]
    d = mat[,2]-mat[,1]
    hist(d,nc=50,xlab=sprintf("%s-%s", colnames(mat)[2], colnames(mat)[1]),main="methylation")
}

for (i in 1:6) {
    plot(assay(cpgSE,"beta")[,i], assay(cpgSE,"entropy")[,i],main=pd$dirname[i])
}
dev.off()






################

compflow <-  tapply(patterns(obj)$abundance, patterns(obj)$cid, sum)
patternp <- patterns(obj)$abundance / compflow[as.character(patterns(obj)$cid)]

tmp <- lapply(4:6, function(subj) {
    nindex <- match(paste("CAP","N",subj,sep="_"), names(objs))
    tindex <- match(paste("CAP","T",subj,sep="_"), names(objs))
    nobj <- objs[[nindex]]
    tobj <- objs[[tindex]]
    
    ncomps <- components(nobj)
    tcomps <- components(tobj)
    
    d <- distanceToNearest(ncomps, tcomps)
    keep <- mcols(d)$distance == 0 & width(ncomps)>100 & width(tcomps[subjectHits(d),]) > 100
    iind <- which(keep)
    jind <- subjectHits(d)[keep]
    
    list(comps=ncomps[keep,],
         npats=cbind(npats[[nindex]][iind], npats[[tindex]][jind]),
         ent=cbind(compent[[nindex]][iind], compent[[tindex]][jind]),
         meth=cbind(compMeth[[nindex]][iind], compMeth[[tindex]][jind]))
})

pdf(file.path(figdir,"compcompare.pdf"),height=4,width=12)
mypar(1,3)

  for (i in 1:3) {
    rng <- range(tmp[[i]]$npats, na.rm=TRUE)
    plot(tmp[[i]]$npats, xlim=rng, ylim=rng, xlab="normal num. patterns", ylab="tumor num. patterns",cex.lab=1.7,cex=1.3,pch=19)
    abline(0,1)
    
    rng <- range(tmp[[i]]$ent, na.rm=TRUE)
    plot(tmp[[i]]$ent, xlim=rng, ylim=rng, xlab="normal entropy", ylab="tumor entropy",cex.lab=1.7,cex=1.3,pch=19)
    abline(0,1)
    
    nent <- tmp[[i]]$ent / log(tmp[[i]]$npats)
    rng <- range(nent,na.rm=TRUE)
    plot(nent, xlim=rng, ylim=rng, xlab="normal entropy (normalized)", ylab="tumor entropy (normalized",cex.lab=1.7,cex=1.3,pch=19)
    abline(0,1)
}
dev.off()
  
pdf(file.path(figdir, "compcompare2.pdf"), height=4, width=12)
mypar(1,3)

  for (i in 1:3) {
    nent <- tmp[[i]]$ent / log(tmp[[i]]$npats)
    deltam <- tmp[[i]]$meth[,2]-tmp[[i]]$meth[,1]
    deltae <- nent[,2] - nent[,1]
    plot(deltam, deltae,xlab="methylation difference", ylab="entropy difference",cex.lab=1.7,cex=1.3,pch=19,main=paste("Subject",i))
  }
  dev.off()
  
})




plotRegion2 <- GRanges("chr1",IRanges(start=102600000,end=102700000))
plotComps <- lapply(objs, function(obj) which(components(obj) %over% plotRegion2))

objs2 <- lapply(objs, processMethylpats)

compMeth <- lapply(objs2, componentAvgMeth)

require(epivizr)
mgr <- startEpiviz()

devs <- lapply(seq(along=objs), function(i) mgr$addDevice(components(objs[[i]]), names(objs)[i]))
