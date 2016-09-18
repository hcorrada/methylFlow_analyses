#!/usr/bin/env Rscript
### run ./sam.r par1

### par1 = input directory
#data[1] = "SRR1015705-10006910-15008000"
library("devtools")
path = "~/Projects/methylFlow/src/methylFlow/methylFlowr"
load_all(path)



dirNew <- "/Users/faezeh/Projects/methylFLow"
dirOld <- "/Users/faezeh/Projects/methylFLow_old"
#install_github("rafalib","ririzarr")
library("rafalib")

resultdir <- file.path(paste(dirNew,"/exps/sam",sep=""))
#resultdirOld <- file.path(dirOld,"exps/sam/auto")



#dir <- dirOld
#resultdir <- resultdirOld

dir <- dirNew
#resultdir <- resultdirNew


print(resultdir)

pd <- data.frame(subject=c("1020523", "1020509", "1015705", "1015434"), status=c("wild B-cell","wild B-cell", "KSL", "KSL"))
j <- "3"

pd$dirname=sprintf("SRR%s-bam", pd$subject)
pd$name=sprintf("SRR%s", pd$subject)
#pd$dirname = datadir
objs <- lapply(seq(len=nrow(pd)), function(i) {
    curdir <- file.path(resultdir,pd$dirname[i])
    read.methylflow.dir(curdir,pd$dirname[i], has.header=TRUE)
})
names(objs) <- pd$name
#names(objs)<-pd$status
objs <- lapply(objs, processMethylpats)
filteredObjs <- lapply(objs, mfFilterBy, minComponentCoverage=50, minComponentWidth=86)
names(filteredObjs) <- names(objs)


figdir <- file.path(resultdir,"B-cell-new")

if (!file.exists(figdir)) dir.create(figdir)

pdf(file.path(figdir,"fragment_length.pdf"),height=9,width=9)
#mypar(1,1)

widths <- lapply(filteredObjs, function(obj) width(components(obj)))
ncpgs <- lapply(filteredObjs, function(obj) ncpgs(obj, level="pattern", summary="max"))


par(cex.axis=1.5)
par(cex.lab=1.5)
names(ncpgs) <- names(widths) <- names(filteredObjs)
boxplot(widths,  cex=1.5)
boxplot(ncpgs, cex=1.5)
dev.off()



####
filteredObjs <- lapply(filteredObjs, mfFilterBy, minNumberOfPatterns=1)

npats <- lapply(filteredObjs, npatterns, by.component=TRUE)
compCoverage <- lapply(filteredObjs, counts, level="component", kind="raw")

pdf(file.path(figdir,"patt_by_coverage.pdf"), height=2, width=8)
mypar(1,4)

for (i in seq(along=filteredObjs)) {
    plot(compCoverage[[i]], npats[[i]], main=names(filteredObjs[i]),xlab="coverage",ylab="num. patterns", cex=0.7)
}
dev.off()



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

#pdf(file.path(figdir,"meth_percentage.pdf"), height=9, width=9)
#mypar(1,1)

for (i in seq(along=betaGR)) {
  pdf(file.path(figdir,paste(i,"_meth_percentage.pdf")), height=8, width=8)
    gr <- betaGR[[i]]
    par(mar=c(5, 6, 3, 2) + 0.1)
    plot(gr$rawBeta,gr$estimatedBeta,
    bty='l',
    main=names(filteredObjs)[i],
    ylab="pattern methyl Percentage",
    xlab="region methyl Precentage",
    cex.axis = 1.5,
    cex.lab=2)
  dev.off()
}

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
#dev.off()


allBeta <- betaGR[[1]]
tmp <- allBeta$estimatedBeta

mcols(allBeta) <- NULL
mcols(allBeta)[[names(filteredObjs)[1]]] <- tmp

for (i in seq(2,length(betaGR))) {
  #i=2
  gr <- betaGR[[i]]
  olaps <- findOverlaps(gr, allBeta)
  
  sname <- names(filteredObjs)[i]
  mcols(allBeta)[[sname]] <- NA
  mcols(allBeta)[[sname]][subjectHits(olaps)] <- gr$estimatedBeta[queryHits(olaps)]
}

library(bumphunter)
clusters <- clusterMaker(seqnames(allBeta), start(allBeta), maxGap=100)

library(doParallel)
registerDoParallel(cores=2)
betamat <- as.matrix(mcols(allBeta))

bumps <- bumphunter(betamat, model.matrix(~1+status,data=pd),
                    pos=start(allBeta), cluster=clusters, B=0, smooth=FALSE, cutoff=.2)

plotRegion <- GRanges("3", IRanges(start=bumps$tab$start[1], bumps$tab$end[1]))
plotRegion <- plotRegion * -8

#plotRegion <- GRanges("3", IRanges(start=38444300,end=38444850))

regionBeta <- subsetByOverlaps(allBeta, plotRegion)

betamat <- as.matrix(mcols(regionBeta))
Index <- splitit(pd$status)
betamns <- sapply(Index, function(ind) rowMeans(betamat[,ind]))

fits <- lapply(1:2, function(i) loess(betamns[,i]~start(regionBeta)))

load_all(path)
pdf(file.path(figdir, "region_sum_v3.pdf"),width=8,height=6)

layoutmat <- matrix(c(9,rep(1,4),10,rep(2,3),rep(3,3),rep(4,3),rep(5,3),rep(6,3),rep(7,3),rep(8,1)),nc=1)
layout(layoutmat)
par(mar=c(1, 5, 0, 1))
par(oma=c(0,0,0,0))

xlim <- c(start(plotRegion), end(plotRegion)-12)
#xlim <- c(min(start(regionBeta)), max(end(regionBeta))+40)


plot(0, type="n", xlim=xlim, ylim=c(0,1), ylab="",xaxt="n",yaxt="n")
lat <- pretty(xlim-20,n=4)
axis(side=1, at=lat, cex.axis=1, mgp=c(1,.3,0))
axis(side=2, at=c(0, 0.5, 1),las=1,cex.axis=1, mgp=c(1,.4,0))
title(xlab="genomic position (chr3)", ylab="%methyl", cex.lab=1.2)


for (i in seq(along=filteredObjs)) {
  points(start(regionBeta), mcols(regionBeta)[[names(filteredObjs[i])]], col=ifelse(pd$status[i]=="wild", 1, 2), cex=.5, pch=19)
}

for (i in 1:2) {
  curve(predict(fits[[i]], x), from=min(start(regionBeta)), to=max(start(regionBeta)), col=3 - i ,add=TRUE,lwd=1.5)
}
legend("bottomright", c("Wild","Cancer"), pch=22,pt.bg=1:2)

indexes <- c(1,3,2,4)
for (i in indexes) {
 # plotPatterns(filteredObjs[[i]], plotRegion,
  #             xaxt="n", yaxt="n", ylab="", xlab="",
   #            bty="n")
  #    mtext(side=3,paste(gsub("CAP_","",names(filteredObjs)[i])), cex=.7)
  keep <- components(filteredObjs[[i]]) %over% plotRegion
  
  npats <- npatterns(filteredObjs[[i]], by.component=TRUE)[keep]
  cids <- components(filteredObjs[[i]])$cid[keep]
  
  npatterns <- max(npats)
  plot(0, type="n", xlim = xlim, ylim=c(0,npatterns+1),  xaxt="n", yaxt="n", ylab="", xlab="",
       bty="n")
  
  for (j in seq(along=cids)) {
    cid <- cids[j]
    patterns <- patterns(filteredObjs[[i]])[patterns(filteredObjs[[i]])$cid == cid,]
    
    o <- order(patterns$abundance)
    patterns <- patterns[o,]
    
    npatterns <- length(patterns)
    
   # segments(min(start(patterns)), seq(along=patterns),
    #         max(end(patterns)), seq(along=patterns))
    segments(min(start(plotRegion))-1, seq(along=patterns),
             max(end(regionBeta))+2, seq(along=patterns))
    
    
    #        text(end(patterns), seq(along=patterns),
    #            labels=sprintf("%.3f", patterns$abundance), pos=4,cex=.8)
    
#    scaledAbundance <- patterns$abundance / max(patterns$abundance)
    scaledAbundance <- patterns$abundance / sum(patterns$abundance)
   
    quantizedAbundance <- cut(scaledAbundance, breaks=seq(0,1,len=11))
    
    cols <- colorRampPalette(brewer.pal(5,"Blues"))(11)
    cols <- colorRampPalette(c("#EFF3FF", "#2163A7", "#08519C"))(10)

    col <- cols[quantizedAbundance]
    
   # rect(max(end(patterns))+10,seq(along=patterns)-.5,max(end(patterns)) + 40,seq(along=patterns)+.5,
    #     col=col)
   
    rect(max(end(regionBeta))+4,seq(along=patterns)-.5,max(end(regionBeta)) + 8,seq(along=patterns)+.5,
        col=col)
   # text(max(end(patterns))+42, (npatterns+1)/2, label="abundance of ", cex=.7,srt=90)
  #  text(max(end(patterns))+47, (npatterns+1)/2, label= names(filteredObjs[i]) ,cex=.7,srt=90)
   
   #text(max(end(regionBeta))+10, (npatterns+1)/2, label="abundance", cex=1,srt=90)
   #   text(max(end(regionBeta))+12, (npatterns+1)/2, label= names(filteredObjs[i]) ,cex=0.9,srt=90)
   text(max(end(regionBeta))+10, (npatterns+1)/2, label= pd$status[i] ,cex=0.7,srt=90)

    #        points(end(patterns)+10, seq(along=patterns), pch=22, bg=col)
    
    locs_all <- rep(start(patterns), patterns$ncpgs) + unlist(patterns$locs)
  
    locs <- locs_all
    locs <- locs_all[locs_all > start(plotRegion) & locs_all < end(plotRegion)]
    
    nlocs <- patterns$ncpgs
    nlocs <- rep(length(locs) / npatterns, npatterns)
  
    ylocs <- rep(seq(along=patterns), nlocs)
    methyl <- unlist(patterns$meth)
    col  <- ifelse(methyl=="M", "black", "white")
    points(locs,ylocs,pch=21,bg=col)
  }
  
  
  
}


#mypar(1,1)
p <- seq(0,1,len=11)
#keyCols <- colorRampPalette(brewer.pal(5,"Blues"))(10)
keyCols <- colorRampPalette(c("#EFF3FF", "#2163A7", "#08519C"))(10)
#keyCols <- rgb(grayrgb[1],grayrgb[2],grayrgb[3],alpha=255*p,max=255)
plot(1,1,type="n",xlim=c(-0.5,length(keyCols)),ylim=c(0,1),bty="n",
     xaxt="n",yaxt="n",xlab="",ylab="")
#rect(xleft=0,xright=1,
 #    ybottom=seq(along=keyCols)-1,
  #   ytop=seq(along=keyCols),col=keyCols)
rect(xleft=seq(along=keyCols)/2+1,xright=seq(along=keyCols)/2+1.5,
         ybottom=0.1,
        ytop=0.4,col=keyCols)
     
     
axis(side=1,at=p*5+1.5,p,las=2, cex=0.8)


dev.off()

