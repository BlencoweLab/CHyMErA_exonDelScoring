### Gonatopoulos-Purnatzis, Aregger, Brown et al.
### Genetic interaction mapping and exon-resolution functional genomics with a hybrid Cas9–Cas12a platform
### Nature Biotechnology (2020)
### https://doi.org/10.1038/s41587-020-0437-z
###
### CHyMErA screen for growth phenotypes of alternative exon deletion in RPE1 cells
### Scoring of growth phenotypes per exon based on drop-out of guide pairs


### Change to fit your local paths ###
countsFile  <- "input/exonDeletionLibrary_normCounts_NovaSeq_18Sept18.txt.gz"
scoresFile1 <- "input/intronic_intronic_guides.tsv.gz"
scoresFile2 <- "input/other_guides.tsv.gz"
annoFile    <- "input/Knockout.exons.human.csv"
vastToolsExprFile <- "input/cRPKM_AND_COUNTS-Hsa1.tab.gz"
vastToolsASFile   <- "input/INCLUSION_LEVELS_FULL-Hsa1-hg19.tab.gz"

outPath <- "output"


### HIT CUTOFF ###
hitCO <- 0.18  # Minimum fraction of exon-deletion guide paris that must be significant for an exon
               # to be deemed to have a phenotype
minClean <- 0  # Minimum number of significant guide pairs of which neither single guide, when paired
               # with an intergenic control guide, drops out significantly on its own.


### Function definitions ###

generateAnnotation <- function(counts,           # Read counts table from the CHyMErA screen at hand
                               annoFile,         # Basic exon annotation from design phase
                               vastToolsASFile,  # Vast-tools AS output for given cell
                               vastToolsExprFile,# Vast-tools expression output for given cell
                               outPath,          # Output Path
                               lfcCO=c(-3,-1),   # Log2-fold change cutoffs definiting essential and non-essential
                                                 # Dropout of TKOv3 Cas9 guides to assess gene depletion
                                                 # Essentiality within-screen
                               exprLo=5,         # Minimum expression RPKM to be considered expressed
                               PSIskipped=5,     # PSI below which an exon is considered skipped
                               PSIlow=15         # PSI below which an exon is considered low
                               ) {
### Supplement the exon annotation with gene essentiality, expression and exon splicing status in the cell line
### Value: data.frame with exon annotation
    
    anno   <- read.csv(annoFile)
    anno   <- anno[anno$EVENT %in% counts$Event,]
    
    ## Determine fitness effect based on TKOv3 guides targeting exons
    geneLFC <- aggregate(counts$RPE1.logFC_T24[counts$Data.Subset == "ExonDeletion_TKOv3"],
                         by=list(gene=counts$Target.gene[counts$Data.Subset == "ExonDeletion_TKOv3"]),
                         FUN=mean, na.rm=T)
    
    pdf(file.path(outPath, "Fitness.RPE1.pdf"), wid=4, hei=4.4)
    plot(dens <- density(geneLFC$x, na.rm=T), type="n",
         xlab="Mean guide depletion log2 (fold change) T24", main="Fitness effect in RPE1")
    lo <- which(dens$x <= lfcCO[1])
    polygon(c(dens$x[lo], rev(dens$x[lo])), c(dens$y[lo], rep(0, length(lo))), col="red", border="red")
    hi <- which(dens$x >= lfcCO[2])
    polygon(c(dens$x[hi], rev(dens$x[hi])), c(dens$y[hi], rep(0, length(hi))), col="dodgerblue1",
            border="dodgerblue1")
    lines(dens, lwd=2)
    abline(v=lfcCO, lty=2)
    legend("topright", legend=c("Essential","Non-essential"), text.col=c("red","dodgerblue1"), bty="n")
    dev.off()
    
    
    ## Load vast-tools AS and expression output
    as <- read.delim(vastToolsASFile)
    as$RPE1.Cas9[!cleanAS(as$RPE1.Cas9.Q, as$COMPLEX)] <- NA
    
    expr <- read.delim(vastToolsExprFile)

    anno <- data.frame(anno[,1:5],
                       PSI.RPE1        = as$RPE1.Cas9[match(anno$EVENT, as$EVENT)],
                       cRPKM.RPE1      = expr$RPE1.Cas9.cRPKM[match(anno$EnsemblID, expr$ID)],
                       ess.RPE1.log2FC = geneLFC$x[match(anno$GENE, geneLFC$gene)],
                       anno[, c("protImpact","selection","cyclic")]
                       )
    
    
    ## Check how expression and PSI agree with essentiality
    pdf(file.path(outPath, "Essentiality_RNAseq.pdf"), wid=5, hei=5.4)
    plot(anno$ess.RPE1.log2FC, log10(0.1 + anno$cRPKM.RPE1), main="TKO v3 essentiality and RNA-seq expr.",
         xlab="TKOv3 log2 fold change", ylab="Expression log10 (0.1 + cRPKM)")
    abline(v=c(-3, -1), lty=2, col=c("red","dodgerblue1"))
    
    plot(anno$ess.RPE1.log2FC, anno$PSI.RPE1, main="TKO v3 essentiality and RNA-seq AS.",
         xlab="TKOv3 log2 fold change", ylab="PSI")
    abline(v=c(-3, -1), lty=2, col=c("red","dodgerblue1"))
    dev.off()
    
    
    ## Add cell-type specific annotation
    tmp <- data.frame(exprOK    = anno$cRPKM.RPE1 > exprLo,
                      frameDisr = anno$protImpact == "ORF disruption upon sequence exclusion",
                      ess       = anno$ess.RPE1.log2FC <= lfcCO[1],
                      nonEss    = anno$ess.RPE1.log2FC >= lfcCO[2],
                      undec     = anno$ess.RPE1.log2FC >  lfcCO[1] & anno$ess.RPE1.log2FC < lfcCO[2] |
                          is.na(anno$ess.RPE1.log2FC)
                      )
    
    sel  <- rep(NA, nrow(anno))
    sel[anno$cRPKM.RPE1 < 1] <- "notExpressed"
    sel[tmp$exprOK & tmp$ess    &  tmp$frameDisr] <- "frameDisr/essential"
    sel[tmp$exprOK & tmp$nonEss &  tmp$frameDisr] <- "frameDisr/nonessential"
    sel[tmp$exprOK & tmp$undec  &  tmp$frameDisr] <- "frameDisr"
    sel[tmp$exprOK & tmp$ess    & !tmp$frameDisr] <- "framePres/essential"
    sel[tmp$exprOK & tmp$nonEss & !tmp$frameDisr] <- "framePres/nonessential"
    sel[tmp$exprOK & tmp$undec  & !tmp$frameDisr] <- "framePres"
    sel[tmp$exprOK & anno$PSI.RPE1 < PSIskipped & (anno$selection != "other" | anno$cyclic)] <- "ctrl-skipped"
    sel[(!tmp$exprOK | anno$PSI.RPE1 < PSIlow)  & is.na(sel)] <- "lowExpressed/lowPSI"
    sel[is.na(sel)] <- "other"


    pdf(file.path(outPath, "Knockout.exons_RPE1.pdf"), wid=5, hei=4.5)
    dat <- table(sel)[c(3,4,2,6,7,5,1,9,8,10)]
    par(mar=c(5,11,4,1))
    ypos <- barplot(rev(dat), horiz=T, las=1,
                    col=c(rep("white",2),rep("dodgerblue1",2),rep("grey70",5),rep("brown1")),
                    xlab="Number of exons", main="Interrogated exons per type\nin RPE1")
    abline(v=0)
    text(0, ypos, rev(dat), pos=4)
    dev.off()
    
    data.frame(anno, type.RPE1 = sel)
}


cleanAS <- function(x,                                        # Quality column of one sample
                    complexCol,                               # Colum 'COMPLEX' of vast-tools output
                    covThresCE    = c("SOK","OK","LOW"),      # Allowable coverage score (score 3) for CE
                    covThresOther = c("SOK","OK","LOW"),      # Allowable coverage score (score 3) for ME, Alt3, Alt5
                    covThresIR    = 15,                       # Minimum coverage (score 3) for IR
                    balThresCE    = c("OK","B1"),             # Allowable balance score (score 4) for CE
                    balThresIR    = 0.05                      # Maximum balance p-value (score 4, bigger is better) for IR
                    ) {
### Filter PSI output of vast-tools AS pipeline using the quality column
### See https://github.com/vastgroup/vast-tools
### Value: logical vector indicating which events passed the criteria
    
    is.ce <- complexCol %in% c("S","C1","C2","C3","ANN")
    is.ir <- grepl("IR", complexCol)
    is.o  <- complexCol %in% c("Alt5", "Alt3", "MIC")

    cov <- sub("[^,]+,[^,]+,([^,]+),.+", "\\1", x)
    bal4 <- sub("[^,]+,[^,]+,[^,]+,([^,]+),.+", "\\1", x)
    bal5 <- sub("[^,]*,[^,]*,[^,]*,[^,]*,([^,@]*)[@,].*", "\\1", x)
    suppressWarnings(
        nreads <- as.integer(sub("([^=]*)=.+", "\\1", bal4)) +
                  as.integer(sub("[^=]*=([^=]*)=.*", "\\1", bal4)) +
                  as.integer(sub("[^=]*=[^=]*=(.*)", "\\1", bal4))
    )
    ok.ce <- is.ce & cov %in% covThresCE & bal4 %in% balThresCE
    suppressWarnings(
        ok.ir <- is.ir & nreads >= covThresIR & as.numeric(bal5) > balThresIR
    )
    ok.o  <- is.o  & cov %in% covThresOther 

    ok.ce | ok.ir | ok.o
}





### Read inputs ###

if (!dir.exists(outPath)) {dir.create(outPath)}

scores1 <- read.delim(scoresFile1)
scores2 <- read.delim(scoresFile2)
scores  <- rbind(scores1[,match(names(scores2), names(scores1))], scores2)

counts  <- read.delim(countsFile)


### Supplement exon annotation with cell type specific expression, PSI,
### and gene essentiality
anno <- generateAnnotation(counts            = counts,
                           annoFile          = annoFile,
                           vastToolsASFile   = vastToolsASFile,
                           vastToolsExprFile = vastToolsExprFile,
                           outPath           = outPath
                           )
                           

### Count significantly depleted guides per exon ###

counts <- counts[match(scores$ID, counts$ID),]

## Remove significant pairs that are increasing during screen
scores$rpe1_sig[scores$RPE1_logFC > 0] <- FALSE


conformCount <- function(set=NA) {
    x <- aggregate(scores$Event[set], by=list(scores$Event[set]), FUN=length)
    out <- x$x[match(anno$EVENT, x[,1])]
    out[is.na(out)] <- 0
    out
}

exscr <- data.frame(Event = anno$EVENT)
exscr$n.exonDel <- conformCount(scores$type == "exon")
exscr$sig.exonDel <- conformCount(scores$type == "exon" & scores$rpe1_sig)
exscr$n.singleCtrl <- conformCount(scores$type == "single")
exscr$sig.singleCtrl <- conformCount(scores$type == "single" & scores$rpe1_sig)
exscr$n.exonic <- conformCount(scores$type == "exonic")
exscr$sig.exonic <- conformCount(scores$type == "exonic" & scores$rpe1_sig)

exscr$n.Cas9 <- { # number of Cas9 guides in exon-deletion pairs
    sel <- grep("intronic",scores$Cas9_Guide_Type)
    x <- aggregate(counts$Cas9.Guide[sel], by=list(scores$Event[sel]),
                   FUN=function(x) {length(unique(x))})
    out <- x$x[match(anno$EVENT, x[,1])]
    out[is.na(out)] <- 0
    out
}

exscr$n.Cpf1 <- { # number of Cpf1 guides in exon-deletion pairs
    sel <- grep("intronic",scores$Cpf1_Guide_Type)
    x <- aggregate(counts$Cpf1.Guide[sel], by=list(scores$Event[sel]),
                   FUN=function(x) {length(unique(x))})
    out <- x$x[match(anno$EVENT, x[,1])]
    out[is.na(out)] <- 0
    out
}

exscr$sig.Cas9.single <- { # Cas9 guides with significant dropout when paired with an intergenic control
    sel <- which(grepl("intronic",scores$Cas9_Guide_Type) & grepl("intergenic",scores$Cpf1_Guide_Type) &
                 scores$rpe1_sig)
    x <- aggregate(counts$Cas9.Guide[sel], by=list(scores$Event[sel]),
                   FUN=function(x) {length(unique(x))})
    out <- x$x[match(anno$EVENT, x[,1])]
    out[is.na(out)] <- 0
    out
}

exscr$sig.Cpf1.single <- { # Cpf1 guides with significant dropout when paired with an intergenic control
    sel <- which(grepl("intergenic",scores$Cas9_Guide_Type) & grepl("intronic",scores$Cpf1_Guide_Type) &
                 scores$rpe1_sig)
    x <- aggregate(counts$Cpf1.Guide[sel], by=list(scores$Event[sel]),
                   FUN=function(x) {length(unique(x))})
    out <- x$x[match(anno$EVENT, x[,1])]
    out[is.na(out)] <- 0
    out
}

exscr$exonDel.clean <- {     # Number of deletion guide pairs, neither of which is significant
                             # on its own when paired with an intergenic guide
    bad.Cas9 <- grepl("intronic",scores$Cas9_Guide_Type) & grepl("intergenic",scores$Cpf1_Guide_Type) &
        scores$rpe1_sig
    bad.Cpf1 <- grepl("intergenic",scores$Cas9_Guide_Type) & grepl("intronic",scores$Cpf1_Guide_Type) &
        scores$rpe1_sig
    sel <- which(grepl("intronic",scores$Cas9_Guide_Type) & grepl("intronic",scores$Cpf1_Guide_Type) &
                 !(counts$Cas9.Guide %in% counts$Cas9.Guide[bad.Cas9]) &
                 !(counts$Cpf1.Guide %in% counts$Cpf1.Guide[bad.Cpf1]))
    x <- aggregate(counts$ID[sel], by=list(scores$Event[sel]),
                   FUN=function(x) {length(unique(x))})
    out <- x$x[match(anno$EVENT, x[,1])]
    out[is.na(out)] <- 0
    out
}

exscr$sig.exonDel.clean <- { # Number of significant deletion guide pairs, neither of which is significant
                             # on its own when paired with an intergenic guide
    bad.Cas9 <- grepl("intronic",scores$Cas9_Guide_Type) & grepl("intergenic",scores$Cpf1_Guide_Type) &
        scores$rpe1_sig
    bad.Cpf1 <- grepl("intergenic",scores$Cas9_Guide_Type) & grepl("intronic",scores$Cpf1_Guide_Type) &
        scores$rpe1_sig
    sel <- which(grepl("intronic",scores$Cas9_Guide_Type) & grepl("intronic",scores$Cpf1_Guide_Type) &
                 scores$rpe1_sig &
                 !(counts$Cas9.Guide %in% counts$Cas9.Guide[bad.Cas9]) &
                 !(counts$Cpf1.Guide %in% counts$Cpf1.Guide[bad.Cpf1]))
    x <- aggregate(counts$ID[sel], by=list(scores$Event[sel]),
                   FUN=function(x) {length(unique(x))})
    out <- x$x[match(anno$EVENT, x[,1])]
    out[is.na(out)] <- 0
    out
}

exscr$frac.singleCtrl <- round(exscr$sig.singleCtrl / exscr$n.singleCtrl, 4)
exscr$frac.all        <- round(exscr$sig.exonDel / exscr$n.exonDel, 4)
exscr$frac.clean      <- round(exscr$sig.exonDel.clean / exscr$exonDel.clean, 4)


exTypes <- c("ctrl-skipped","frameDisr/essential","frameDisr/nonessential",
             "framePres/essential","framePres/nonessential")
typeCols <- c("grey80","cyan2","grey60","dodgerblue","red")


### Numbers of off-target guides and significant guides per hit ###

hits <- exscr$frac.all > hitCO & exscr$sig.exonDel.clean >= minClean
exscr$hit <- hits

write.csv(data.frame(anno, exscr[,-1]), row.names=F,
          file=file.path(outPath, "ExonDeletionScores.RPE1.csv"))


### Plot screen performance based on exon deletion pairs ###

scoreCO <- seq(0, 0.5, 0.01)
fore <- c("frameDisr/essential")
back <- c("notExpressed","ctrl-skipped")
fdis <- c(                      "frameDisr/nonessential","frameDisr")
fpre <- c("framePres/essential","framePres/nonessential","framePres")

typeCol <- c("red","dodgerblue1","grey40","grey80")


dat1 <- matrix(nrow=length(scoreCO), ncol=4)
for (i in 1:length(scoreCO)) {
    dat1[i,1] <- length(which(anno$type.RPE1 %in% fore &
                              exscr$frac.all > scoreCO[i] &
                              exscr$sig.exonDel.clean > minClean)) /
        length(which(anno$type.RPE1 %in% fore))
    dat1[i,2] <- length(which(anno$type.RPE1 %in% back &
                              exscr$frac.all > scoreCO[i] &
                              exscr$sig.exonDel.clean > minClean)) /
        length(which(anno$type.RPE1 %in% back))
    dat1[i,3] <- length(which(anno$type.RPE1 %in% fdis &
                              exscr$frac.all > scoreCO[i] &
                              exscr$sig.exonDel.clean > minClean)) /
        length(which(anno$type.RPE1 %in% fdis))
    dat1[i,4] <- length(which(anno$type.RPE1 %in% fpre &
                              exscr$frac.all > scoreCO[i] &
                              exscr$sig.exonDel.clean > minClean)) /
        length(which(anno$type.RPE1 %in% fpre))
}
dat1 <- 100 * dat1
dat1[is.na(dat1)] <- 0


pdf(file.path(outPath, "ScreenPerformance.RPE1.exonDelPairs.pdf"), wid=5, hei=5.4)
plot(1, 1, type="n", xlim=c(0, max(scoreCO)), ylim=c(0, max(dat1)),
     xlab="Cutoff for calling exon phenotype (frac. of signif. pairs)", ylab="% of exons with phenotype",
     main="Screen performance (RPE1)")
for (i in 1:ncol(dat1)) {
    lines(scoreCO, dat1[,i], col=typeCol[i], lwd=2)
}
abline(h=0)
abline(v=hitCO, lty=2)
text(hitCO + 0.005, dat1[min(which(scoreCO >= hitCO)),] + 1,
     paste0(round(dat1[min(which(scoreCO >= hitCO)),], 1), "%"), adj=c(0,0), col=typeCol)

legend("topright",
       legend=c(fore, paste(back, collapse="\n"), paste(fdis, collapse="\n"), paste(fpre, collapse="\n")),
       lwd=2, col=typeCol, bty="n")
dev.off()


### Similar, but for intronic-intergenic control guide pairs ###

dat0 <- matrix(nrow=length(scoreCO), ncol=4)
for (i in 1:length(scoreCO)) {
    dat0[i,1] <- length(which(anno$type.RPE1 %in% fore &
                              exscr$frac.singleCtrl > scoreCO[i])) /
        length(which(anno$type.RPE1 %in% fore))
    dat0[i,2] <- length(which(anno$type.RPE1 %in% back &
                              exscr$frac.singleCtrl > scoreCO[i] )) /
        length(which(anno$type.RPE1 %in% back))
    dat0[i,3] <- length(which(anno$type.RPE1 %in% fdis &
                              exscr$frac.singleCtrl > scoreCO[i])) /
        length(which(anno$type.RPE1 %in% fdis))
    dat0[i,4] <- length(which(anno$type.RPE1 %in% fpre &
                              exscr$frac.singleCtrl > scoreCO[i])) /
        length(which(anno$type.RPE1 %in% fpre))
}
dat0 <- 100 * dat0
dat0[is.na(dat0)] <- 0


pdf(file.path(outPath, "ScreenPerformance.cutoffs.RPE1.singleCutPairs.pdf"), wid=5, hei=5.4)
plot(1, 1, type="n", xlim=c(0, max(scoreCO)), ylim=c(0, max(dat0)),
     xlab="Cutoff for calling exon phenotype (frac. of signif. pairs)", ylab="% of exons with phenotype",
     main="Screen performance (RPE1)\n(intronic-intergenic pairs)")
for (i in 1:ncol(dat0)) {
    lines(scoreCO, dat0[,i], col=typeCol[i], lwd=2)
}
abline(h=0)
abline(v=hitCO, lty=2)
text(hitCO + 0.005, dat0[min(which(scoreCO >= hitCO)),] + 1,
     paste0(round(dat0[min(which(scoreCO >= hitCO)),], 1), "%"), adj=c(0,0), col=typeCol)

legend("topright",
       legend=c(fore, paste(back, collapse="\n"), paste(fdis, collapse="\n"), paste(fpre, collapse="\n")),
       lwd=2, col=typeCol, bty="n")
dev.off()



### Bar plots showing hits at chosen cutoff ###

mat1 <- mat0 <- matrix(nrow=2, ncol=2,
                       dimnames=list(c("Hits","No hits"),
                                     c(fore, paste(back, collapse="\n")))
                       )
mat1[1,1] <- length(which(anno$type.RPE1 %in% fore &
                          exscr$frac.all  > hitCO))
mat1[2,1] <- length(which(anno$type.RPE1 %in% fore &
                          exscr$frac.all <= hitCO)) 
mat1[1,2] <- length(which(anno$type.RPE1 %in% back &
                          exscr$frac.all  > hitCO)) 
mat1[2,2] <- length(which(anno$type.RPE1 %in% back &
                          exscr$frac.all <= hitCO)) 
mat0[1,1] <- length(which(anno$type.RPE1 %in% fore &
                          exscr$frac.singleCtrl  > hitCO)) 
mat0[2,1] <- length(which(anno$type.RPE1 %in% fore &
                          exscr$frac.singleCtrl <= hitCO)) 
mat0[1,2] <- length(which(anno$type.RPE1 %in% back &
                          exscr$frac.singleCtrl  > hitCO)) 
mat0[2,2] <- length(which(anno$type.RPE1 %in% back &
                          exscr$frac.singleCtrl <= hitCO)) 

datComb <- cbind(t(t(mat1[1,]) / colSums(mat1)), t(t(mat0[1,]) / colSums(mat0)))
p.1 <- signif(fisher.test(mat1)$p.value, 2)
p.0 <- signif(fisher.test(mat0)$p.value, 2)

pdf(file.path(outPath, paste0("HitCalling_", hitCO, ".pdf")),
    wid=4, hei=4.4)
xpos <- barplot(100 * datComb, beside=T, ylab="% of exons with phenotype", main="Hit calling",
                col=c("red", "dodgerblue1"), names.arg=c("Exon deletion", "Single cut"),
                legend.text=c(fore, paste(back, collapse="\n")),
                args.legend=list(bty="n")
                )
par(xpd=T)
text(colMeans(xpos), apply(datComb, MAR=2, max), paste("p =", c(p.1, p.0)), pos=3)
par(xpd=F)
abline(h=0)
dev.off()


### Plot hit properties ###

typeOK <- anno$type.RPE1 %in% c("framePres/nonessential","framePres/essential","framePres",
                                    "frameDisr/nonessential","frameDisr/essential","frameDisr") 
hitGroups <- list("Hits"    = which(exscr$frac.all  > hitCO & exscr$sig.exonDel.clean > minClean & typeOK),
                  "No hits" = which(exscr$frac.all <= hitCO & typeOK)
                  )
hitCols <- c("orange","grey")


pdf(file.path(outPath, paste0("Hits_RPE1_", hitCO, ".properties.pdf")),
    wid=5, hei=5.4)

## Essentiality score
essDat <- lapply(hitGroups, FUN=function(x) {anno$ess.RPE1.log2FC[x]})

plot(density(essDat[[2]], na.rm=T), lwd=2, col=hitCols[2],
     main="RPE1 dropout rate (TKOv3 guides)", xlab="Log2 (fold change)")
lines(density(essDat[[1]], na.rm=T), col=hitCols[1], lwd=2)
abline(v=sapply(essDat, median, na.rm=T), lty=2, col=hitCols)
legend("topleft", "Hit", text.col="orange", bty="n")
legend("left", paste("p =", signif(wilcox.test(essDat[[1]], essDat[[2]])$p.value, 2)), bty="n")

## Length
lenDat <- lapply(hitGroups, FUN=function(x) {log10(anno$LENGTH[x])})

plot(density(lenDat[[2]], na.rm=T), lwd=2, col=hitCols[2],
     main="Exon length", xlab="Log10 (bp)")
lines(density(lenDat[[1]], na.rm=T), col=hitCols[], lwd=2)
abline(v=sapply(lenDat, median, na.rm=T), lty=2, col=hitCols)
legend("topleft", "Hit", text.col="orange", bty="n")
legend("left", paste("p =", signif(wilcox.test(lenDat[[1]], lenDat[[2]])$p.value, 2)), bty="n")

dev.off()

