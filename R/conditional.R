#' Conditional analysis of CNVs
#' 
#' This function performs the conditional analysis of fCNVs
#'
#' @param x Object of class diggit
#' @param ... Additional parameters to pass to the function
#' @export
#' @docType methods
#' @rdname conditional-methods
setGeneric("conditional", function(x, ...) standardGeneric("conditional"))

#' @param pheno Character string indicating the feature for sample groups
#' @param group1 Character string indicating the treatment group
#' @param group2 Optional character string indicating the reference group
#' @param cnv Single number or vector of two numbers indicating the thresholds for CNVs
#' @param mr Either vector of character strings indicating the MR genes, or number indicating the corrected p-value threshold for selecting the MRs
#' @param mr.adjust Character string indicating the multiple-hypothesis correction to apply to the MR p-values
#' @param modul Number indicating the p-value threshold for a modulator to be considered associated with the MR activity
#' @param modul.adjust Character string indicating the multiple-hypothesis correction to apply to the aQTL results
#' @param fet.pval Number indicating the FET p-value threshold for the association between CNVs and sample groups
#' @param cores Integer indicating the number of cores to use (1 for Windows-based systems)
#' @param verbose Logical, whether progress should be reported
#' @return Object of class diggit with conditional analysis results
#' @export
#' @docType method
#' @rdname conditional-methods
#' @examples
#' data(gbm.expression, package="diggitdata")
#' data(gbm.cnv, package="diggitdata")
#' data(gbm.aracne, package="diggitdata")
#' dobj <- diggitClass(expset=gbmExprs, cnv=gbmCNV, regulon=gbmTFregulon)
#' dobj <- fCNV(dobj)
#' dobj <- aqtl(dobj, mr=c("CEBPD", "STAT3"), fcnv.adjust="fdr", verbose=FALSE)
#' dobj <- conditional(dobj, pheno="subtype", group1="MES", group2="PN", mr="STAT3", verbose=FALSE)
#' dobj
setMethod("conditional", c(x="diggit"), function(x, pheno="cond", group1, group2=NULL, cnv=.2, mr=.01, mr.adjust=c("none", "fdr", "bonferroni"), modul=.01, modul.adjust=c("none", "fdr", "bonferroni"), fet.pval=0.05, cores=1, verbose=TRUE) {
    mr.adjust <- match.arg(mr.adjust)
    modul.adjust <- match.arg(modul.adjust)
    if (length(cnv)==1) cnv <- c(-abs(cnv), abs(cnv))
    if (is.numeric(mr)) {
        if (abs(mr)<1) {
            tmp <- p.adjust(pnorm(abs(diggitMR(x)), lower.tail=FALSE)*2, mr.adjust)
            mr <- names(tmp)[tmp<mr]
        }
    }
    mr <- mr[mr %in% rownames(diggitAqtl(x))]
    if (length(mr)==0) stop("No MR satisfying conditions for the analysis", call.=FALSE)
    pb <- NULL
    if (cores==1 | length(mr)<cores) {
        if (verbose) pb <- txtProgressBar(max=length(mr), style=3)
        res <- lapply(1:length(mr), function(i, mr, x, pheno, group1, group2, cnv, modul, modul.adjust, fet.pval, pb) {
            if (!is.null(pb)) setTxtProgressBar(pb, i)
            mr <- mr[i]
            conditionalInternal(x, mr=mr, pheno=pheno, group1=group1, group2=group2, cnv=cnv, modul=modul, modul.adjust=modul.adjust, fet.pval=fet.pval)
        }, mr=mr, x=x, pheno=pheno, group1=group1, group2=group2, cnv=cnv, modul=modul, modul.adjust=modul.adjust, fet.pval=fet.pval, pb=pb)
    }
    else {
        res <- mclapply(1:length(mr), function(i, mr, x, pheno, group1, group2, cnv, modul, modul.adjust, fet.pval) {
            mr <- mr[i]
            conditionalInternal(x, mr=mr, pheno=pheno, group1=group1, group2=group2, cnv=cnv, modul=modul, modul.adjust=modul.adjust, fet.pval=fet.pval)
        }, mr=mr, x=x, pheno=pheno, group1=group1, group2=group2, cnv=cnv, modul=modul, modul.adjust=modul.adjust, fet.pval=fet.pval, mc.cores=cores)
    }
    names(res) <- mr
    x@conditional <- res
    return(x)
})

conditionalInternal <- function(x, mr, pheno, group1, group2, cnv, modul, modul.adjust, fet.pval) {
    # Clustering
    aqtemp <- diggitAqtl(x)
    aqtemp <- filterRowMatrix(aqtemp, rownames(aqtemp) %in% mr)
    aqtemp <- filterColMatrix(aqtemp, !is.na(colSums(aqtemp)))
    if (is.numeric(modul)) {
        if (abs(modul)<1) {
            tmp <- p.adjust(aqtemp, modul.adjust)
            modul <- colnames(aqtemp)[tmp<modul]
        }
    }
    modul <- modul[modul %in% colnames(aqtemp)]
    if (length(modul)==0) stop("No modulators satisfying conditions for the analysis", call.=FALSE)
    if (is.null(group2)) {
        pos <- pData(exprs(x))[[pheno]] %in% group1
        if (length(pos)==0) stop(paste(pheno, " was not found in the ExpressionSet Object", sep=""), call.=FALSE)
        if (length(which(pos))==0) stop(paste(group1, " was not found in ", pheno, sep=""), call.=FALSE)
        samp <- colnames(exprs(exprs(x)))
        samp1 <- samp[pos]
    }
    else  {
        pos1 <- pData(exprs(x))[[pheno]] %in% group1
        pos2 <- pData(exprs(x))[[pheno]] %in% group2
        if (length(pos1)==0) stop(paste(pheno, " was not found in the ExpressionSet Object", sep=""), call.=FALSE)
        if (length(which(pos1))==0) stop(paste(group1, " was not found in ", pheno, sep=""), call.=FALSE)
        if (length(which(pos2))==0) stop(paste(group2, " was not found in ", pheno, sep=""), call.=FALSE)
        samp1 <- colnames(exprs(exprs(x)))[pos1]
        samp2 <- colnames(exprs(exprs(x)))[pos2]
        samp <- c(samp1, samp2)
    }
    cnvmat <- diggitCNV(x)
    cnvmat <- filterRowMatrix(cnvmat, match(modul, rownames(cnvmat)))
    samp <- intersect(colnames(cnvmat), samp)
    cnvmat <- filterColMatrix(cnvmat, match(samp, colnames(cnvmat)))
    samp1 <- samp1[samp1 %in% samp]
    if (length(samp1)==0 | length(samp1)==length(samp)) stop("Not enough samples per condition", call.=FALSE)

# Test for significantly associated CNVs
# deletions
    tmp <- cnvmat < cnv[1]
    tmp <- filterRowMatrix(tmp, rowSums(tmp)>0)
    fetr <- apply(tmp, 1, function(x, test) {
        fisher.test(test, x, alternative="greater")$p.value
    }, test=samp %in% samp1)
    mod.del <- fetr[fetr<fet.pval]
# amplifications
    tmp <- cnvmat > cnv[2]
    tmp <- filterRowMatrix(tmp, rowSums(tmp)>0)
    fetr <- apply(tmp, 1, function(x, test) {
        fisher.test(test, x, alternative="greater")$p.value
    }, test=samp %in% samp1)
    mod.amp <- fetr[fetr<fet.pval]

# Clustering and conditional analysis
#deletions
res.del <- NULL
if (length(mod.del)>2) {
    tmp <- t(filterRowMatrix(cnvmat, rownames(cnvmat) %in% names(mod.del)))
    tmp <- tmp < cnv[1]
    ji <- jaccard(filterColMatrix(tmp, colSums(tmp)>0))
    tmp <- ji
    diag(tmp) <- 0
    pos <- rowSums(tmp)>0
    ji <- filterColMatrix(filterRowMatrix(ji, pos), pos)
    hc <- hclust(as.dist(1-ji), "average")
    groups <- cutree(hc, h=1-distMode(ji[ji>0 & ji<1]))
    res.del <- tapply(names(groups), groups, function(gene, samp2, cnv, cnvmat) {
        tmp1 <- filterRowMatrix(cnvmat, match(gene, rownames(cnvmat)))
        tmp1 <- tmp1 < cnv
        res <- apply(tmp1, 1, function(pos, tmp1, samp2) {
            apply(tmp1, 1, function(tmp1, pos, samp2) {
                tmp <- NA
                try(tmp <- fisher.test(samp2[!pos], tmp1[!pos], alternative="greater")$p.value, silent=TRUE)
                tmp
            }, samp2=samp2, pos=pos)
        }, tmp1=tmp1, samp2=samp2)
        if (!is.null(nrow(res))) {
            pos <- order(rowMeans(res, na.rm=TRUE))
            res <- res[pos, ][, pos]
        }
        return(res)
    }, samp2=samp %in% samp1, cnv=cnv[1], cnvmat=cnvmat)
}
else if (length(mod.del)>0) {
    tmp1 <- filterRowMatrix(cnvmat, match(names(mod.del), rownames(cnvmat)))
    tmp1 <- tmp1 < cnv[1]        
    res <- apply(tmp1, 1, function(pos, tmp1, samp2) {
        apply(tmp1, 1, function(tmp1, pos, samp2) {
            tmp <- NA
            try(tmp <- fisher.test(samp2[!pos], tmp1[!pos], alternative="greater")$p.value, silent=TRUE)
            tmp
        }, samp2=samp2, pos=pos)
    }, tmp1=tmp1, samp2=samp %in% samp1)
    if (!is.null(nrow(res))) {
        pos <- order(rowMeans(res, na.rm=TRUE))
        res <- res[pos, ][, pos]
    }
    res.del <- list("1"=res)
}
# Amplifications
    res.amp <- NULL
    if (length(mod.amp)>2) {
        tmp <- t(filterRowMatrix(cnvmat, rownames(cnvmat) %in% names(mod.amp)))
        tmp <- tmp > cnv[2]
        ji <- jaccard(filterColMatrix(tmp, colSums(tmp)>0))
        tmp <- ji
        diag(tmp) <- 0
        pos <- rowSums(tmp)>0
        ji <- filterColMatrix(filterRowMatrix(ji, pos), pos)
        hc <- hclust(as.dist(1-ji), "average")
        groups <- cutree(hc, h=1-distMode(ji[ji>0 & ji<1]))
        res.amp <- tapply(names(groups), groups, function(gene, samp2, cnv, cnvmat) {
            tmp1 <- filterRowMatrix(cnvmat, match(gene, rownames(cnvmat)))
            tmp1 <- tmp1 > cnv
            res <- apply(tmp1, 1, function(pos, tmp1, samp2) {
                apply(tmp1, 1, function(tmp1, pos, samp2) {
                    tmp <- NA
                    try(tmp <- fisher.test(samp2[!pos], tmp1[!pos], alternative="greater")$p.value, silent=TRUE)
                    tmp
                }, samp2=samp2, pos=pos)
            }, tmp1=tmp1, samp2=samp2)
            if (!is.null(nrow(res))) {
                pos <- order(rowMeans(res, na.rm=TRUE))
                res <- res[pos, ][, pos]
            }
            return(res)
        }, samp2=samp %in% samp1, cnv=cnv[2], cnvmat=cnvmat)
    }
    else if (length(mod.amp)>0) {
        tmp1 <- filterRowMatrix(cnvmat, match(names(mod.amp), rownames(cnvmat)))
        tmp1 <- tmp1 > cnv[2]
        res <- apply(tmp1, 1, function(pos, tmp1, samp2) {
            apply(tmp1, 1, function(tmp1, pos, samp2) {
                tmp <- NA
                try(tmp <- fisher.test(samp2[!pos], tmp1[!pos], alternative="greater")$p.value, silent=TRUE)
                tmp
            }, samp2=samp2, pos=pos)
        }, tmp1=tmp1, samp2=samp %in% samp1)
        if (!is.null(nrow(res))) {
            pos <- order(rowMeans(res, na.rm=TRUE))
            res <- res[pos, ][, pos]
        }
        res.amp <- list("1"=res)
    }
    tmp <- lapply(c(res.del, res.amp), function(x) {
        if (length(x)>1) return(apply(x, 1, max, na.rm=TRUE))
        return(x)
    })
    integ <- unlist(tmp, use.names=FALSE)
    names(integ) <- unlist(lapply(tmp, names), use.names=FALSE)
    mod.all <- c(mod.del, mod.amp)
    integ <- cbind(integ, mod.all[match(names(integ), names(mod.all))])
    integ <- apply(integ, 1, max, na.rm=TRUE)
    list(deleted=res.del, amplified=res.amp, original=sort(mod.all), integrated=sort(integ))    
}

