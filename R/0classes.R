# Clases
setOldClass("regulon")

#'The diggit class
#'
#'This class stores parameters and results of the diggit algorithm
#'  
#'@section Slots:
#'    \describe{
#'      \item{\code{expset}:}{ExpressionSet object containing the gene expression data}
#'      \item{\code{cnv}:}{Matrrix containing the CNV data}
#'      \item{\code{regulon}:}{Regulon object containing the transcriptional interactome}
#'      \item{\code{mindy}:}{Regulon object containing the post-translational interactome}
#'      \item{\code{fcnv}:}{Numeric vector containing the p-values for functional CNVs}
#'      \item{\code{mr}:}{Numeric vector of normalized enrichment scores for the MARINa analysis}
#'      \item{\code{viper}:}{Numeric matrix of normalized enrichment scores for the VIPER analysis}
#'      \item{\code{aqtl}:}{Numeric matrix of association p-values for the aQTL analysis}
#'      \item{\code{conditional}:}{List containing the conditional analysis results}
#'}
#'@rdname diggit-class
#'@details see \code{\link[=print.diggit]{diggit-methods} for related methods}
setClass("diggit", slots=c(expset="ExpressionSet", cnv="matrix", regulon="regulon", mindy="regulon", fcnv="vector", mr="vector", viper="matrix", aqtl="matrix", conditional="list"))

#' The diggit class constructor
#' 
#' This function generates diggit class objects
#' 
#' @param expset ExpressionSet object or numeric matrix of expression data, with features in rows and samples in columns
#' @param cnv Numeric matrix of CNV data
#' @param regulon Regulon class object containing the transcriptional interactome
#' @param mindy Regulon class object containing the post-translational interactome
#' @param fcnv Vector of F-CNV p-values
#' @param mr Vector of master regulator Z-score (NES)
#' @param viper Numeric matrix of VIPER results
#' @param aqtl Numeric matrix of aQTL p-values
#' @param conditional List containing the conditional analysis results
#' @return Object of class diggit
#' @examples
#' data(gbm.expression, package="diggitdata")
#' data(gbm.aracne, package="diggitdata")
#' dobj <- diggitClass(expset=gbmExprs, regulon=gbmTFregulon)
#' print(dobj)
#' @rdname diggit-class
#' @aliases diggitClass
#' @export
diggitClass <- function(expset=NULL, cnv=NULL, regulon=NULL, mindy=NULL, fcnv=NULL, mr=NULL, viper=NULL, aqtl=NULL, conditional=NULL) {
    if (is.null(expset)) expset <- ExpressionSet()
    if (is.null(cnv)) cnv <- matrix(, 0, 0)
    if (is.null(regulon)) {
        regulon <- list()
        class(regulon) <- "regulon"
    }
    if (is.null(mindy)) {
        mindy <- list()
        class(mindy) <- "regulon"
    }
    if (is.null(fcnv)) fcnv <- numeric(0)
    if (is.null(mr)) mr <- numeric(0)
    if (is.null(viper)) viper <- matrix(, 0, 0)
    if (is.null(aqtl)) aqtl <- matrix(, 0, 0)
    if (is.null(conditional)) conditional=list()
    new("diggit", expset=expset, cnv=cnv, regulon=regulon, mindy=mindy, fcnv=fcnv, mr=mr, viper=viper, aqtl=aqtl, conditional=conditional)
}

#' Basic methods for class diggit
#' 
#' This document lists a series of basic methods for the class diggit
#'
#' @param x Object of class diggit
#' @param object Object of class diggit
#' @param pval P-value threshold for the conditional analysis
#' @return print returns summary information about the diggit object
#' @rdname diggit-methods
#' @aliases print.diggit
#' @docType methods
#' @examples
#' data(gbm.expression, package="diggitdata")
#' data(gbm.cnv, package="diggitdata")
#' data(gbm.aracne, package="diggitdata")
#' dobj <- diggitClass(expset=gbmExprs, cnv=gbmCNV, regulon=gbmTFregulon)
#' print(dobj)
#' show(dobj)
#' exprs(dobj)
#' diggitCNV(dobj)[1:3, 1:3]
#' diggitRegulon(dobj)
#' diggitMindy(dobj)
#' diggitFcnv(dobj)
#' diggitMR(dobj)
#' diggitViper(dobj)
#' diggitAqtl(dobj)
#' diggitConditional(dobj)
#' head(dobj)
#' @export
setMethod("print", "diggit", function(x, pval=.05) {
    cat("\nAn object of class diggit\n")
    cat("Slot expset:\n")
    if (nrow(exprs(exprs(x)))>0) cat("Expression data of", nrow(exprs(exprs(x))), "features by", ncol(exprs(exprs(x))), "samples\n")
    else cat("Empty\n")
    cat("\nSlot cnv:\n")
    if (nrow(diggitCNV(x))>0) cat("CNV data of", nrow(diggitCNV(x)), "features by", ncol(diggitCNV(x)), "samples\n")
    else cat("Empty\n")
    cat("\nSlot regulon:\n")
    if (length(diggitRegulon(x))>0) print(diggitRegulon(x))
    else cat("Empty\n")
    cat("\nSlot mindy:\n")
    if (length(diggitMindy(x))>0) print(diggitMindy(x))
    else cat("Empty\n")
    cat("\nSlot fcnv:\n")
    if (length(diggitFcnv(x))>0) cat("F-CNV statistical significance for", length(diggitFcnv(x)), "features\n")
    else cat("Empty\n")
    cat("\nSlot mr:\n")
    if (length(diggitMR(x))>0) cat("Master regulator NES for", length(diggitMR(x)), "features\n")
    else cat("Empty\n")
    cat("\nSlot viper:\n")
    if (nrow(diggitViper(x))>0) cat("VIPER protein activity matrix for", nrow(diggitViper(x)), "features by", ncol(diggitViper(x)), "samples\n")
    else cat("Empty\n")
    cat("\nSlot aqtl:\n")
    if (nrow(diggitAqtl(x))>0) cat("aQTL matrix for", nrow(diggitAqtl(x)), "regulators by", ncol(diggitAqtl(x)), "genetic alterations\n")
    else cat("Empty\n")
    cat("\nSlot conditional:\n")
    if (length(diggitConditional(x))>0) {
        tmp <- diggitConditional(x)
        cat("Conditional analysis results for ", length(tmp), " MRs and ", length(unlist(lapply(tmp, function(x) x$integ), use.names=FALSE)), " modulators:\n\n", sep="")
        for (i in 1:length(tmp)) {
            cat(names(tmp)[i], "\n")
            cat(length(which(tmp[[i]]$integrated<pval)), " modulators at p<", pval, ":\n", sep="")
            cat(paste(names(tmp[[i]]$integrated)[tmp[[i]]$integrated<pval], collapse=", "), "\n\n")
        }
    }
    else cat("Empty\n")
})

#' @rdname diggit-methods
#' @aliases show.diggit
#' @return show returns summary information about the object of class diggit
#' @export
#' 
setMethod("show", "diggit", function(object) {
    print(object)
})


#' @rdname diggit-methods
#' @aliases exprs.diggit
#' @return exprs returns the ExpressionSet object containing the expression profile data
#' @export
setMethod("exprs", "diggit", function(object) {
    object@expset
})

setGeneric("diggitCNV", function(x, ...) standardGeneric("diggitCNV"))

#' @rdname diggit-methods
#' @aliases diggitCNV
#' @return diggitCNV returns a matrix containing the CNV data
#' @export
setMethod("diggitCNV", "diggit", function(x) {
    x@cnv
})

setGeneric("diggitRegulon", function(x, ...) standardGeneric("diggitRegulon"))

#' @rdname diggit-methods
#' @aliases diggitRegulon
#' @return diggitRegulon returns a regulon object containing the transcriptional interactome
#' @export
setMethod("diggitRegulon", "diggit", function(x) {
    x@regulon
})

setGeneric("diggitMindy", function(x, ...) standardGeneric("diggitMindy"))

#' @rdname diggit-methods
#' @aliases diggitMindy
#' @return diggitMindy returns a regulon object containing the post-translational interactome
#' @export
setMethod("diggitMindy", "diggit", function(x) {
    x@mindy
})

setGeneric("diggitFcnv", function(x, ...) standardGeneric("diggitFcnv"))

#' @rdname diggit-methods
#' @aliases diggitFcnv
#' @return diggitFcnv returns a vector of p-values for the F-CNVs
#' @export
setMethod("diggitFcnv", "diggit", function(x) {
    x@fcnv
})

setGeneric("diggitMR", function(x, ...) standardGeneric("diggitMR"))

#' @rdname diggit-methods
#' @aliases diggitMR
#' @return diggitMR returns a vector of master regulators NES
#' @export
setMethod("diggitMR", "diggit", function(x) {
    x@mr
})

setGeneric("diggitViper", function(x, ...) standardGeneric("diggitViper"))

#' @rdname diggit-methods
#' @aliases diggitViper
#' @return diggitViper returns a matrix of VIPER results
#' @export
setMethod("diggitViper", "diggit", function(x) {
    x@viper
})

setGeneric("diggitAqtl", function(x, ...) standardGeneric("diggitAqtl"))

#' @rdname diggit-methods
#' @aliases diggitAqtl
#' @return diggitAqtl returns a matrix of aQTLs (p-value)
#' @export
setMethod("diggitAqtl", "diggit", function(x) {
    x@aqtl
})

setGeneric("diggitConditional", function(x, ...) standardGeneric("diggitConditional"))

#' @rdname diggit-methods
#' @aliases diggitConditional
#' @return diggitConditional returns a list containing the conditional analysis results
#' @export
setMethod("diggitConditional", "diggit", function(x) {
    x@conditional
})

#' @rdname diggit-methods
#' @aliases summary.diggit
#' @return summary returns the integrated results from the conditional analysis
#' @export
setMethod("summary", "diggit", function(object) {
    x <- diggitConditional(object)
    if (length(x)==0) stop("No conditional analysis results found", call.=FALSE)
    tmp <- lapply(x, function(x) x$integrated)
    return(tmp)
})

#' @param rows Integer indicating the maximum number of rows to show
#' @param cols Integer indicating the maximum number of columns to show
#' @rdname diggit-methods
#' @aliases head.diggit
#' @return head returns a list containing a reduced view for an object of class diggit
#' @export
setMethod("head", "diggit", function(x, rows=4, cols=4) {
    res <- list()
    if (nrow(exprs(exprs(x)))>0) {
        tmp <- exprs(exprs(x))
        res$expset <- filterRowMatrix(filterColMatrix(tmp, 1:min(ncol(tmp), cols)), 1:min(nrow(tmp), rows))
    }
    if (nrow(diggitCNV(x))>0) {
        tmp <- diggitCNV(x)
        res$cnv <- filterRowMatrix(filterColMatrix(tmp, 1:min(ncol(tmp), cols)), 1:min(nrow(tmp), rows))
    }
    if (length(diggitFcnv(x))>0) {
        tmp <- sort(diggitFcnv(x))
        res$fcnv <- tmp[1:min(length(tmp), rows)]
    }
    if (length(diggitMR(x))>0) {
        tmp <- sort(diggitMR(x), decreasing=TRUE)
        res$mr <- tmp[1:min(length(tmp), rows)]
    }
    if (nrow(diggitViper(x))>0) {
        tmp <- diggitViper(x)
        res$viper <- filterRowMatrix(filterColMatrix(tmp, 1:min(ncol(tmp), cols)), 1:min(nrow(tmp), rows))
    }    
    if (nrow(diggitAqtl(x))>0) {
        tmp <- diggitAqtl(x)
        res$aqtl <- filterRowMatrix(filterColMatrix(tmp, 1:min(ncol(tmp), cols)), 1:min(nrow(tmp), rows))
    }    
    if (length(diggitConditional(x))>0) {
        tmp <- diggitConditional(x)
        tmp <- lapply(x, function(x) x$integrated)
        tmp <- tmp[min(length(tmp), rows)]
        tmp <- lapply(tmp, function(x, cols) {
            x[1:min(length(x), cols)]
        }, cols=cols)
        res$conditional <- tmp
    }
    return(res)
})

setGeneric("mindyFiltering", function(x, ...) standardGeneric("mindyFiltering"))

#' @rdname diggit-methods
#' @aliases mindyFiltering
#' @param mr Either a numerical value between 0 and 1 indicating the p-value threshold for the Master Regulator (MR) selection, or a vector of character strings listing the MRs
#' @param mr.adjust Character string indicating the multiple hypothesis test correction for the MRs
#' @examples
#' data(gbm.expression, package="diggitdata")
#' data(gbm.cnv, package="diggitdata")
#' data(gbm.mindy, package="diggitdata")
#' dobj <- diggitClass(expset=gbmExprs, cnv=gbmCNV, mindy=gbmMindy)
#' dobj <- fCNV(dobj)
#' dobj
#' dobj <- mindyFiltering(dobj, mr=c("STAT3", "CEBPD"))
#' dobj
#' @return mindyFiltering returns a diggit class object with CNV and aQTL slots filtered to contain only MINDy post-translational modulators of the MRs
#' @export
setMethod("mindyFiltering", c(x="diggit"), function(x, mr=.01, mr.adjust=c("none", "fdr", "bonferroni")) {
    mr.adjust <- match.arg(mr.adjust)
    tmp <- diggitMindy(x)
    mindy1 <- cbind(rep(names(tmp), sapply(tmp, function(x) length(x$tfmode))), unlist(lapply(tmp, function(x) names(x$tfmode)), use.names=FALSE))
    mindy <- names(diggitFcnv(x))
    if (is.numeric(mr)) {
        if (length(mr)==1) {
            if (mr >= 0 & mr <= 1 & length(diggitMR(x))>0) {
                tmp1 <- p.adjust(diggitMR(x), method=mr.adjust)
                mr <- names(tmp1)[tmp1 < mr]
            }
        }
    }
    if (length(which(mindy1[, 2] %in% mr))>0) {
        mindy <- unique(mindy1[mindy1[, 2] %in% mr, 1])
    }    
    if (length(diggitFcnv(x))>0) {
        tmp <- diggitFcnv(x)
        tmp <- tmp[names(tmp) %in% mindy]
        x@fcnv <- tmp
    }
    if (ncol(diggitAqtl(x))>0) {
        tmp <- diggitAqtl(x)
        tmp <- filterColMatrix(tmp, colnames(tmp) %in% mindy)
        tmp1 <- t(sapply(1:nrow(tmp), function(i, tmp, mindy) {
            mod <- mindy[mindy[, 2]==(rownames(tmp)[i]), 1]
            tmp1 <- tmp[i, ]
            tmp1[!(colnames(tmp) %in% mod)] <- NA
            tmp1
        }, tmp=tmp, mindy=mindy1))
        rownames(tmp1) <- rownames(tmp)
        x@aqtl <- tmp1
    }
    return(x)
})

