#' Inference of aQTL
#' 
#' This function infers aQTLs from F-CNVs and VIPER activity
#' 
#' @param x Object of class diggit
#' @param ... Additional parameters to pass to the function
#' @return Updated diggit object with viper and aqtl slots
#' @export
#' @docType methods
#' @rdname aqtl-methods
setGeneric("aqtl", function(x, ...) standardGeneric("aqtl"))

#' @param mr Either a numerical value between 0 and 1 indicating the p-value threshold for the Master Regulator (MR) selection, or a vector of character strings listing the MRs
#' @param mr.adjust Character string indicating the multiple hypothesis test correction for the MRs
#' @param fcnv Either a numerical value between 0 and 1 indicating the p-value threshold for the F-CNV, or a vector of character strings listing the F-CNVs
#' @param fcnv.adjust Character string indicating the multiple hypothesis test correction for the F-CNVs
#' @param method Character string indicating the method for computing the association between F-CNV and regulator activity (aQTL analysis)
#' @param mindy Logical, whether only post-translational modulators of each evaluated TF should be considered as putative genetic driver
#' @param cores Integer indicating the number of cores to use (1 for Windows-based systems)
#' @param verbose Logical, whether progress should be reported
#' @examples
#' data(gbm.expression, package="diggitdata")
#' data(gbm.cnv, package="diggitdata")
#' data(gbm.aracne, package="diggitdata")
#' dobj <- diggitClass(expset=gbmExprs, cnv=gbmCNV, regulon=gbmTFregulon)
#' dobj <- fCNV(dobj)
#' dobj <- aqtl(dobj, mr=c("CEBPD", "STAT3"), fcnv.adjust="fdr")
#' dobj
#' diggitAqtl(dobj)[, 1:4]
#' @rdname aqtl-methods
#' @export
setMethod("aqtl", c(x="diggit"), function(x, mr=.01, mr.adjust=c("none", "fdr", "bonferroni"), fcnv=.01, fcnv.adjust=c("none", "fdr", "bonferroni"), method=c("spearman", "mi", "pearson", "kendall"), mindy=FALSE, cores=1, verbose=TRUE) {
    mr.adjust <- match.arg(mr.adjust)
    fcnv.adjust <- match.arg(fcnv.adjust)
    if (mindy) x <- mindyFiltering(x, mr=mr, mr.adjust=mr.adjust)
    method <- match.arg(method)
    if (is.numeric(mr)) {
        if (abs(mr)<1) {
            tmp <- p.adjust(pnorm(abs(diggitMR(x)), lower.tail=FALSE)*2, mr.adjust)
            mr <- names(tmp)[tmp<mr]
        }
    }
    mr <- mr[mr %in% names(diggitRegulon(x)) | mr %in% rownames(diggitViper(x))]
    if (is.numeric(fcnv)) {
        if (abs(fcnv)<1) {
            tmp <- p.adjust(diggitFcnv(x), fcnv.adjust)
            fcnv <- names(tmp)[tmp<fcnv]
        }
    }
    fcnv <- fcnv[fcnv %in% names(diggitFcnv(x))]    
    if (length(mr)==0) stop("No Master Regulator fulfills the conditions", call.=FALSE)
    if (length(fcnv)==0) stop("No F-CNV fulfills the conditions", call.=FALSE)
# running VIPER
    if (all(mr %in% rownames(diggitViper(x)))) {
        act <- diggitViper(x)
        act <- filterRowMatrix(act, match(mr, rownames(act)))
    }
    else {
        regul <- diggitRegulon(x)
        regul <- regul[names(regul) %in% mr]
        act <- viper(exprs(exprs(x)), regul, method="scale", eset.filter=FALSE, cores=cores, verbose=verbose)
        x@viper <- act
    }
# aQTLs
    cnv <- diggitCNV(x)
    cnv <- cnv[rownames(cnv) %in% fcnv, ]
    samples <- intersect(colnames(act), colnames(cnv))
    cnv <- filterColMatrix(cnv, match(samples, colnames(cnv)))
    act <- filterColMatrix(act, match(samples, colnames(act)))
    if (method %in% c("mi")) {
        mi <- mutualInfo(t(act), t(cnv), per=100, pairwise=FALSE, cores=cores, verbose=verbose)
    }
    else mi <- correlation(t(act), t(cnv), method=method, pairwise=FALSE)
# Generate the output diggit object
    x@aqtl <- mi$p.value
    if (mindy) x <- mindyFiltering(x, mr=mr, mr.adjust=mr.adjust)
    return(x)
})
