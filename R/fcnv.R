#' Inference of functional CNVs
#' 
#' This function infers functional CNVs by computing their association with gene expression
#' 
#' @param x Object of class diggit, expressionSet object or numeric matrix of expression data, with features in rows and samples in columns
#' @param ... Additional arguments
#' @return Objet of class diggit with updated fCNV slot
#' @export 
#' @docType methods
#' @rdname fcnv-methods
setGeneric("fCNV", function(x, ...) standardGeneric("fCNV"))

#' @param expset Optional numeric matrix of expression data
#' @param cnv Optional numeric matrix of CNVs
#' @param method Character string indicating the method for computing the association between CNVs and expression
#' @param cores Integer indicating the number of cores to use (1 for Windows-based systems)
#' @param verbose Logical, whether to report analysis progress
#' @examples
#' data(gbm.expression, package="diggitdata")
#' data(gbm.cnv, package="diggitdata")
#' genes <- intersect(rownames(gbmExprs), rownames(gbmCNV))[1:100]
#' gbmCNV <- gbmCNV[match(genes, rownames(gbmCNV)), ]
#' dgo <- diggitClass(expset=gbmExprs, cnv=gbmCNV)
#'
#' dgo <- fCNV(dgo)
#' dgo
#' diggitFcnv(dgo)[1:5]
#' @rdname fcnv-methods
#' @export
setMethod("fCNV", c(x="diggit"), function(x, expset=NULL, cnv=NULL, method=c("spearman", "mi", "pearson", "kendall"), cores=1, verbose=TRUE) {
    method <- match.arg(method)
    if (is.null(expset)) expset <- exprs(x)
    if (is.null(cnv)) cnv <- diggitCNV(x)
    tmp <- fCNV(x=exprs(expset), cnv=cnv, method=method, cores=cores, verbose=verbose)
    x@fcnv <- diggitFcnv(tmp)
    x@expset <- expset
    x@cnv <- cnv
    return(x)
})

#' @examples
#'
#' dgo <- fCNV(gbmExprs, gbmCNV)
#' print(dgo)
#' diggitFcnv(dgo)[1:5]
#' @rdname fcnv-methods
#' @export
setMethod("fCNV", c(x="ExpressionSet"), function(x, cnv, method=c("spearman", "mi", "pearson", "kendall"), cores=1, verbose=TRUE) {
    method <- match.arg(method)
    tmp <- fCNV(x=exprs(x), cnv=cnv, method=method, cores=cores, verbose=verbose)
    tmp@expset <- x
    return(tmp)
})

#' @examples
#'
#' dgo <- fCNV(exprs(gbmExprs), gbmCNV)
#' dgo
#' diggitFcnv(dgo)[1:5]
#' @rdname fcnv-methods
#' @export
setMethod("fCNV", c(x="matrix"), function(x, cnv, method=c("spearman", "mi", "pearson", "kendall"), cores=1, verbose=TRUE) {
    method <- match.arg(method)
# Compatible structure for expset and cnv
    genes <- intersect(rownames(x), rownames(cnv))
    samples <- intersect(colnames(x), colnames(cnv))
    x <- x[match(genes, rownames(x)), ][, match(samples, colnames(x))]
    cnv <- cnv[match(genes, rownames(cnv)), ][, match(samples, colnames(cnv))]
# Compute MI
    if (method %in% c("mi")) {
        mi <- mutualInfo(t(x), t(cnv), per=100, pairwise=TRUE, cores=cores, verbose=verbose)
    }
    else mi <- correlation(t(x), t(cnv), method=method, pairwise=TRUE)
# Generate the output diggit object
    diggitClass(expset=ExpressionSet(assayData=x), cnv=cnv, fcnv=mi$p.value)
})

#' @examples
#'
#' dgo <- fCNV(as.data.frame(exprs(gbmExprs)), gbmCNV)
#' dgo
#' diggitFcnv(dgo)[1:5]
#' @rdname fcnv-methods
#' @export
setMethod("fCNV", c(x="data.frame"), function(x, cnv, method=c("spearman", "mi", "pearson", "kendall"), cores=1, verbose=TRUE) {
    method <- match.arg(method)
    tmp <- fCNV(x=as.matrix(x), cnv=cnv, method=method, cores=cores, verbose=verbose)
    return(tmp)
})


