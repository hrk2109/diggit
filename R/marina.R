#' Inference of Master Regulators
#' 
#' This function infers the master regulators for the transition between two phenotypes
#' 
#' @param x Object of class diggit, expressionSet object or numerical matrix containing the test samples
#' @param ... Additional arguments
#' @return Updated diggit object with Master Regulator results
#' @export
#' @docType methods
#' @rdname marina-methods
setGeneric("marina", function(x, ...) standardGeneric("marina"))

#' @param y Numerical matrix containing the control samples
#' @param mu Number indicating the control mean when \code{y} is ommited
#' @param regulon Transcriptional interactome
#' @param per Interger indicating the number of permutations to compute the marina null model
#' @param cores Integer indicating the number of cores to use (1 for Windows-based systems)
#' @param verbose Logical, whether progress should be reported
#' @examples
#' cores <- 3*(Sys.info()[1] != "Windows")+1
#' data(gbm.expression, package="diggitdata")
#' data(gbm.aracne, package="diggitdata")
#'
#' eset <- exprs(gbmExprs)
#' samples <- pData(gbmExprs)[["subtype"]]
#' x <- eset[, samples=="MES"]
#' y <- eset[, samples=="PN"]
#' dgo <- marina(x, y, regulon=gbmTFregulon, per=100, cores=cores)
#' dgo
#' diggitMR(dgo)[1:5]
#' @rdname marina-methods
#' @export
setMethod("marina", c(x="matrix"), function(x, y=NULL, mu=0, regulon, per=1000, cores=1, verbose=TRUE) {
    dnull <- NULL
    tt <- rowTtest(x=x, y=y, mu=mu)
    tt <- (qnorm(tt$p.value/2, lower.tail=FALSE)*sign(tt$statistic))[, 1]
    if (!is.null(y)) {
        if (ncol(x)>6 & ncol(y)>6) {
            dnull <- ttestNull(x, y, repos=FALSE, seed=0, per=per, cores=cores, verbose=verbose)
        }
        else {
            dnull <- ttestNull(x, y, repos=TRUE, seed=0, per=per, cores=cores, verbose=verbose)
        }
    }
    mr <- msviper(ges=tt, regulon=regulon, nullmodel=dnull, cores=cores, verbose=verbose)$es$nes
    if (is.null(y)) {
        d1 <- ExpressionSet(assayData=x)
    }
    else{
        d1 <- cbind(x, y)
        d2 <- data.frame(cond=rep(c("treat", "ctrl"), c(ncol(x), ncol(y))))
        rownames(d2) <- colnames(d1)
        d1 <- ExpressionSet(assayData=d1, phenoData=AnnotatedDataFrame(d2))
    }
    tmp <- diggitClass(expset=d1, regulon=regulon, mr=mr)
    return(tmp)
})

#' @param pheno Character string indicating the phenotype data to use
#' @param group1 Vector of character strings indicating the category from phenotype \code{pheno} to use as test group
#' @param group2 Vector of character strings indicating the category from phenotype \code{pheno} to use as control group
#' @examples
#'
#' dgo <- marina(gbmExprs, pheno="subtype", group1="MES", group2="PN", regulon=gbmTFregulon, per=100, cores=cores)
#' dgo
#' diggitMR(dgo)[1:5]
#' @rdname marina-methods
#' @export
setMethod("marina", c(x="ExpressionSet"), function(x, pheno="cond", group1, group2=NULL, mu=0, regulon, per=1000, cores=1, verbose=TRUE) {
    if (is.null(group2)) {
        pos <- pData(x)[[pheno]] %in% group1
        if (length(pos)==0) stop(paste(pheno, " was not found in the ExpressionSet Object", sep=""), call.=FALSE)
        if (length(which(pos))==0) stop(paste(group1, " was not found in ", pheno, sep=""), call.=FALSE)
        tmp <- marina(exprs(x)[, pos], mu=mu, regulon=regulon, per=per, cores=cores, verbose=verbose)
        tmp@expset <- x
        return(tmp)
    }
    pos1 <- pData(x)[[pheno]] %in% group1
    pos2 <- pData(x)[[pheno]] %in% group2
    if (length(pos1)==0) stop(paste(pheno, " was not found in the ExpressionSet Object", sep=""), call.=FALSE)
    if (length(which(pos1))==0) stop(paste(group1, " was not found in ", pheno, sep=""), call.=FALSE)
    if (length(which(pos2))==0) stop(paste(group2, " was not found in ", pheno, sep=""), call.=FALSE)
    tmp <- marina(exprs(x)[, pos1], exprs(x)[, pos2], regulon=regulon, per=per, cores=cores, verbose=verbose)
    tmp@expset <- x
    return(tmp)
})

#' @examples
#'
#' x <- diggitClass(expset=gbmExprs, regulon=gbmTFregulon)
#' dgo <- marina(x, pheno="subtype", group1="MES", group2="PN", per=100, cores=cores)
#' dgo
#' diggitMR(dgo)[1:5]
#' @rdname marina-methods
#' @export
setMethod("marina", c(x="diggit"), function(x, pheno, group1, group2=NULL, mu=0, regulon=NULL, per=1000, cores=1, verbose=TRUE) {
    if (is.null(regulon)) regulon <- diggitRegulon(x)
    tmp <- marina(exprs(x), pheno=pheno, group1=group1, group2=group2, mu=mu, regulon=regulon, per=per, cores=cores, verbose=verbose)
    x@mr <- diggitMR(tmp)
    x@regulon <- diggitRegulon(tmp)
    return(x)
})
