#' Diggit plot
#' 
#' This function generate plots for the diggit conditional analysis
#' 
#' @param x Diggit class object
#' @param mr Optional vector of character strings indicating the MR names
#' @param cluster Optional vector of cluster names
#' @param sub Optional sub-title for the plot
#' @param ... Additional parameters to pass to the plot function
#' @return Nothing, plots are generated in the default output device
#' @export
#' @rdname diggit-plot
#' @aliases plot
#' @docType methods
#' @examples
#' data(gbm.expression, package="diggitdata")
#' data(gbm.cnv, package="diggitdata")
#' data(gbm.aracne, package="diggitdata")
#' dobj <- diggitClass(expset=gbmExprs, cnv=gbmCNV, regulon=gbmTFregulon)
#' dobj <- fCNV(dobj)
#' dobj <- aqtl(dobj, mr=c("CEBPD", "STAT3"), fcnv.adjust="fdr", verbose=FALSE)
#' dobj <- conditional(dobj, pheno="subtype", group1="MES", group2="PN", mr="STAT3", verbose=FALSE)
#' plot(dobj, cluster="3")
setMethod("plot", "diggit", function(x, mr=NULL, cluster=NULL, sub=NULL, ...) {
    x <- diggitConditional(x)
    x1 <- lapply(x, function(x) {
        tmp <- c(x$deleted, x$amplified)
        names(tmp) <- 1:length(tmp)
        return(tmp)
    })
    if (is.null(mr)) mr <- names(x)
    if (length(which(mr %in% names(x)))==0) stop("No conditional results for the selected MRs were found", call.=FALSE)
    for (i in mr) {
        if (is.null(cluster)) cluster <- names(x1[[i]])
        if (length(which(cluster %in% names(x1[[i]])))==0) stop("No conditional results for the selected clusters were found", call.=FALSE)
        orig <- x[[i]]$original
        integ <- x[[i]]$integrated
        for (ii in cluster) {
            d <- x1[[i]][[ii]]
            if (length(d)>3) {
                mtmp <- sub
                if (is.null(mtmp)) mtmp <- paste(i, " cluster ", ii, sep="")
                diag(d) <- orig[match(colnames(d), names(orig))]
                d <- cbind(d, Integ=integ[match(rownames(d), names(integ))])
                d <- -log10(d)
                plothm(d, sub=mtmp, font.sub=2, cex.sub=1.2, ...)
                abline(v=ncol(d)-.5, lwd=4)
                abline(v=ncol(d)-.5, lwd=2, col="white")
                text(rep(1:ncol(d), rep(nrow(d), ncol(d))), rep(nrow(d):1, ncol(d)), round(d, 1), cex=.8)
                axis(3, 1:ncol(d), colnames(d), las=2, tick=F, line=-.5)
                axis(2, nrow(d):1, rownames(d), tick=FALSE, line=-.5, las=2)
                axis(1, ncol(d)/2+.5, "Condition", tick=FALSE, font=2, cex=1.2)
                axis(4, nrow(d)/2+.5, "Association", tick=FALSE, font=2, cex=1.2)
            }
        }
    }
})

