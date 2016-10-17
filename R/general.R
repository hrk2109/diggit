#' @import methods
#' @import Biobase
#' @importFrom parallel mclapply
#' @importFrom ks Hpi
#' @importFrom ks kde
#' @importFrom viper viper
#' @importFrom viper msviper
#' @importFrom viper rowTtest
#' @importFrom viper ttestNull
NULL

#' Mutual information
#' 
#' This function estimates the mutual information between continuous variables using a fix bandwidth implementation
#' 
#' @param x Numeric vector or matrix
#' @param y Optional numeric vector or matrix
#' @param per Integer indicating the number of permutations to compute p-values
#' @param pairwise Logical, wether columns of x and y should be compared in a pairwise maner. x and y must have the same number of columns
#' @param bw Integer indicating the grid size for integrating the joint probability density
#' @param cores Integer indicating the number of cores to use (1 for Windows-based systems)
#' @param verbose Logical, whether progression bars should be shown
#' @return Numeric value, vector or matrix of results
#' @description This function estimates the mutual information between x and y given both are numeric vectors, between the columns of x if it is a numeric matrix, or between the columns of x and y if both are numeric matrixes
#' @export
#' @examples
#' x <- seq(0, pi, length=100)
#' y <- 5*sin(x)+rnorm(100)
#' cor.test(x, y)
#' mutualInfo(x, y, per=100)
mutualInfo <- function(x, y=NULL, per=0, pairwise=FALSE, bw=100, cores=1, verbose=TRUE){
    pb <- NULL
    if (verbose) message("Estimating kernel bandwidth...")
    if (is.matrix(x)) {
        if (ncol(x)==1) x <- x[, 1]
    }
    if (is.matrix(y)) {
        if (ncol(y)==1) y <- y[, 1]
    }
    if (is.matrix(x)) {
        x <- qnorm(apply(x, 2, rank)/(nrow(x)+1))
        if (is.null(y)) {
            h <- min(choose(ncol(x), 2), bw)
            tmp <- combn(min(100, ncol(x)), 2)
            tmp <- tmp[, sample(ncol(tmp), h)]
            if (cores==1) {
                h <- apply(tmp, 2, function(i, x) {
                    Hpi(x[, i])
                }, x=x)
            }
            else {
                h <- mclapply(1:ncol(tmp), function(i, tmp, x) {
                    Hpi(x[, tmp[, i]])
                }, tmp=tmp, x=x, mc.cores=cores)
                h <- sapply(h, function(x) x)
            }
            tmp <- apply(h, 1, rank)-(ncol(h)/2)
            tmp <- which.min(abs(rowMeans(tmp)))
            h <- matrix(h[, tmp], 2, 2)
            if (verbose) message("Computing MI...")
            if (cores==1) {
                if (verbose) {
                    pb <- txtProgressBar(max=ncol(x)-1, style=3)
                }
                tmp <- lapply(1:(ncol(x)-1), function(i, x, h, pb) {
                    if (!is.null(pb)) setTxtProgressBar(pb, i)
                    apply(filterColMatrix(x, (i+1):ncol(x)), 2, function(x1, x2, h) {
                        internalMIfb(x1, x2, h=h)
                    }, x2=x[, i], h=h)
                }, x=x, h=h, pb=pb)
            }
            else {
                tmp <- mclapply(1:(ncol(x)-1), function(i, x, h) {
                    apply(filterColMatrix(x, (i+1):ncol(x)), 2, function(x1, x2, h) {
                        internalMIfb(x1, x2, h=h)
                    }, x2=x[, i], h=h)
                }, x=x, h=h, mc.cores=cores)
            }
            tmp <- unlist(tmp, use.names=FALSE)
            mi <- matrix(NA, ncol(x), ncol(x))
            mi[lower.tri(mi)] <- tmp
            mi <- t(mi)
            mi[lower.tri(mi)] <- tmp
            rownames(mi) <- colnames(mi) <- colnames(x)
        }
        else {
            if (is.matrix(y)) {
                y <- qnorm(apply(y, 2, rank)/(nrow(y)+1))
                if (ncol(x)==ncol(y) & pairwise) {
                    h <- min(ncol(x), bw)
                    tmp <- sample(ncol(x), h)
                    if (cores==1) {
                        h <- sapply(tmp, function(i, x, y) {
                            Hpi(cbind(x[, i], y[, i]))
                        }, x=x, y=y)
                    }
                    else {
                        h <- mclapply(tmp, function(i, x, y) {
                            Hpi(cbind(x[, i], y[, i]))
                        }, x=x, y=y, mc.cores=cores)
                        h <- sapply(h, function(x) x)
                    }
                    tmp <- apply(h, 1, rank)-(ncol(h)/2)
                    tmp <- which.min(abs(rowMeans(tmp)))
                    h <- matrix(h[, tmp], 2, 2)
                    if (verbose) message("Computing MI...")
                    if (cores==1) {
                        if (verbose) {
                            pb <- txtProgressBar(max=ncol(x), style=3)
                        }                        
                        mi <- sapply(1:ncol(x), function(i, x, y, h, pb) {
                            if (!is.null(pb)) setTxtProgressBar(pb, i)
                            internalMIfb(x[, i], y[, i], h=h)
                        }, x=x, y=y, h=h, pb=pb)
                    }
                    else {
                        mi <- mclapply(1:ncol(x), function(i, x, y, h) {
                            internalMIfb(x[, i], y[, i], h=h)
                        }, x=x, y=y, h=h, mc.cores=cores)
                        mi <- unlist(mi, use.names=FALSE)
                    }
                    names(mi) <- colnames(x)
                }
                else {
                    h <- min(ncol(x)*ncol(y), bw)
                    tmp <- rbind(rep(1:ncol(x), rep(ncol(y), ncol(x))), rep(1:ncol(y), ncol(x)))[, sample(h)]
                    if (cores==1) {
                        h <- apply(tmp, 2, function(i, x, y) {
                            Hpi(cbind(x[, i[1]], y[, i[2]]))
                        }, x=x, y=y)
                    }
                    else {
                        h <- mclapply(1:ncol(tmp), function(i, tmp, x, y) {
                            i <- tmp[, i]
                            Hpi(cbind(x[, i[1]], y[, i[2]]))
                        }, x=x, y=y, tmp=tmp, mc.cores=cores)
                        h <- sapply(h, function(x) x)
                    }
                    tmp <- apply(h, 1, rank)-(ncol(h)/2)
                    tmp <- which.min(abs(rowMeans(tmp)))
                    h <- matrix(h[, tmp], 2, 2)
                    if (verbose) message("Computing MI...")
                    if (cores==1) {
                        if (verbose) {
                            pb <- txtProgressBar(max=ncol(y), style=3)
                        }
                        mi <- sapply(1:ncol(y), function(i, x, y, h, pb) {
                            y <- y[, i]
                            if (!is.null(pb)) setTxtProgressBar(pb, i)
                            apply(x, 2, function(x, y, h) {
                                internalMIfb(x, y, h=h)
                            }, y=y, h=h)
                        }, x=x, y=y, h=h, pb=pb)
                    }
                    else {
                        mi <- mclapply(1:ncol(y), function(i, x, y, h) {
                            y <- y[, i]
                            apply(x, 2, function(x, y, h) {
                                internalMIfb(x, y, h=h)
                            }, y=y, h=h)
                        }, x=x, y=y, h=h, mc.cores=cores)
                        mi <- sapply(mi, function(x) x)
                    }
                    colnames(mi) <- colnames(y)
                    rownames(mi) <- colnames(x)
                }
            }
            else {
                y <- qnorm(rank(y)/(length(y)+1))
                h <- min(ncol(x), bw)
                tmp <- sample(ncol(x), h)
                if (cores==1) {
                    h <- apply(x[, tmp], 2, function(x, y) {
                        Hpi(cbind(x, y))
                    }, y=y)
                }
                else {
                    h <- mclapply(tmp, function(i, x, y) {
                        x <- x[, i]
                        Hpi(cbind(x, y))
                    }, x=x, y=y, mc.cores=cores)
                    h <- sapply(h, function(x) x)
                }
                tmp <- apply(h, 1, rank)-(ncol(h)/2)
                tmp <- which.min(abs(rowMeans(tmp)))
                h <- matrix(h[, tmp], 2, 2)
                if (verbose) message("Computing MI...")
                if (cores==1) {
                    if (verbose) {
                        pb <- txtProgressBar(max=ncol(x), style=3)
                    }
                    mi <- sapply(1:ncol(x), function(i, x, y, h, pb) {
                        x <- x[, i]
                        if (!is.null(pb)) setTxtProgressBar(pb, i)
                        internalMIfb(x, y, h=h)
                    }, y=y, h=h, x=x, pb=pb)
                }
                else {
                    mi <- mclapply(1:ncol(x), function(i, x, y, h) {
                        x <- x[, i]
                        internalMIfb(x, y, h=h)
                    }, y=y, h=h, x=x, mc.cores=cores)
                    mi <- unlist(mi, use.names=FALSE)
                }
                names(mi) <- colnames(x)
            }
        }
    }
    else {
        if (is.null(y)) stop("y is required when x is a numeric vector", call.=FALSE)
        x <- qnorm(rank(x)/(length(x)+1))
        if (is.matrix(y)) {
            y <- qnorm(apply(y, 2, rank)/(nrow(y)+1))
            h <- min(ncol(y), bw)
            tmp <- sample(ncol(y), h)
            if (cores==1) {
                h <- apply(y[, tmp], 2, function(y, x) {
                    Hpi(cbind(x, y))
                }, x=x)
            }
            else {
                h <- mclapply(tmp, function(i, y, x) {
                    y <- y[, i]
                    Hpi(cbind(x, y))
                }, x=x, y=y, mc.cores=cores)
                h <- sapply(h, function(x) x)
            }
            tmp <- apply(h, 1, rank)-(ncol(h)/2)
            tmp <- which.min(abs(rowMeans(tmp)))
            h <- matrix(h[, tmp], 2, 2)
            if (verbose) message("Computing MI...")
            if (cores==1) {
                if (verbose) {
                    pb <- txtProgressBar(max=ncol(y), style=3)
                }
                mi <- sapply(1:ncol(y), function(i, x, y, h, pb) {
                    y <- y[, i]
                    if (!is.null(pb)) setTxtProgressBar(pb, i)
                    internalMIfb(x, y, h=h)
                }, y=y, h=h, x=x, pb=pb)
            }
            else {
                mi <- mclapply(1:ncol(y), function(i, x, y, h) {
                    y <- y[, i]
                    internalMIfb(x, y, h=h)
                }, y=y, h=h, x=x, mc.cores=cores)
                mi <- unlist(mi, use.names=FALSE)
            }
            names(mi) <- colnames(y)
        }
        else {
            y <- qnorm(rank(y)/(length(y)+1))
            h <- Hpi(cbind(x, y))
            mi <- internalMIfb(x, y, h=h)
        }
    }
    if (!is.null(pb)) message("\n")    
    if (per>0) {
        if (is.matrix(x)) dnull <- qnorm((1:nrow(x))/(nrow(x)+1))
        else dnull <- qnorm((1:length(x))/(length(x)+1))
        if (verbose) message("Computing MI null model...")
        if (cores==1) {
            if (verbose) {
                pb <- txtProgressBar(max=per, style=3)
            }
            nullmi <- sapply(1:per, function(i, dnull, h, pb) {
                if (!is.null(pb)) setTxtProgressBar(pb, i)
                internalMIfb(sample(dnull), sample(dnull), h=h)
            }, dnull=dnull, h=h, pb=pb)
        }
        else {
            nullmi <- mclapply(1:per, function(i, dnull, h) {
                internalMIfb(sample(dnull), sample(dnull), h=h)
            }, dnull=dnull, h=h, mc.cores=cores)
            nullmi <- sapply(nullmi, function(x) x)
        }
        tmp <- aecdf(nullmi)
        tmp <- lapply(as.vector(mi), function(x, tmp) tmp(x, alternative="greater"), tmp=tmp)
        z <- sapply(tmp, function(x) x$nes)
        p <- sapply(tmp, function(x) x$p.value)
        dim(z) <- dim(p) <- dim(mi)
        if (!is.null(nrow(mi))) {
            colnames(z) <- colnames(p) <- colnames(mi)
            rownames(z) <- rownames(p) <- rownames(mi)
        }
        else names(z) <- names(p) <- names(mi)
        return(list(mi=mi, z=z, p.value=p))
    }
    if (is.null(pb)) message("\n")
    return(mi)
}

internalMIfb <- function(x, y, h, n=10) {
    x1 <- rep(seq(min(x), max(x), length=n), rep(n, n))
    y1 <- rep(seq(min(y), max(y), length=n), n)
    k1 <- dnorm(x1)
    k2 <- dnorm(y1)
    k12 <- kde(cbind(x, y), h, xmin=c(min(x), min(y)), xmax=c(max(x), max(y)), eval.points=cbind(x1, y1))$estimate
    sum(k12*log10(k12/k1/k2), na.rm=T)/n^2*(max(x)-min(x))
}

#' Correlation test
#' 
#' This function computes correlation and associated p-values
#' 
#' @param x Numeric vector or matrix
#' @param y Optional numeric vector or matrix
#' @param method Character string indicating the correlation method
#' @param pairwise Logical, wether columns of x and y should be compared in a pairwise manner. x and y must have the same number of columns
#' @return Numeric value, vector or matrix of results
#' @description This function computes the correlation between x and y given both are numeric vectors, between the columns of x if it is a numeric matrix, or between the columns of x and y if both are numeric matrixes
#' @export
#' @examples
#' x <- seq(0, 10, length=50)
#' y <- x+rnorm(length(x), sd=2)
#' correlation(x, y)
correlation <- function(x, y=NULL, method=c("pearson", "spearman", "kendall"), pairwise=FALSE){
    method <- match.arg(method)
    if (is.matrix(x)) {
        n <- nrow(x)
        if (is.matrix(y)) {
            if (ncol(x)==ncol(y) & pairwise) {
                mi <- sapply(1:ncol(x), function(i, x, y) {
                    cor(x[, i], y[, i], method=method, use="pairwise.complete.obs")
                }, x=x, y=y)
                names(mi) <- colnames(x)
                if ((sum(is.na(x))+sum(is.na(y)))>0) {
                    n <- sapply(1:ncol(x), function(i, x, y) {
                        length(which(rowSums(!is.na(cbind(x[, i], y[, i])))==2))
                    }, x=x, y=y)
                }
            }
            else {
                mi <- cor(x, y, method=method, use="pairwise.complete.obs")
                if ((sum(is.na(x))+sum(is.na(y)))>0) {
                    n <- sapply(1:ncol(y), function(i, x, y) {
                        sapply(1:ncol(x), function(ii, x, y, i) {
                            length(which(rowSums(!is.na(cbind(x[, ii], y[, i])))==2))
                        }, x=x, y=y, i=i)
                    }, x=x, y=y)
                }
            }
        }
        else {
            y <- matrix(y, length(y), 1, dimnames=list(names(y), "y"))
            mi <- cor(x, y, method=method, use="pairwise.complete.obs")
            if ((sum(is.na(x))+sum(is.na(y)))>0) {
                n <- sapply(1:ncol(x), function(i, x, y) {
                    length(which(rowSums(!is.na(cbind(x[, i], y)))==2))
                }, x=x, y=y)
            }
        }
    }
    else {
        n <- length(x)
        if (is.null(y)) stop("y is required when x is a numeric vector", call.=FALSE)
        if (is.matrix(y)) {
            x <- matrix(x, length(x), 1, dimnames=list(names(x), "x"))
        }
        else {
            x <- matrix(x, length(x), 1, dimnames=list(names(x), "x"))
            y <- matrix(y, length(x), 1, dimnames=list(names(x), "x"))
        }
        mi <- cor(x, y, method=method, use="pairwise.complete.obs")
        if ((sum(is.na(x))+sum(is.na(y)))>0) {
            n <- sapply(1:ncol(y), function(i, x, y) {
                length(which(rowSums(!is.na(cbind(x, y[, i])))==2))
            }, x=x, y=y)
        }
    }
    tmp <- cortest(mi, n)
    return(list(r=mi, z=tmp$z, p.value=tmp$p.value))
}

cortest <- function (r, n, alternative = "two.sided") 
{
    z <- log((1 + r)/(1 - r))/2 * sqrt(n - 3)
    switch(match.arg(alternative, c("two.sided", "less", "greater")), 
           two.sided = {
               p <- pnorm(abs(z), lower.tail = F) * 2
           }, less = {
               p <- pnorm(z, lower.tail = T)
           }, greater = {
               p <- pnorm(z, lower.tail = F)
           })
    return(list(z = z, p.value = p))
}

#' Approximate empirical commulative distribution function
#'
#' This function generates an empirical null model that computes a normalized statistics and p-value
#' 
#' @param dnull Numerical vector representing the null model
#' @param symmetric Logical, whether the distribution should betreated as symmetric around zero and only one tail should be approximated
#' @param n Integer indicating the number of points to evaluate the empirical cummulative probability function
#' @return function with two parameters, \code{x} and \code{alternative}

aecdf <- function(dnull, symmetric=FALSE, n=100) {
    dnull <- dnull[is.finite(dnull)]
    if (symmetric) {
        iqr <- quantile(abs(dnull), c(.5, 1-5/length(dnull)))
        epd <- ecdf(abs(dnull))
        a <- list(x=knots(epd), y=epd(knots(epd)))
        fit <- lm(y~0+x, data=list(x=a$x[length(a$x)-(15:4)]-iqr[2], y=log(1-epd(iqr[2]))-log(1-a$y[length(a$x)-(15:4)])))
        val <- seq(0, iqr[2], length=n)
        pd <- approxfun(val, epd(val), method="linear", yleft=0, rule=2)
        dnull <- function(x, alternative=c("two.sided", "greater", "less")) {
            alternative <- match.arg(alternative)
            x1 <- abs(x)
            p <- exp(log(1-pd(iqr[2]))-predict(fit, list(x=x1-iqr[2])))
            p[!is.finite(p)] <- 1
            p <- p * (x1>iqr[2]) + (1-pd(x1)) * (x1<=iqr[2])
            nes <- qnorm(p/2, lower.tail=F)*sign(x)
            switch(alternative,
                   two.sided={p <- p},
                   greater={p <- p/2; p[x<0] <- 1-p[x<0]},
                   less={p <- p/2; p[x>0] <- 1-p[x>0]}
            )
            names(nes) <- names(p) <- names(x)
            list(nes=nes, p.value=p)
        }
        return(dnull)
    }
    iqr <- quantile(dnull, c(5/length(dnull), .5, 1-5/length(dnull)))
    epd <- ecdf(dnull)
    a <- list(x=knots(epd), y=epd(knots(epd)))
    fit1 <- lm(y~0+x, data=list(x=a$x[5:14]-iqr[1], y=log(epd(iqr[1]))-log(a$y[5:14])))
    fit2 <- lm(y~0+x, data=list(x=a$x[length(a$x)-(15:4)]-iqr[3], y=log(1-epd(iqr[3]))-log(1-a$y[length(a$x)-(15:4)])))
    val <- seq(iqr[1], iqr[3], length=n)
    pd <- approxfun(val, epd(val), method="linear", rule=2)
    dnull <- function(x, alternative=c("two.sided", "greater", "less")) {
        alternative <- match.arg(alternative)
        p1 <- exp(log(pd(iqr[1]))-predict(fit1, list(x=x-iqr[1])))
        p2 <- exp(log(1-pd(iqr[3]))-predict(fit2, list(x=x-iqr[3])))
        p1[!is.finite(p1)] <- 1
        p2[!is.finite(p2)] <- 1
        p <- p1*(x<iqr[1]) + p2*(x>iqr[3]) + pd(x)*(x>=iqr[1] & x<iqr[2]) + (1-pd(x))*(x>=iqr[2] & x<=iqr[3])
        nes <- qnorm(p, lower.tail=F)*sign(x-iqr[2])
        switch(alternative,
               two.sided={p <- p*2},
               greater={p[x<iqr[2]] <- 1-p[x<iqr[2]]},
               less={p[x>=iqr[2]] <- 1-p[x>=iqr[2]]}
        )
        names(nes) <- names(p) <- names(x)
        list(nes=nes, p.value=p)
    }
    return(dnull)
}

jaccard <- function(x, y=NULL) {
    if (is.matrix(x)) {
        tmp <- lapply(1:(ncol(x)-1), function(i, x) {
            apply(filterColMatrix(x, (i+1):ncol(x)), 2, function(x1, x2) {
                jaccard(x1, x2)
            }, x2=x[, i])
        }, x=x)
        tmp <- unlist(tmp, use.names=FALSE)
        ji <- matrix(1, ncol(x), ncol(x))
        ji[lower.tri(ji)] <- tmp
        ji <- t(ji)
        ji[lower.tri(ji)] <- tmp
        rownames(ji) <- colnames(ji) <- colnames(x)
        ji[is.na(ji)] <- 0
        return(ji)
    }
    sum(x&y)/sum(x|y)
}

fet <- function(x, y=NULL) {
    if (is.matrix(x)) {
        tmp <- lapply(1:(ncol(x)-1), function(i, x) {
            apply(filterColMatrix(x, (i+1):ncol(x)), 2, function(x1, x2) {
                fet(x1, x2)
            }, x2=x[, i])
        }, x=x)
        tmp <- unlist(tmp, use.names=FALSE)
        ji <- matrix(1, ncol(x), ncol(x))
        ji[lower.tri(ji)] <- tmp
        ji <- t(ji)
        ji[lower.tri(ji)] <- tmp
        rownames(ji) <- colnames(ji) <- colnames(x)
        ji[is.na(ji)] <- 0
        return(ji)
    }
    fisher.test(x, y, alternative="greater")$p.value
}

distMode <- function (x, adj = 1) 
{
    tmp <- density(x, adjust = adj)
    return(tmp$x[which.max(tmp$y)])
}

plothm  <- function(x, color=c("cornflowerblue","salmon"), gama=1, grid=T, scmax=0, ...) {
    if (scmax==0) scmax <- max(abs(x), na.rm=T)
    pos <- which(abs(x) > scmax)
    if (length(pos)>0) x[pos] <- scmax*sign(x[pos])
    x <- abs(x/scmax)^gama*sign(x)
    x <- filterRowMatrix(x, nrow(x):1)
    color <- rgb2hsv(col2rgb(color))
    satval <- color[3,]
    color <- color[1,]
    x1 <- x
    x1[is.na(x1)] <- 0
    coli <- hsv(ifelse(x1<0, color[1], color[2]), abs(x1), 1)
    coli[is.na(x)] <- hsv(0, 0, .5)
    image(1:ncol(x), 1:nrow(x), t(matrix(1:(ncol(x)*nrow(x)), nrow(x), ncol(x))), col=coli, ylab="", xlab="", axes=F, ...)
    box()
    if (grid) grid(ncol(x), nrow(x), col="black", lty=1)
}

filterRowMatrix <- function (x, filter) {
    if (is.logical(filter)) 
        largo <- length(which(filter))
    else largo <- length(filter)
    matrix(x[filter, ], largo, ncol(x), dimnames = list(rownames(x)[filter], colnames(x)))
}

filterColMatrix <- function(x, filter) t(filterRowMatrix(t(x), filter))
