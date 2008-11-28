####
#### Phylogenetic ordination tools
####
#### Thibaut Jombart 2008 (jombart@biomserv.univ-lyon1.fr)

################
# Function ppca
################
ppca <- function(x, prox=NULL, method=c("patristic","nNodes","Abouheif","sumDD"), a=1,
                 center=TRUE, scale=TRUE, scannf=TRUE, nfposi=1, nfnega=0){

    ## handle arguments
    if(!require(ade4)) stop("The package ade4 is not installed.")
    if (is.character(chk <- check_phylo4d(x))) stop("Invalid phylo4d object: \n",chk)
    tre <- as(x, "phylo4")
    method <- match.arg(method)

    ## proximity matrix
    if(is.null(prox)){ # have to compute prox
        x <- as(x, "phylo4")
        if (is.character(checkval <- check_phylo4(x))) stop(checkval)
        W <- proxTips(x, tips="all", method=method, a=a, normalize="row", symmetric=TRUE)
    } else { # prox is provided
        W <- as.matrix(prox)
        if(!is.matrix(W)) stop("W is not a matrix")
        if(ncol(W) != nrow(W)) stop("W is not a square matrix")
        diag(W) <- 0
        W <- 0.5 * (t(W) + W) # re-symmetrization
    }

    N <- nTips(x)

    ## data matrix X
    X <- tdata(x, which="tip")
    X.colnames <- names(X)
    X.rownames <- row.names(X)
    temp <- sapply(X, is.numeric)
    if(!all(temp)) {
        warning(paste("non-numeric data are removed:", X.colnames[!temp]))
        X <- X[,temp]
        X.colnames <- X.colnames[!temp]
        X.rownames <- X.rownames[!temp]
    }

    ## replace NAs
    f1 <- function(vec){
        if(any(is.na(vec))) {
            warning("Replacing missing values (NA) by mean values")
            m <- mean(vec,na.rm=TRUE)
            vec[is.na(vec)] <- m
        }

        return(vec)
    }

    X <- as.data.frame(lapply(X, f1))
    X <- scalewt(X, center=center, scale=scale) # centring/scaling of traits

    ##


    ## main computation ##

    ## make a skeleton of dudi
    res <- dudi.pca(X, center=FALSE, scale=FALSE, scannf=FALSE,nf=2)
    Upca <- as.matrix(res$c1)

    ## computations of the ppca
    X <- as.matrix(X)
    decomp <- eigen((t(X) %*% W %*% X)/n, sym=TRUE)
    U <- decomp$vectors # U: principal axes
    p <- ncol(U)
    lambda <- U$values

    if(scannf){ # interactive part
        barplot(eig[1:rank])
        cat("Select the number of global axes: ")
        nfposi <- as.integer(readLines(n = 1))
        cat("Select the number of local axes: ")
        nfnega <- as.integer(readLines(n = 1))
    }

    nfposi <- max(nfposi, 1)
    nfnega <- max(nfposi, 0)
    posi.idx <- 1:nfposi
    if(nfnega<1) {
        nega.idx <- NULL
    } else {
        nega.idx <- (p-nfnega+1):p
    }

    axes.idx <- unique(c(posi.idx, nega.idx)) # index of kept axes
    U <- U[, axes.idx, drop=FALSE]

    S <- X %*% U # S: scores (=princ. components)
    LS <- W %*% S # LS: lagged scores
    A <- t(Upca) %*% U # A: pca princ. axes onto ppca princ. axes.

    ## build the output
    axes.lab <- paste("PA",axes.idx, sep="")
    scores.lab <- paste("PC",axes.idx, sep="")

    res$eig <- lambda # eigenvalues
    res$nf <- NULL
    res$nfposi <- nfposi
    res$nfnega <- nfnega
    res$kept.axes <- axes.idx

    res$c1 <- as.data.frame(U) # principal axes
    names(res$c1) <- axes.lab
    row.names(res$c1) <- X.colnames

    res$li <-  as.data.frame(S) # scores (princ. components)
    names(res$li) <- scores.lab
    row.names(res$li) <- X.rownames

    res$ls <-  as.data.frame(S) # lagged scores
    names(res$li) <- scores.lab
    row.names(res$li) <- X.rownames

    res$as <- as.data.frame(A) # PCA axes onto pPCA axes
    names(res$li) <- axes.lab
    row.names(res$li) <- paste("PCA axis", 1:nrow(A))

    res$tre <- as(tre,"phylo4") # tree

    res$prox <- prox # proximity matrix

    class(res) <- "ppca"

    return(res)
} # end ppca



#####################
# Function plot.ppca
#####################
plot.ppca <- function(x,laged=FALSE, ...){
    if(laged){
        df <- as.data.frame(x$ls)
    } else{
        df <- as.data.frame(x$li)
    }

    obj <- phylo4d(x$tre,df)
    args <- list(...)
    if(is.null(args$ratio.tree)){
        args$ratio.tree <- 0.5
    }
    args <- c(obj,args)
    do.call(plot, args)
}


### testing
## obj <- phylo4d(read.tree(text=mjrochet$tre),mjrochet$tab)
## x@edge.length= rep(1,length(x@edge.label))
## M = cophenetic.phylo(as(x,"phylo"))
## M = 1/M
## diag(M) <- 0


## ppca1 <- ppca(obj,scannf=FALSE,nfp=1,nfn=0)

## plot(ppca1)
