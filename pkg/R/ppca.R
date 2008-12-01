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
    if (is.character(chk <- check_phylo4(x))) stop("bad phylo4d object: ",chk)
    if (is.character(chk <- check_data(x))) stop("bad phylo4d object: ",chk)

    tre <- as(x, "phylo4")
    method <- match.arg(method)

    ## proximity matrix
    if(is.null(prox)){ # have to compute prox
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
    res <- dudi.pca(X, center=center, scale=scale, scannf=FALSE,nf=2)
    Upca <- as.matrix(res$c1)

    ## computations of the ppca
    X <- as.matrix(X)
    decomp <- eigen((t(X) %*% W %*% X)/N, sym=TRUE)
    U <- decomp$vectors # U: principal axes
    p <- ncol(U)
    lambda <- decomp$values

    if(scannf){ # interactive part
        barplot(lambda[1:res$rank])
        cat("Select the number of global axes: ")
        nfposi <- as.integer(readLines(n = 1))
        cat("Select the number of local axes: ")
        nfnega <- as.integer(readLines(n = 1))
    }

    nfposi <- max(nfposi, 1)
    nfnega <- max(nfnega, 0)
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

    res$cent <- res$norm <- res$co <- NULL # cleaning

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
    names(res$ls) <- scores.lab
    row.names(res$ls) <- X.rownames

    res$as <- as.data.frame(A) # PCA axes onto pPCA axes
    names(res$as) <- axes.lab
    row.names(res$as) <- paste("PCA axis", 1:nrow(A))

    res$tre <- as(tre,"phylo4") # tree

    res$prox <- W # proximity matrix

    res$call <- match.call() # call

    class(res) <- "ppca"

    return(res)
} # end ppca





#####################
# Function plot.ppca
#####################
plot.ppca <- function(x, axis=1:ncol(x$li), useLag=FALSE, ...){
    if(useLag){
        df <- as.data.frame(x$ls)
    } else{
        df <- as.data.frame(x$li)
    }

    if(any(axis < 1 | axis > ncol(x$li)) ) stop("Wrong axis specified.")
    df <- df[, axis, drop=FALSE]

    obj <- phylo4d(x$tre,df)
    args <- list(...)
    if(is.null(args$ratio.tree)){
        args$ratio.tree <- 0.5
    }
    args <- c(obj,args)
    do.call(s.phylo4d, args)
}





######################
# Function print.ppca
######################
print.ppca <- function(x, ...){
  cat("\t#############################################\n")
  cat("\t# phylogenetic Principal Component Analysis #\n")
  cat("\t#############################################\n")
  cat("class: ")
  cat(class(x))
  cat("\n$call: ")
  print(x$call)
  cat("\n$nfposi:", x$nfposi, "axis-components saved")
  cat("\n$nfnega:", x$nfnega, "axis-components saved")
  cat("\n$kept.axes: index of kept axes")

  cat("\nPositive eigenvalues: ")
  l0 <- sum(x$eig >= 0)
  cat(signif(x$eig, 4)[1:(min(5, l0))])
  if (l0 > 5)
    cat(" ...\n")
  else cat("\n")
  cat("Negative eigenvalues: ")
  l0 <- sum(x$eig <= 0)
  cat(sort(signif(x$eig, 4))[1:(min(5, l0))])
  if (l0 > 5)
    cat(" ...\n")
  else cat("\n")
  cat('\n')
  sumry <- array("", c(1, 4), list(1, c("vector", "length",
                                        "mode", "content")))
  sumry[1, ] <- c('$eig', length(x$eig), mode(x$eig), 'eigenvalues')
  class(sumry) <- "table"
  print(sumry)
  cat("\n")
  sumry <- array("", c(4, 4), list(1:4, c("data.frame", "nrow", "ncol", "content")))
  sumry[1, ] <- c("$c1", nrow(x$c1), ncol(x$c1), "principal axes: scaled vectors of traits loadings")
  sumry[2, ] <- c("$li", nrow(x$li), ncol(x$li), "principal components: coordinates of taxa ('scores')")
  sumry[3, ] <- c("$ls", nrow(x$ls), ncol(x$ls), 'lag vector of principal components')
  sumry[4, ] <- c("$as", nrow(x$as), ncol(x$as), 'pca axes onto ppca axes')

  class(sumry) <- "table"
  print(sumry)

  cat("\n$tre: a phylogeny (class phylo4)")
  cat("\n$prox: a matrix of phylogenetic proximities")

  cat("\n\nother elements: ")
  if (length(names(x)) > 16)
    cat(names(x)[17:(length(names(x)))], "\n")
  else cat("NULL\n")
} #end print.ppca





### testing
## obj <- phylo4d(read.tree(text=mjrochet$tre),mjrochet$tab)
## x@edge.length= rep(1,length(x@edge.label))
## M = cophenetic.phylo(as(x,"phylo"))
## M = 1/M
## diag(M) <- 0


## ppca1 <- ppca(obj,scannf=FALSE,nfp=1,nfn=0)

## plot(ppca1)
