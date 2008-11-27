###################
# orthobasis.phylo
###################
orthobasis.phylo <- function(x=NULL, prox=NULL,
                             method=c("brlength","nNodes","Abouheif","sumDD"), a=1){
    if(!require(phylobase)) stop("phylobase package is not installed")
    if(!require(ade4)) stop("ade4 package is not installed")

    ## handle arguments
    method <- match.arg(method)

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

    n <- nrow(W)


    ## main computation
    W0 <- bicenter.wt(W)
    decomp <- eigen(W0, sym=TRUE)
    E <- decomp$vectors # Moran eigenvectors (ME)
    ## must be re-orthogonalized
    temp <- cbind(rep(1,n) , E)
    temp <- qr.Q(qr(temp))
    E <- as.data.frame(temp[,-1])*sqrt(n)
    row.names(E) <- rownames(W)
    names(E) <- paste("ME", 1:ncol(E), sep=".")

    ## retrieve Moran' I for each ME
    f1 <- function(vec){
        res <- mean( vec * (W0 %*% vec) )
        return(res)
    }

    Ival <- apply(E, 2, f1)

    ## build output
    res <- E
    attr(res,"values") <- Ival
    attr(res,"weights") <- rep(1/n,n)
    attr(res,"call") <- match.call()
    attr(res,"class") <- c("orthobasis","data.frame")

    return(res)
} # end orthobasis.phylo





###########
# me.phylo
###########
me.phylo <- orthobasis.phylo
