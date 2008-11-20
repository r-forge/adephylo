##
## Functions to obtain partitions of tips from a tree.
## For instance to obtain dummy vectors used in the orthogram.
##





############
# listTips
############
listTips <- function(x){
    if(!require(phylobase)) stop("phylobase package is not installed")

    ## conversion from phylo, phylo4 and phylo4d
    x <- as(x, "phylo4")

    ## check phylo4 object
    if (is.character(checkval <- check_phylo4(x))) stop(checkval)

    ## computations
    res <- lapply(1:nNodes(x), function(i) descendants(x, i))

    if(hasNodeLabels(x)) {names(res) <- nodeLabels(x)}

    return(res)
} # end listTips





############
# listNodes
############
listNodes <- function(x){
    if(!require(phylobase)) stop("phylobase package is not installed")

    ## conversion from phylo, phylo4 and phylo4d
    x <- as(x, "phylo4")

    ## check phylo4 object
    if (is.character(checkval <- check_phylo4(x))) stop(checkval)

    ## computations
    nodIdx <- nTips(x)+1
    nodIdx <- nodIdx:(nodIdx+nNodes(x)-1)
    res <- lapply(nodIdx, function(i) children(x, i))

    if(hasNodeLabels(x)) {names(res) <- nodeLabels(x)}

    return(res)
} # end listNodes





###########
# treePart
###########
treePart <- function(x){
    if(!require(phylobase)) stop("phylobase package is not installed")

    ## conversion from phylo, phylo4 and phylo4d
    x <- as(x, "phylo4")

    ## check phylo4 object
    if (is.character(checkval <- check_phylo4(x))) stop(checkval)

    n <- nTips(x)

    ## function coding one dummy vector
    fDum <- function(vec){ # vec is a vector of tip numbers
        dum <- integer(n)
        dum[vec] <- 1
        return(dum)
    }

    ## main computations
    temp <- listTips(x)
    res <- data.frame(lapply(temp,fDum))
    row.names(res) <- x@tip.label
    res <- res[,-1]

    return(res)
} # end treePart





############
# path2tips
############
path2tips <- function(x, tip1, tip2){
    if(!require(phylobase)) stop("phylobase package is not installed")

    ## conversion from phylo, phylo4 and phylo4d
    x <- as(x, "phylo4")

    ## check phylo4 object
    if (is.character(checkval <- check_phylo4(x))) stop(checkval)

    ## main computations
    t1 <- getnodes(x,tip1)
    t2 <- getnodes(x,tip2)

    comAnc <- MRCA(x, t1, t2) # common ancestor
    desComAnc <- descendants(x, comAnc, which="all")
    ancT1 <- ancestors(x, t1, which="all")
    path1 <- intersect(desComAnc, ancT1) # path: common anc -> t1

    ancT2 <- ancestors(x, t2, which="all")
    path2 <- intersect(desComAnc, ancT2) # path: common anc -> t2

    res <- union(path1, path2) # union of the path
    res <- getnodes(res)

    return(res)
} # end path2tips





###########
# distRoot
###########
distRoot <- function(x, tip){
    if(!require(phylobase)) stop("phylobase package is not installed")

    ## conversion from phylo, phylo4 and phylo4d
    x <- as(x, "phylo4")

    ## check phylo4 object
    if (is.character(checkval <- check_phylo4(x))) stop(checkval)

    ## main computations
    tip <- getnodes(x, tip)
    root <- getnodes(x, nTips(x)+1)
    ancTip <- ancestors(x, tip, which="all")

    pathNodes <- setdiff(ancTip, root) # only internal nodes, without root
    pathNodes <- c(tip, pathNodes)

    pathNodes <- getnodes(x, pathNodes)

    res <- sumEdgeLength(x, pathNodes)

    return(res)
} # end distRoot
