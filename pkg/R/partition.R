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
    nodIdx <- nTips(x)+1
    nodIdx <- nodIdx:(nodIdx+nNodes(x)-1)
    res <- lapply(nodIdx, function(i) descendants(x, i))

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
