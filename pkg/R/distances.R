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
