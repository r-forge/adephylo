##
## Functions to obtain partitions of tips from a tree.
## For instance to obtain dummy vectors used in the orthogram.
##

############
# listTips
############
listTips <- function(tree){
    if(!require(phylobase)) stop("phylobase package is not installed")

    x <- tree
    ## conversion from phylo4 and phylo4d
    x <- as(x, "phylo4")

    ## check phylo4 object
    if (is.character(checkval <- check_phylo4(res))) stop(checkval)

    ## computations
    res <- lapply(nodeLabels(x), function(e) descendants(x, e))

    names(res) <- nodeLabels(x)

    return(res)
}
