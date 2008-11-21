############
# distNodes
############
distNodes <- function(x, node1, node2, method=c("brlength","nAncestors","Abouheif")){
    if(!require(phylobase)) stop("phylobase package is not installed")

    ## conversion from phylo, phylo4 and phylo4d
    x <- as(x, "phylo4")
    method <- match.arg(method)

    ## some checks
    if (is.character(checkval <- check_phylo4(x))) stop(checkval)
    t1 <- getnodes(x, node1)
    t2 <- getnodes(x, node2)
    if(any(is.na(c(t1,t2)))) stop("wrong node specified")
    if(t1==t2) return(0)

    ## get the path between node1 and node2
    path <- shortestPath(x, node1, node2)

    ## compute distances
    if(method=="brlength"){
        if(!hasEdgeLength(x)) stop("x does not have branch length")
        path <- c(node1, node2, path)
        edge.idx <- getedges(x, path)
        res <- sum(edgeLength(x)[edge.idx])
        return(res)
    } # end brlength

    if(method=="nAncestors"){
        res <- length(path)
        return(res)
    } # end nAncestors

    if(method=="Abouheif"){

    } # end Abouheif



    return(res)
} # end distNodes




###########
# distRoot
###########
distRoot <- function(x, method=c("brlength","nAncestors","Abouheif")){
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
