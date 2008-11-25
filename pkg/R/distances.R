############
# distNodes
############
distNodes <- function(x, node1, node2,
                      method=c("brlength","nNodes","Abouheif","sumDD")){

    if(!require(phylobase)) stop("phylobase package is not installed")

    ## conversion from phylo, phylo4 and phylo4d
    x <- as(x, "phylo4")
    method <- match.arg(method)

    ## some checks
    if (is.character(checkval <- check_phylo4(x))) stop(checkval)
    node1 <- getnodes(x, node1)
    node2 <- getnodes(x, node2)
    if(any(is.na(c(node1,node2)))) stop("wrong node specified")
    if(node1==node2) return(0)

    ## get the path between node1 and node2
    path <- shortestPath(x, node1, node2)

    ## compute distances
    if(method=="brlength"){
        if(!hasEdgeLength(x)) stop("x does not have branch length")
        path <- c(node1, node2, path)
        path <- path[path != MRCA(x, node1, node2)]
        edge.idx <- getedges(x, path)
        res <- sum(edgeLength(x)[edge.idx], na.rm=TRUE)
        return(res)
    } # end brlength

    if(method=="nNodes"){
        res <- length(path)
        return(res)
    } # end nNodes

    if(method=="Abouheif"){
        E <- x@edge
        temp <- table(E[,1])[as.character(path)] # number of dd per node
        res <- prod(temp)
        return(res)
    } # end Abouheif

    if(method=="sumDD"){
        E <- x@edge
        temp <- table(E[,1])[as.character(path)] # number of dd per node
        res <- sum(temp)
        return(res)
 } # end sumDD

} # end distNodes





###########
# distRoot
###########
distRoot <- function(x, method=c("brlength","nNodes","Abouheif")){
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
