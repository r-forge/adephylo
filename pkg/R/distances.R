############
# distNodes
############
distTips <- function(x, tips="all",
                      method=c("brlength","nNodes","Abouheif","sumDD")){

    if(!require(phylobase)) stop("phylobase package is not installed")

    ## handle arguments
    x <- as(x, "phylo4")
    method <- match.arg(method)
    tips <- getnodes(x, tips)
    N <- nTips(x)
    if(tips="all") { tips <- 1:N }

    ## some checks
    if (is.character(checkval <- check_phylo4(x))) stop(checkval)
    if(any(is.na(tips))) stop("wrong tips specified")

    ## create all couples of observations
    findAllPairs <- function(vec){
        res <- list(i=NULL,j=NULL)
        k <- 0
        for(i in 1:(length(vec)-1)){
            for(j in 2:length(vec)){
                k <- k+1
                res[[1]][k] <- i
                res[[2]][k] <- j
            }
        }
        res <- data.frame(res)
        return(res)
    }

    allPairs <- findAllPairs(tips) # this contains all possible pairs of tips

    ## get the shortest path between all pairs of tips
    if(method != "brlength") {
        allPath <- sp.tips(x, allPairs$i, allPairs$j, useTipNames=TRUE, quiet=TRUE)
    } else {
        allPath <- sp.tips(x, allPairs$i, allPairs$j, useTipNames=TRUE, quiet=TRUE, include.mrca=FALSE)
    }

    ## compute distances
    if(method=="brlength"){
        if(!hasEdgeLength(x)) stop("x does not have branch length")
        ## add tip1 and tip2 to the paths, so that these edges are counted
        allPath <- as.data.frame(allPath)
        tip1 <- allPairs$i
        tip2 <- allPairs$j
        for(i in 1:length(allPath)){
            allPath[[i]] <- c(allPath[[i]], tip1, tip2)
        }

        edge.idx <- lapply(allPath, function(e) getedges(x, e) )
        allEdgeLength <- edgeLength(x)
        res <- lapply(edge.idx, function(idx) sum(allEdgeLength[idx], na.rm=TRUE) )
        return(res)
    } # end brlength

    if(method=="nNodes"){
        res <- length(allPath)
        return(res)
    } # end nNodes

    if(method=="Abouheif"){
        E <- x@edge
        temp <- table(E[,1])[as.character(allPath)] # number of dd per node
        res <- prod(temp)
        return(res)
    } # end Abouheif

    if(method=="sumDD"){
        E <- x@edge
        temp <- table(E[,1])[as.character(allPath)] # number of dd per node
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
