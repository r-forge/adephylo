##########
# sp.tips
##########
sp.tips <- function(x, tip1, tip2, useTipNames=FALSE, quiet=FALSE, include.mrca=TRUE){
    if(!require(phylobase)) stop("phylobase package is not installed")

    ## conversion from phylo, phylo4 and phylo4d
    x <- as(x, "phylo4")

    ## some checks
    if (is.character(checkval <- check_phylo4(x))) stop(checkval)
    t1 <- getnodes(x, tip1)
    t2 <- getnodes(x, tip2)
    if(any(is.na(c(t1,t2)))) stop("wrong tip specified")
    if(any(c(t1,t2) > nTips(x))) stop("specified nodes are internal nodes")
    if(length(t1) != length(t2)) { # recycle tip1 and tip2
        maxLength <- max(length(t1), length(t2))
        t1 <- rep(t1, length.out=maxLength)
        t2 <- rep(t2, length.out=maxLength)
    }
    toRemove <- (t1==t2)
    if(sum(toRemove)>0) {
        t1 <- t1[!toRemove]
        t2 <- t2[!toRemove]
        if(length(t1)==0) stop("tip1 and tip2 are the same vectors")
        if(!quiet) warning("tip1 and tip2 are sometimes the same; erasing these cases")
    }


    ## some global variables
    root <- rootNode(x)
    N <- nTips(x)
    E <- x@edge
    allTips <- unique(c(t1,t2))


    ## function tipToRoot
    tipToRoot <- function(E, tip){
        path <- NULL
        curNode <- tip
        while(curNode != root){
            curNode <- E[(curNode==E[,2]),1] # one node <- its ancestor
            path <- c(path, curNode)
        } # end while

        path <- getnodes(x, path)
        return(path)
    } # end tipToRoot


    ## function pathTwoTips (takes two path-to-root as args)
    pathTwoTips <- function(path1, path2){
        cpath <- c(path1, path2)
        temp <- factor(cpath, level=unique(cpath))
        CA <- temp[table(temp)==2][1] # most recent common ancestor (MRCA)
        CA <- as.integer(as.character(CA)) # retrieve integer type
        path1 <- path1[1:(which(path1==CA))] # cut path1 after MRCA (keep MRCA)
        temp <- which(path2==CA)
        if(temp==1) return(path1)
        path2 <- path2[1:(temp-1)] # cut path2 after MRCA (erase MRCA)
        return(c(path1,path2))
    } # end pathTwoTips


    pathTwoTips.no.mrca <- function(path1, path2){
        cpath <- c(path1, path2)
        temp <- intersect(path1, path2)
        res <- setdiff(cpath, temp)
        return(res)
    } # end pathTwoTips



    ## main computations
    allPathToRoot <- lapply(allTips, function(i) tipToRoot(E, i))
    names(allPathToRoot) <- allTips

    allPath1 <- allPathToRoot[as.character(t1)]
    allPath2 <- allPathToRoot[as.character(t2)]

    if(include.mrca) {
        res <- lapply(1:length(allPath1), function(i) pathTwoTips(allPath1[[i]], allPath2[[i]]) )
    } else {
        res <- lapply(1:length(allPath1), function(i) pathTwoTips.no.mrca(allPath1[[i]], allPath2[[i]]) )
        temp.names <- names(res)
        temp <- sapply(res, function(vec) length(vec)>0)
        res[temp] <- lapply(res[temp], function(vec) getnodes(x, vec) ) # name the nodes
        names(res) <- temp.names
    }

    if(useTipNames) {
        names(res) <- paste(names(t1), names(t2), sep="-")
    } else {
        names(res) <- paste(t1,t2,sep="-")
    }

    return(res)
} # end sp.tips



# examples
# source("/home/master/dev/adephylo/pkg/R/utils.R")
#phy <- as(rtree(15),"phylo4")
plot(phy,show.n=T)
tip1 <- "t1"
tip2 <- "t2"


sp.tips(phy, "t1", "t2")
sp.tips(phy, rep(1,15), 1:15)
sp.tips(phy, rep(1, 15), 1:15, TRUE)

## heavier tree
# x <- as(rtree(1000), "phylo4")
# system.time(sp.tips(x,1,1:1000))