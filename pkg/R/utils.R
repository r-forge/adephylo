###############
# shortestPath -> moved to phylobase
###############
## shortestPath <- function(x, node1, node2){
##     if(!require(phylobase)) stop("phylobase package is not installed")

##     ## conversion from phylo, phylo4 and phylo4d
##     x <- as(x, "phylo4")

##     ## come checks
##     if (is.character(checkval <- check_phylo4(x))) stop(checkval)
##     t1 <- getnodes(x, node1)
##     t2 <- getnodes(x, node2)
##     if(any(is.na(c(t1,t2)))) stop("wrong node specified")
##     if(t1==t2) return(NULL)

##     ## main computations
##     comAnc <- MRCA(x, t1, t2) # common ancestor
##     desComAnc <- descendants(x, comAnc, which="all")
##     ancT1 <- ancestors(x, t1, which="all")
##     path1 <- intersect(desComAnc, ancT1) # path: common anc -> t1

##     ancT2 <- ancestors(x, t2, which="all")
##     path2 <- intersect(desComAnc, ancT2) # path: common anc -> t2

##     res <- union(path1, path2) # union of the path
##     ## add the common ancestor if it differs from t1 or t2
##     if(!comAnc %in% c(t1,t2)){
##         res <- c(comAnc,res)
##     }

##     res <- getnodes(x, res)

##     return(res)
## } # end shortestPath

