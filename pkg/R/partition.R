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





###########
# treePart
###########
treePart <- function(x, res=c("basis", "orthobasis")){
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

    if(res=="basis"){
        return(res)
    }



    ## If orthobasis is required ##

    ## Find values 'w' for all nodes
    ##
    ## Notations:
    ## - n: an internal node (HTU)
    ## - Dn: the set of all internal nodes descending from 'n'
    ## - En: the set 'n U Dn' (that is, Dn plus n itself)
    ## - ndd(e): the number of direct descendants from a node 'e'
    ##
    ## Then the values 'w' are computed as:
    ##
    ## w(n) = sum_{e \in En} lgamma( ndd(e) + 1)
    ##

    nbOfDD <- sapply(listDD(x), length) # nb of DD for each node
    HTU.idx <- (n+1):(n+nNodes(x)) # index of internal nodes (HTU)
    names(nbOfDD) <- HTU.idx # used to match the results of Dn

    findAlldHTU <- function(node){ # find all HTU descending from a node
        res <- descendants(x, node, which="all") # tips and HTU
        res <- res[res > n] # only HTU (here, just node numbers are kept
        if(length(res)==0) return(NULL)
        return(res)
    }


    listAlldHTU <- lapply(HTU.idx, function(node) c(node,findAlldHTU(node))) # ='Dn': for each HTU, list all HTU descending from it

    w <- sapply(listAlldHTU, function(e) sum(lgamma(nbOfDD[as.character(e)]+1))) # w(n)
    ## w stores the w(n) values.



    ## sorting of dummy vectors according to val

    ## discard dummy vectors for each node with smallest value

} # end treePart
