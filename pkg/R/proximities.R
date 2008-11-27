###########
# proxTips
###########
proxTips <- function(x, tips="all",
                      method=c("brlength","nNodes","Abouheif","sumDD"),
                     a=1, normalize=c("row","col","none"), symmetric=TRUE){

    if(!require(phylobase)) stop("phylobase package is not installed")

    ## handle arguments
    x <- as(x, "phylo4")
    method <- match.arg(method)
    normalize <- match.arg(normalize)
    N <- nTips(x)
    if(tips[1]=="all") { tips <- 1:N }
    tips <- getnodes(x, tips)
    tips.names <- names(tips)

    ## some checks
    if (is.character(checkval <- check_phylo4(x))) stop(checkval)
    if(any(is.na(tips))) stop("wrong tips specified")

    ## compute distances
    D <- distTips(x, tips=tips, method=method)
    D <- as.matrix(D)

    ## compute proximities
    res <- (1/D)^a
    diag(res) <- 0

    ## standardization
    if(normalize=="row") {
        res <- prop.table(res, 1)
    }

    if(normalize=="col") {
        res <- prop.table(res, 2)
    }

    ## re-symmetrize
    if(symmetric){
        D <- 0.5 * (D + t(D))
    }

    ## set the output
    return(res)

} # end proxTips
