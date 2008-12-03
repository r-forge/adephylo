gearymoran <- function (x, W, nrepet=999,alter=c("greater", "less", "two-sided")) {

    ## some checks
    if(!require(ade4)) stop("The ade4 package is not installed.")
    alter <- match.arg(alter)

    ## checks for W
    if (any(W<0)) stop ("term <0 found in 'W'")
    if (nrow(W) != nobs) stop ("'W' is not squared")
    W <- as.matrix(W)
    nobs <- ncol(W)

    ## W has to be symmetric
    W <- (W + t(W))/2

    ## main computations
    test.names <- names(x)
    x <- data.matrix(x)
    if (nrow(x) != nobs) stop ("non convenient dimension")
    nvar <- ncol(x)
    res <- .C("gearymoran",
        param = as.integer(c(nobs,nvar,nrepet)),
        data = as.double(x),
        W = as.double(W),
        obs = double(nvar),
        result = double (nrepet*nvar),
        obstot = double(1),
        restot = double (nrepet),
        PACKAGE="adephylo"
    )
    res <- as.krandtest(obs=res$obs,sim=matrix(res$result,ncol=nvar, byr=TRUE),names=test.names,alter=alter)
    return(res)
} # end gearymoran
