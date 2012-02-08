treeGroup <- function(x, grp, dat, FUN, n.boot=10,
                      n.dens=4096,
                      plot=TRUE, ...){
    if(!require(ape)) stop("ape package is required")
    if(!inherits(x,"phylo")) stop("x is not a phylo object")
    grp <- factor(grp)
    k <- length(lev <- levels(grp))
    n <- nrow(dat)
    D <- cophenetic.phylo(x)


    ## FUNCTION TO ESTIMATE A DENSITY AT A GIVEN POINT ##
    find.dens <- function(dens, value){
        if(value<=min(dens$x)) return((dens$y[1]+0) / 2)
        if(value>=max(dens$x)) return((dens$y[length(dens$y)]+0) / 2)
        idx <- min(which(c(dens$x > value)))
        return(mean(dens$y[c(idx,idx+1)]))
    }


    ## FUNCTION TO GET A VECTOR OF PAIRWISE DISTANCE FOR ONE GROUP ##
    getdist.grp <- function(M, g){
        temp <- M[grp==g,grp==g]
        return(temp[lower.tri(temp)])
    }


    ## FUNCTION TO GET PROBA FOR ONE INDIV / ONE GROUP ##
    getprob.grp <- function(i, g){ # g indicates a group
        find.dens(list.dens[[g]] D[i,grp==lev[g]])
        temp <- D[,grp==g]

    }


    ## PERFORM BOOTSTRAP TREES ##
    list.trees <- lapply(1:n.boot, function(i) FUN(dat[sample(1:n,replace=TRUE),]))


    ## GET WITHIN-GROUP DISTANCES FOR EACH BOOTSTRAP SAMPLE ##
    list.D <- lapply(list.trees, cophenetic.phylo)
    list.D <- lapply(lev, function(g) unlist(lapply(listD, function(e) getdist.grp(e,g))))


    ## COMPUTE DENSITIES ##
    list.dens <- lapply(list.D, density, n=n.dens, ...)
    if(plot){
        par(mfrow = c(ceiling(sqrt(k)),ceiling(sqrt(k))) )
        for(i in 1:k){
            plot(list.dens[[i]], main=paste("Group:",lev[i]),xlab="Phylogenetic pairwise distance",ylab="Density", col="blue")
            points(list.D[[i]], rep(0,length(list.D[[i]])), pch="|", col="blue")
        }
    }




}
