treeGroup <- function(x, grp, dat, FUN, n.boot=10,
                      n.dens=4096, plot=TRUE, check.lab=TRUE, ...){
    if(!require(ape)) stop("ape package is required")
    if(!inherits(x,"phylo")) stop("x is not a phylo object")
    if(check.lab && !identical(x$tip.label,rownames(dat))) warning("Tip labels in x and rownames of dat differ \nplease make sure the same order is used in x, grp, and dat")
    grp <- factor(grp)
    K <- length(lev <- levels(grp))
    N <- nrow(dat)
    D <- cophenetic.phylo(x)


    #### AUXILIARY FUNCTIONS ####
    ## FUNCTION TO ESTIMATE A DENSITY AT A SERIES OF POINTS ##
    compute.dens <- function(dens, values){
        pred.y <- double(n <- length(values))
        return(.C("predict_density", dens$x, dens$y, length(dens$x), as.double(values), pred.y, n, PACKAGE="adephylo")[[5]])
    }


    ## FUNCTION TO GET A VECTOR OF PAIRWISE DISTANCE FOR ONE GROUP ##
    getdist.grp <- function(M, g){
        temp <- M[grp==g,grp==g]
        return(temp[lower.tri(temp)])
    }


    ## FUNCTION TO GET PROBA FOR ONE INDIV / ONE GROUP ##
    getprob.ind <- function(i, g){ # i: idx of indiv; g: idx of a group
        temp <- 1:ncol(D)
        find.dens(list.dens[[g]], D[i,grp==lev[g] & temp!=i])
        res <- exp(sum(log(find.dens(list.dens[[g]], D[i,grp==lev[g] & temp!=i]))))
        return(res)
    }

    ## FUNCTION TO GET PROBA FOR ALL INDIV / ONE GROUP ##
    getprob.grp <- function(g){ # g: idx of a group
        return(sapply(1:N, function(i) getprob.ind(i,g)))
    }



    ## GET BOOTSTRAPPED TREES ##
    list.trees <- lapply(1:n.boot, function(i) FUN(dat[sample(1:N,replace=TRUE),]))


    ## GET WITHIN-GROUP DISTANCES FOR EACH BOOTSTRAP SAMPLE ##
    list.D <- lapply(list.trees, cophenetic.phylo)
    list.D <- lapply(lev, function(g) unlist(lapply(list.D, function(e) getdist.grp(e,g)))) # within dist. pooled across trees for each grp


    ## COMPUTE DENSITIES ##
    list.dens <- lapply(list.D, density, n=n.dens, ...)
    if(plot){
        find.mfrow <- function(i) {
            nrow <- ceiling(sqrt(i))
            ncol <- ceiling(i/ceiling(sqrt(i)))
            return(c(nrow,ncol))
        }
        par(mfrow = find.mfrow(K))
        for(i in 1:K){
            plot(list.dens[[i]], main=paste("Group:",lev[i]),xlab="Phylogenetic pairwise distance",ylab="Density", col="blue")
            points(list.D[[i]], rep(0,length(list.D[[i]])), pch="|", col="blue")
        }
    }


    ## COMPUTE MEMBERSHIP PROBABILITIES ##
    prob <- matrix(unlist(lapply(1:K, getprob.grp)), ncol=K)
    prob <- prop.table(prob,1)
    colnames(prob) <- lev
    rownames(prob) <- x$tip.label


    ## FIND GROUP ASSIGNMENTS ##
    temp <- factor(colnames(prob)[apply(prob,1, which.max)])
    annot <- rep(" ", N)
    annot[as.character(grp)!=as.character(temp)] <- "!"
    groups <- data.frame(observed=grp, inferred=temp, annot=annot)
    rownames(groups) <- x$tip.labels


    ## BUILD / RETURN RESULT ##
    res <- list(prob,groups)
    return(res)
} # end treeGroup
