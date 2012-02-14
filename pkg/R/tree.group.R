##############
## treeGroup
##############
treeGroup <- function(x, grp, dat=NULL, FUN=NULL, boot=FALSE, n.boot=10,
                      n.dens=4096, plot=TRUE, warn.lab=FALSE, ...){
    if(!require(ape)) stop("ape package is required")
    if(!inherits(x,"phylo")) stop("x is not a phylo object")
    if(boot && (is.null(dat) || is.null(FUN))) stop("dat and FUN must be provided for the bootstrap procedure")
    if(warn.lab && !is.null(dat) && !identical(x$tip.label,rownames(dat))) warning("Tip labels in x and rownames of dat differ \nplease make sure the same order is used in x, grp, and dat")
    grp <- factor(grp)
    K <- length(LEV <- levels(grp))
    N <- length(x$tip.label)
    D <- cophenetic.phylo(x)
    THRES <- 1e-320 # densities < THRES will be set to THRES to avoid log(x)=-Inf


    ## RE-ORDER GRP AND DATA MATRIX AFTER TIP LABELS ##
    if(!is.null(dat)){
        if(is.null(rownames(dat))) rownames(dat) <- x$tip.label
        if(!all(x$tip.label %in% rownames(dat))) stop("some tips do not have data matching their label")
        grp <- grp[match(x$tip.label, rownames(dat))] # grp is assumed to be in the same order as 'dat'
        dat <- dat[x$tip.label,,drop=FALSE]
    }

    #### AUXILIARY FUNCTIONS ####
    ## FUNCTION TO ESTIMATE A DENSITY AT A SERIES OF POINTS ##
    compute.dens <- function(dens, values){
        pred.y <- double(n <- length(values))
        return(.C("predict_density", dens$x, dens$y, length(dens$x), as.double(values), pred.y, n, PACKAGE="adephylo")[[5]])
    }


    ## FUNCTION TO GET A VECTOR OF PAIRWISE DISTANCE WITHIN ONE GROUP ##
    getdist.within.grp <- function(M, fac, val){ # val is one level of fac
        temp <- M[fac==val,fac==val]
        return(temp[lower.tri(temp)])
    }


    ## FUNCTION TO GET A VECTOR OF PAIRWISE DISTANCES WITHIN GROUP, FOR ALL GROUPS ##
    getdist.within.allgrp <- function(M, fac){
        res <- lapply(LEV, function(e) getdist.within.grp(M, fac, e))
        names(res) <- LEV
        return(res)
    }


    ## FUNCTION TO GET PROBA FOR ONE INDIV / ONE GROUP ##
    getprob.ind <- function(i, g){ # i: idx of indiv; g: idx of a group
        temp <- 1:ncol(D)
        dens <- compute.dens(list.dens[[g]], D[i,grp==LEV[g] & temp!=i])
        dens[dens < THRES] <- THRES
        res <- exp(mean(log(dens)))
        return(res)
    }


    ## FUNCTION TO GET PROBA FOR ALL INDIV / ONE GROUP ##
    getprob.grp <- function(g){ # g: idx of a group
        return(sapply(1:N, function(i) getprob.ind(i,g)))
    }


    ## FUNCTION TO GET A BOOTSTRAPPED TREE AND MATCHING GRP ##
    getboot.tree.grp <- function(){
        samp <- sample(1:N,replace=TRUE)
        tre <- FUN(dat[samp,,drop=FALSE])
        newgrp <- factor(grp[samp])
        return(list(tre=tre, grp=newgrp))
    }


    if(boot){ # USE BOOTSTRAP
        ## GET BOOTSTRAPPED TREES ##
        list.trees.grp <- lapply(1:n.boot, function(i) getboot.tree.grp())


        ## GET WITHIN-GROUP DISTANCES FOR EACH BOOTSTRAP SAMPLE ##
        list.D <- lapply(list.trees.grp, function(e) cophenetic.phylo(e$tre))
        temp <- lapply(1:n.boot, function(i) getdist.within.allgrp(list.D[[i]], list.trees.grp[[i]]$grp)) # for each replicate, list of distances within for each grp


        ## GET DENSITIES FOR EACH GROUP ##
        dens.dat <- lapply(LEV, function(onegroup) unlist(lapply(temp, function(e) e[[onegroup]]))) # density data for each group
        list.dens <- lapply(dens.dat, density, n=n.dens, ...) # densities for each group

    } else { # DO NOT USE BOOTSTRAP
         dens.dat <- getdist.within.allgrp(D, grp) # density data for each group
         list.dens <- lapply(dens.dat, density, n=n.dens, ...) # densities for each group
    }


    ## PLOT DENSITIES ##
    if(plot){
        find.mfrow <- function(i) {
            nrow <- ceiling(sqrt(i))
            ncol <- ceiling(i/ceiling(sqrt(i)))
            return(c(nrow,ncol))
        }
        par(mfrow = find.mfrow(K))
        for(i in 1:K){
            plot(list.dens[[i]], main=paste("Group:",LEV[i]),xlab="Within-group pairwise distance",ylab="Density", col="blue")
            points(dens.dat[[i]], rep(0,length(dens.dat[[i]])), pch="|", col="blue")
        }
    }


    ## COMPUTE MEMBERSHIP PROBABILITIES ##
    prob <- matrix(unlist(lapply(1:K, getprob.grp)), ncol=K)
    prob <- prop.table(prob,1)
    colnames(prob) <- LEV
    rownames(prob) <- x$tip.label


    ## FIND GROUP ASSIGNMENTS ##
    temp <- factor(colnames(prob)[apply(prob,1, which.max)])
    annot <- rep(" ", N)
    annot[as.character(grp)!=as.character(temp)] <- "!"
    groups <- data.frame(observed=grp, inferred=temp, annot=annot)
    rownames(groups) <- rownames(prob)


    ## BUILD / RETURN RESULT ##
    propcorrect <- mean(annot==" ")
    ## propcorrect.bygroup <- tapply(annot==" ", grp, mean)
    assignability <- mean((apply(prob,1,max)-.5)/.5)
    ##res <- list(prob=prob,groups=groups, mean.correct=propcorrect, prop.correct=propcorrect.bygroup)
    res <- list(prob=prob, groups=groups, assigndex=assignability, mean.correct=propcorrect)

    return(res)
} # end treeGroup








##############################
## simulate data with groups
##############################
simDatGroups <- function(k=2, p=1000, n=10, mu=0, sigma=1){
    ## RECYCLE ARGUMENTS ##
    n <- rep(n, length=k)
    mu <- rep(mu, length=k)
    sigma <- rep(sigma, length=k)


    ## GENERATE DATA ##
    dat <- list()
    for(i in 1:k){
        dat[[i]] <- replicate(p, rnorm(n[i], mu[i], sigma[i]))
    }

    dat <- Reduce(rbind,dat)
    rownames(dat) <- paste("ind", 1:nrow(dat))

    ## SHAPE OUTPUT ##
    grp <- factor(paste("grp", rep(1:k, n)))
    res <- list(dat=dat, grp=grp)
    return(res)
}





