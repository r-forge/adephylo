\encoding{UTF-8}
\name{dibas}
\alias{dibas}
\alias{dibas.matrix}
\alias{dibas.dist}
\alias{dibas.phylo}
\alias{simDatGroups}
\title{DIstance-Based Assignment}
\description{
  These functions are under development. Please do not use them unless
  asked by the author.
}
\usage{
\method{dibas}{phylo}(x, grp, method = c("default", "leaveOneOut"),
      metric = c("nNodes", "patristic", "Abouheif", "sumDD"), \ldots)

\method{dibas}{matrix}(x, grp, method = c("default", "leaveOneOut"), ldots)

\method{dibas}{dist}(x, grp, method = c("default", "leaveOneOut"), \ldots)

}
\arguments{
  \item{x}{a \code{phylo} object, or a symmetric matrix of pairwise distances
    of class \code{matrix} or \code{dist}.}
  \item{grp}{a \code{factor} indicating the groups of observations.}
  \item{method}{a character string indicating the method to be used for
    estimating the distribution of pairwise distances within
    groups. The default method ("default") uses all observations, while
    the "leaveOneOut" estimates separate group distributions for each
    individual, leaving this one out in the estimation process.}
  \item{metric}{a character string matching "nNodes", "patristic", "Abouheif", or
    "sumDD" indicating the distance measure to be used. See
    \code{\link{distTips}} for more information. Note that patristic
    distances should be avoided in presence of one or more highly
    diverse group because of the 'hand fan' syndrome (see example).}
  \item{\dots}{further arguments passed to other methods. Can be used to
    provide arguments to \code{\link{table.phylo4d}} in \code{plot} method.}
}
\author{ Thibaut Jombart \email{tjombart@imperial.ac.uk} }
\examples{
if(require(ape)){
#### SIMPLE SIMULATED DATA ####
## 50 variables, 2 groups, 30 indiv
dat <- simDatGroups(k=2, p=50, n=c(15,15), mu=c(0,1))
names(dat)

## make a tree
tre <- nj(dist(dat$dat))
plot(tre,type="unr", tip.col=c("blue","red")[as.integer(dat$grp)],main="simulated data - tree")
 
## use dibas method
res <- dibas(tre, dat$grp, metric="nNodes")
res

barplot(t(res$prob), main="group membership probabilities")



#### NON-PARAMETRIC TEST BASED ON MEAN SUCCESSFUL ASSIGNMENT ####
## use dibas method
distHo <- replicate(100, dibas(tre, sample(dat$grp), metric="patristic")$mean.ok)
pval <- mean(res$mean.ok<=c(distHo,res$mean.ok))
pval

hist(c(distHo,res$mean.ok), col="grey", main="Mean successful assignement - permuted values")
abline(v=res$mean.ok, col="red")
mtext(side=3, text="Observed value in red")



#### HAND FAN SYNDROME ####
## 50 variables, 2 groups, 30 indiv
dat <- simDatGroups(k=2, p=50, n=c(15,15), mu=c(0,1), sigma=c(2,4))
names(dat)

## make a tree
tre <- nj(dist(dat$dat))
plot(tre,type="unr", tip.col=c("blue","red")[as.integer(dat$grp)],main="simulated data - tree")

## use dibas method
res <- dibas(tre, dat$grp, metric="nNodes")
res$grp.tab # results OK

res.patri <- dibas(tre, dat$grp, metric="patristic")
res.patri$grp.tab # poor results





#### MORE COMPLEX DATASET ####
if(require(adegenet)){

dat <- simDatGroups(k=5, p=50, n=c(5,10,10,30,60), mu=1:5, sigma=sample(1:5)/2)
names(dat)

## make a tree
tre <- nj(dist(dat$dat))
plot(tre,type="unr", tip.col=fac2col(dat$grp),main="simulated data - tree")

## use dibas method
res <- dibas(tre, dat$grp, metric="Abouheif")
res

}

}








}
\keyword{multivariate}