\name{Moran's eigenvectors}
\alias{orthobasis.phylo}
\alias{me.phylo}
\title{Computes Moran's eigenvectors from a tree or a phylogenetic
  proximity matrix}
\description{
  The function \code{orthobasis.phylo} (also nicknamed \code{me.phylo})
  computes Moran's eigenvectors (ME) associated to a tree. If the tree has
  'n' tips, (n-1) vectors will be produced. These vectors form an
  orthonormal basis: they are centred to mean zero, have unit variance,
  and are uncorrelated. Each vector models a different pattern of
  phylogenetic autocorrelation. The first vectors are those with maximum
  positive autocorrelation, while the last vectors are those with
  maximum negative autocorrelation. ME can be used, for instance, as
  regressors to remove phylogenetic autocorrelation from data (see
  references).\cr

  ME can be obtained from a tree, specifying the
  phylogenetic proximity to be used. Alternatively, they can be obtained
  directly from a matrix of phylogenetic proximities as constructed by
  \code{\link{proxTips}}.
}
\usage{
me.phylo(x=NULL, prox=NULL, method=c("patristic","nNodes","oriAbouheif","Abouheif","sumDD"), a=1)
orthobasis.phylo(x=NULL, prox=NULL, method=c("patristic","nNodes","oriAbouheif","Abouheif","sumDD"), a=1)
}
\arguments{
  \item{x}{A tree of  class \code{\link[pkg:ape]{phylo}},
    \linkS4class{phylo4} or \linkS4class{phylo4d}.}
  \item{prox}{a matrix of phylogenetic proximities as returned by \code{\link{proxTips}}.}
  \item{method}{a character string (full or abbreviated without
    ambiguity) specifying the method used to compute proximities;
    possible values are:\cr
    - \code{patristic}: (inversed sum of) branch lengths \cr
    - \code{nNodes}: (inversed) number of nodes on the path between the
    nodes \cr
    - \code{oriAbouheif}: original Abouheif's proximity, with diagonal
  (see details in \code{\link{proxTips}}) \cr
    - \code{Abouheif}: Abouheif's proximity (see details in \code{\link{proxTips}}) \cr
    - \code{sumDD}: (inversed) sum of direct descendants of all nodes on the path
    (see details in \code{\link{proxTips}}).}
  \item{a}{the exponent used to compute the proximity (see ?\code{\link{proxTips}}).}
}
\value{
  An object of class \code{orthobasis}. This is a data.frame with
  Moran's eigenvectors in column, with special attributes:\cr
  - attr(...,"values"): Moran's index for each vector
  - attr(...,"weights"): weights of tips; current implementation uses
  only uniform weights
}
\author{ Thibaut Jombart \email{jombart@biomserv.univ-lyon1.fr} }
\seealso{
  - \code{\link{proxTips}} which computes phylogenetic proximities
  between tips.\cr

  - \code{\link{treePart}} which can compute an orthobasis based on
    the topology of a phylogeny.\cr

}
\examples{

## make a tree
x <- rtree(50)

## compute Moran's eigenvectors
ME <- me.phylo(x, met="Abouheif")
ME

## plot the 10 first vectors
obj <- phylo4d(x, as.data.frame(ME[,1:10]))
table.phylo4d(obj, cex.sym=.7, cex.lab=.7)


## removing phylogenetic autocorrelation in a model ##
example(ungulates, package="adephylo")

lm1
resid1 <- residuals(lm1)
orthogram(resid1, tre) # residuals are autocorrelated

## compute Moran's eigenvectors (ME)
myME <- me.phylo(tre, method="Abou")
lm2 <- lm(neonatw ~ myME[,1] + afbw) # use for ME as covariable
resid2 <- residuals(lm2)
orthogram(resid2, tre) # there is no longer phylogenetic autocorrelation 

## see the difference
table.phylo4d(phylo4d(tre, cbind.data.frame(resid1, resid2)))

}
\keyword{manip}
