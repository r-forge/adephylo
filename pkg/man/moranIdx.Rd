\name{Moran's index}
\alias{moran.idx}
\title{Computes Moran's index for a variable}
\description{
  This simple function computes Moran's index of autocorrelation given a
  variable and a matrix of proximities among observations.
}
\usage{
moran.idx(x, prox)
}
\arguments{
  \item{x}{a numeric vector whose autocorrelation is computed.}
  \item{prox}{a matrix of proximities between observations, as computed
   by the \code{\link{proxTips}}. Off-diagonal terms must be positive or
   null, while diagonal terms must all equal zero.}
}
\value{
  The numeric value of Moran's index.
}
\author{ Thibaut Jombart \email{jombart@biomserv.univ-lyon1.fr} }
\seealso{\code{\link{proxTips}} which computes phylogenetic proximities
  between tips of a phylogeny.
}
\examples{

## example(maples)
W <- proxTips(tre, met="Abouheif")
moran.idx(dom, W)
moran.idx(bif, W)
moran.idx(rnorm(nTips(tre)), W)

## build a simple permutation test for 'bif'
sim <- replicate(499, moran.idx(rnorm(nTips(tre)), W)) # permutations
sim <- c(sim, moran.idx(bif, W))

pval <- mean(sim>=sim[1]) # right-tail p-value
pval

hist(sim, col="grey", main="Moran's I Monte Carlo test for 'bif'") # plot
mtext("Histogram of permutations and observation (in red)")
abline(v=sim[1], col="red", lwd=3)

}
\keyword{manip}
