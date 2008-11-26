\name{distRoot}
\alias{distRoot}
\title{Compute the distance of tips to the root}
\description{
  The function \code{distRoot} computes the distance of each tip to the
  root. Several distances can be used, defaulting to the sum of branch lengths.
}
\usage{
distRoot(x, method=c("brlength","nAncestors","Abouheif"))
}
\arguments{
  \item{x}{a tree of  class \code{\link[pkg:ape]{phylo}},
    \linkS4class{phylo4} or \linkS4class{phylo4d}.}
  \item{method}{a character string (full or abbreviated without
    ambiguity) specifying the method used to compute distances ;
    possible values are:\cr
    - \code{brlength}: branch length \cr
    - \code{nAncestors}: number of Ancestors \cr
    - \code{Abouheif}: Abouheif's distance \cr
  }
}
\value{
  A numeric vector containing one distance value for each tip.
}
\author{ Thibaut Jombart \email{jombart@biomserv.univ-lyon1.fr} }
\seealso{\code{\link{distTips}} which computes several phylogenetic
  distances between tips.
}
\examples{
if(require(ape) & require(phylobase)){
## make a tree
x <- as(rtree(10),"phylo4")
D <- distRoot(x)
D

## plot these distances along with the tree
temp <- phylo4d(x, D)
s.phylo4d(temp, cent=FALSE, scale=FALSE)
}
}
\keyword{manip}
