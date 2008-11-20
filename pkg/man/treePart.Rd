\name{treePart}
\alias{treePart}
\title{Define partitions of tips according from a tree}
\description{
  The function \code{treePart} defines partitions of tips reflecting the
  topology of a tree. Each node but the root is translated into a dummy
  vector having one value for each tip: this value is '1' if the tip
  descends from this node, and '0' otherwise.\cr
}
\usage{
treePart(x)
}
\arguments{
  \item{x}{A tree of  class \code{\link[pkg:ape]{phylo}},
    \linkS4class{phylo4} or \linkS4class{phylo4d}.}
}
\value{
  A data.frame of dummy vectors (in columns) having one value for each
  tip (rows).
}
\author{ Thibaut Jombart \email{jombart@biomserv.univ-lyon1.fr} }
\seealso{\code{\link{listNodes}} which is called by \code{treePart}.
}
\examples{
if(require(ape) & require(phylobase)){
## make a tree
x <- as(rtree(10),"phylo4")
partition <- treePart(x)
partition

## plot the dummy vectors with the tree
temp <- phylo4d(x, partition)
s.phylo4d(temp, cent=FALSE, scale=FALSE)
}
}
\keyword{manip}
