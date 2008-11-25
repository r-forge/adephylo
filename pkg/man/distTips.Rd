\name{distNodes}
\alias{distNodes}
\title{Compute some phylogenetic distance between two nodes}
\description{
  The function \code{distNodes} computes a given distance between two
  nodes of a phylogeny. This distance corresponds to the shortest path
  between the two nodes. Several distances can be used, defaulting to the
  sum of branch lengths (see argument \code{method}).
}
\usage{
distNodes(x, node1, node2, method=c("brlength","nNodes","Abouheif","sumDD"))
}
\arguments{
  \item{x}{a tree of  class \code{\link[pkg:ape]{phylo}},
    \linkS4class{phylo4} or \linkS4class{phylo4d}.}
  \item{method}{a character string (full or abbreviated without
    ambiguity) specifying the method used to compute distances ;
    possible values are:\cr
    - \code{brlength}: branch length \cr
    - \code{nNodes}: number of nodes on the path between the nodes \cr
    - \code{Abouheif}: Abouheif's distance (see details) \cr
    - \code{sumDD}: sum of direct descendants of all nodes on the path
    (see details) \cr
   }
}
\value{
  A numeric vector containing one distance value for each tip.
}
\details{
  \code{Abouheif} distance refers to the phylogenetic distance
  underlying the test of Abouheif (see references). Let P be the set of
  all the nodes in the path going from \code{node1} to \code{node2}. Let
  DDP be the number of direct descendants from each node in P. Then, the
  so-called 'Abouheif' distance is the product of all terms in DDP.\cr

   \code{sumDD} refers to a phylogenetic distance quite similar to that
  of Abouheif. We consider the same sets P and DDP. But instead of
  computing the product of all terms in DDP, this distance computes the
  sum of all terms in DDP.
}
\author{ Thibaut Jombart \email{jombart@biomserv.univ-lyon1.fr} }
\seealso{\code{\link{distTips}} which computes several phylogenetic
  distances between tips.
}
\references{
  Pavoine, S.; Ollier, S.; Pontier, D. & Chessel, D. (2008) Testing for
  phylogenetic signal in life history variable: Abouheif's test
  revisited. \emph{Theoretical Population Biology}: \bold{73}, 79-91.
}
\examples{
if(require(ape) & require(phylobase)){
## make a tree
x <- as(rtree(10),"phylo4")
plot(x, show.node=TRUE)
axisPhylo()
## compute different distances
distNodes(x, "t1", "t2")
distNodes(x, "t1", "t2", "nNodes")
distNodes(x, "t1", "t2", "Abouheif")
distNodes(x, "t1", "t2", "sumDD")
}
}
\keyword{manip}
