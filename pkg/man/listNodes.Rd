\name{listNodes}
\alias{listNodes}
\title{List direct descendants for all nodes of a tree}
\description{
  The function \code{listNodes} lists the direct descendants from each node
  of a tree. The tree can be of class \code{\link[pkg:ape]{phylo}},
  \linkS4class{phylo4} or \linkS4class{phylo4d}.
}
\usage{
listNodes(x)
}
\arguments{
  \item{x}{A tree of  class \code{\link[pkg:ape]{phylo}},
    \linkS4class{phylo4} or \linkS4class{phylo4d}.}
}
\value{
  A list whose components are vectors of named nodes (or tips) for a
  given internal node.
}
\author{ Thibaut Jombart \email{jombart@biomserv.univ-lyon1.fr} }
\seealso{\code{\link{listTips}} which lists the tips descending from
  each node. \cr
  
  \code{\link{treePart}} which defines partitions of tips
  according to the tree topology.}
\examples{
if(require(ape) & require(phylobase)){
## make a tree
x <- as(rtree(20),"phylo4")
plot(x,show.node=TRUE)
listNodes(x)
}
}
\keyword{manip}
