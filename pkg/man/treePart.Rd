\name{treePart}
\alias{treePart}
\title{Define partitions of tips according from a tree}
\description{
  The function \code{treePart} defines partitions of tips reflecting the
  topology of a tree. There are two possible outputs (handled by the
  argument \code{result}):\cr
  - \code{basis} mode: each node but the root is translated into a dummy
  vector having one value for each tip: this value is '1' if the tip
  descends from this node, and '0' otherwise.\cr
  - \code{orthobasis}: in this mode, an orthonormal basis is derived
  from the basis previously mentionned. This orthobasis was proposed in
  the orthogram (Ollier \emph{et al.} 2006).
}
\usage{
treePart(x, result=c("basis", "orthobasis"))
}
\arguments{
  \item{x}{a tree of  class \code{\link[pkg:ape]{phylo}},
    \linkS4class{phylo4} or \linkS4class{phylo4d}.}
  \item{result}{a character string specifying the type of result: either
  a basis of dummy vectors (\code{basis}), or an orthobasis derived from these dummy vectors (\code{orthobasis}).}
}
\value{
  A matrix of dummy vectors (in columns) having one value for each
  tip (rows).
}
\author{ Thibaut Jombart \email{jombart@biomserv.univ-lyon1.fr} }
\seealso{
  - \code{\link{listDD}} which is called by \code{treePart}.\cr
  - \code{\link{orthogram}}, which uses by default the orthobasis
  produced by \code{treePart}.\cr
}
\details{
  Orthobasis produced by this function are identical to those stored in
  the \$Bscores component of deprecated \link[pkg:ade4]{phylog} objects,
  from the ade4 package.
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

if(require(ade4)){
## former orthobasis in ade4
data(ungulates)
ung.Bscores <- newick2phylog(ungulates$tre)$Bscores

## version using treePart
ung <- read.tree(text=ungulates$tre)
ung.orthobas <- treePart(ung, res="orthobasis")

## comparison of the first 3 vectors
cbind(ung.Bscores[,1], ung.orthobas[,1],ung.Bscores[,2], ung.orthobas[,2], ung.Bscores[,3], ung.orthobas[,3])
}
}
\keyword{manip}
