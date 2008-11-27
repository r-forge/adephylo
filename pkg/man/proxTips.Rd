\name{proxTips}
\alias{proxTips}
\title{Compute some phylogenetic proximities between tips}
\description{
  The function \code{proxTips} computes a given proximity between a set
  of tips of a phylogeny. A vector of tips is supplied: proximities
  between all possible pairs of these tips are computed.
  The proximities are computed from the shortest path
  between the tips. \cr

  Proximities are computed as the inverse (to the power \code{a}) of a
  phylogenetic distance (computed by \code{\link{distTips}}. Denoting
  \eqn{D=[d_{ij}]} a matrix of phylogenetic distances, the proximity matrix \eqn{M=[m_{ij}]} is
  computed as:
  \deqn{m_{ij} = \frac{1}{d_{ij}^a} \forall i \neq j}{ m_{ij} = (1/d_{ij})^a for all i different from j}
  and
   \deqn{m_{ii} = 0}

   ~\cr
  Several distances can be used, defaulting to the 
  sum of branch lengths (see argument \code{method}).
  Proximities are not true similarity measures, since the proximity of a tip with
  itself is always set to zero.\cr

  The obtained matrix of phylogenetic proximities (M) defines a bilinear
  symmetric form when M is symmetric (default):\cr
  \deqn{f(x,y) = x^{T}My}

  ~\cr
  In general, M is not a metric because it is not positive-definite.
  Such a matrice can be used to measure phylogenetic autocorrelation
  (using Moran's index):
  \deqn{I(x) = \frac{x^TMx}{var(x)}}{I(x) = (x^{T}Mx)/(var(x))  }
  
  or to compute lag vectors (Mx) used in autoregressive
  models, like:
  \deqn{x = Mx + ... + e}
  where '...' is the non-autoregressive part of the model, and 'e' are residuals.
}
\usage{
proxTips(x, tips, method=c("brlength","nNodes","Abouheif","sumDD"),
         a=1, normalize=c("row","col","none"), symmetric=TRUE)
}
\arguments{
  \item{x}{a tree of  class \code{\link[pkg:ape]{phylo}},
    \linkS4class{phylo4} or \linkS4class{phylo4d}.}
  \item{tips}{A vector of integers identifying tips by their numbers, or
    a vector of characters identifying tips by their names. Distances
    will be computed between all possible pairs of tips.}
  \item{method}{a character string (full or abbreviated without
    ambiguity) specifying the method used to compute proximities;
    possible values are:\cr
    - \code{brlength}: (inversed sum of) branch length \cr
    - \code{nNodes}: (inversed) number of nodes on the path between the nodes \cr
    - \code{Abouheif}: Abouheif's proximity (see details) \cr
    - \code{sumDD}: (inversed) sum of direct descendants of all nodes on the path
    (see details) \cr
  }
  \item{a}{the exponent used to compute the proximity}
  \item{normalize}{a character string specifying whether the matrix must
  be normalized by row (\code{row}), column (\code{col}), or not
  (\code{none}). Normalization amounts to dividing each row (or column)
  so that the marginal sum is 1. Hence, default is matrix with each row
  summing to 1.}
\item{symmetric}{a logical stating whether M must be coerced to be
  symmetric (TRUE, default) or not. This is achieved by taking (denoting N the
  matrix of proximities before re-symmetrization):
  \deqn{M = 0.5 * (N + N^{T})}
  Note that \eqn{x^{T}Ny = x^{T}My}, but the latter has the advantage of
  using a bilinear symmetric form (more appropriate for optimization
  purposes).
}
}
\value{
  A matrix of phylogenetic proximities.
}
\details{
  \code{Abouheif} proximity refers to the phylogenetic proximity
  underlying the test of Abouheif (see references). Let P be the set of
  all the nodes in the path going from \code{node1} to \code{node2}. Let
  DDP be the number of direct descendants from each node in P. Then, the
  so-called 'Abouheif' distance is the inverse of the product of all terms in DDP.\cr

   \code{sumDD} refers to a phylogenetic proximity quite similar to that
  of Abouheif. We consider the same sets P and DDP. But instead of
  taking the inverse of the product of all terms in DDP, this proximity computes the
  inverse of the sum of all terms in DDP.
}
\author{ Thibaut Jombart \email{jombart@biomserv.univ-lyon1.fr} }
\seealso{\code{\link{proxTips}} which computes several phylogenetic
  distances between tips.
}
\references{
  - For Moran's index with various proximities:\cr
  Pavoine, S.; Ollier, S.; Pontier, D.; Chessel, D. (2008) Testing for
  phylogenetic signal in life history variable: Abouheif's test
  revisited. \emph{Theoretical Population Biology}: \bold{73}, 79-91.

  - For regression on phylogenetic lag vector:\cr
  Cheverud, J. M.; Dow, M. M.; Leutenegger, W. (1985) The quantitative
  assessment of phylogentic constaints in comparative analyses: sexual
  dimorphism in body weights among primates. \emph{Evolution} \bold{39},
  1335-1351.\cr

  Cheverud, J. M.; Dow, M. M. (1985) An autocorrelation analysis of genetic
  variation due to lineal fission in social groups of Rhesus macaques.
  \emph{American Journal of Phyisical Anthropology} \bold{67}, 113-121.
}
\examples{
if(require(ape) & require(phylobase)){
## make a tree
x <- as(rtree(10),"phylo4")
plot(x, show.node=TRUE)
axisPhylo()
## compute different distances
proxTips(x, 1:5)
proxTips(x, 1:5, "nNodes")
proxTips(x, 1:5, "Abouheif")
proxTips(x, 1:5, "sumDD")
}
}
\keyword{manip}
