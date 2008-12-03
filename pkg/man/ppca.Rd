\encoding{UTF-8}
\name{ppca}
\alias{ppca}
\alias{print.ppca}
\alias{plot.ppca}
\title{Phylogenetic principal component analysis}
\description{These functions are designed to perform a phylogenetic principal
  component analysis (pPCA, Jombart et al. in prep) and to display the
  results. As this method is not published yet, please email the author
  before using it.

  \code{ppca} performs the phylogenetic component analysis. Other
  functions are:\cr

  - \code{print.ppca}: prints the ppca content\cr

  - \code{plot.ppca}: plot principal components using
  \code{\link{s.phylo4d}}\cr

}
\usage{
ppca(x, prox=NULL, method=c("patristic","nNodes","oriAbouheif","Abouheif","sumDD"), a=1,
                 center=TRUE, scale=TRUE, scannf=TRUE, nfposi=1, nfnega=0)

\method{print}{ppca}(x, \dots)

\method{plot}{ppca}(x, axis = 1:ncol(x$li), useLag=FALSE, \dots)
}
\arguments{
  \item{x}{a \linkS4class{phylo4d} object (for \code{ppca}) or a ppca
    object (for \code{print} and \code{plot} methods).}
  \item{prox}{a marix of phylogenetic proximities as returned by
    \code{\link{proxTips}}. If not provided, this matrix will be
    constructed using the arguments \code{method} and \code{a}.}
  \item{method}{a character string (full or abbreviated without
    ambiguity) specifying the method used to compute proximities;
    possible values are:\cr
    - \code{patristic}: (inversed sum of) branch lengths \cr
    - \code{nNodes}: (inversed) number of nodes on the path between the
    nodes \cr
    - \code{oriAbouheif}: original Abouheif's proximity, with diagonal (see details in \code{\link{proxTips}}) \cr
    - \code{Abouheif}: Abouheif's proximity (see details in \code{\link{proxTips}}) \cr
    - \code{sumDD}: (inversed) sum of direct descendants of all nodes on the path
    (see details in \code{\link{proxTips}}).}
  \item{a}{the exponent used to compute the proximity (see ?\code{\link{proxTips}}).}
  \item{center}{a logical indicating whether traits should be centred to
    mean zero (TRUE, default) or not (FALSE).}
  \item{scale}{a logical indicating whether traits should be scaled to
    unit variance (TRUE, default) or not (FALSE).}
  \item{scannf}{a logical stating whether eigenvalues should be chosen
    interactively (TRUE, default) or not (FALSE).}
  \item{nfposi}{an integer giving the number of positive eigenvalues retained
    ('global structures').}
  \item{nfnega}{an integer giving the number of negative eigenvalues retained
    ('local structures').}
  \item{\dots}{further arguments passed to other methods. Can be used to
  provide arguments to \code{\link{s.phylo4d}} in \code{plot} method.}
  \item{axis}{the index of the principal components to be represented.}
  \item{useLag}{a logical stating whether the lagged components
    (\code{x\$ls}) should be used instead of the components (\code{x\$li}).}
}
\details{
  The phylogenetic Principal Component Analysis (pPCA, Jombart et al.,
  in prep) is derived from the spatial Principal Component Analysis
  (spca, Jombart et al. 2008), implemented in the adegenet package (see
  \code{\link[pkg:adegenet]{spca}}).\cr

  pPCA is designed to investigate phylogenetic patterns in the variability of a set of
  traits. The analysis returns principal components maximizing the
  product of variance and phylogenetic autocorrelation (Moran's
  I), therefore reflecting biodemographic strategies that are linked
  to the phylogeny. Large positive and large negative eigenvalues
  correspond to global and local structures.\cr
}
\value{The class \code{ppca} are given to lists with the following
  components:\cr
  \item{eig}{a numeric vector of eigenvalues.}
  \item{nfposi}{an integer giving the number of global structures
    retained.}
  \item{nfnega}{an integer giving the number of local structures retained.}
  \item{c1}{a data.frame of loadings of traits for each axis.}
  \item{li}{a data.frame of coordinates of taxa onto the ppca axes
  (i.e., principal components).}
  \item{ls}{a data.frame of lagged prinpal components; useful to
  represent of global scores.}
  \item{as}{a data.frame giving the coordinates of the axes of an 'ordinary' PCA onto the
    ppca axes.}
  \item{call}{the matched call.}
  \item{tre}{a phylogenetic tre with class \linkS4class{phylo4}.}
  \item{prox}{a matrix of phylogenetic proximities.}
  
  Other functions have different outputs:\cr

  - \code{plot.ppca} returns the matched call.\cr

  - \code{screeplot.ppca} returns the matched call.  
}
\references{
  Jombart, T.; Pavoine, S.; Dufour, A.-B. & Pontier, D. (in prep)
  Exploring phylogeny as a source of ecological variation: a
  methodological approach \cr

  Jombart, T., Devillard, S., Dufour, A.-B. and Pontier, D.
  Revealing cryptic phylogenetic patterns in genetic variability by a new
  multivariate method. \emph{Heredity}, \bold{101}, 92--103.

  Wartenberg, D. E. (1985) Multivariate phylogenetic correlation: a method for
  exploratory geographical analysis. \emph{Geographical Analysis},
  \bold{17}, 263--283.

  Moran, P.A.P. (1948) The interpretation of statistical
  maps. \emph{Journal of the Royal Statistical Society, B}
  \bold{10}, 243--251.

  Moran, P.A.P. (1950) Notes on continuous stochastic
  phenomena. \emph{Biometrika}, \bold{37}, 17--23.

  de Jong, P. and Sprenger, C. and van Veen, F. (1984) On extreme values
  of Moran's I and Geary's c. \emph{Geographical Analysis}, \bold{16}, 17--24.

}
\seealso{The implementation of \code{\link[pkg:adegenet]{spca}} in the
  adegenet package (?\code{\link[pkg:adegenet]{adegenet}}) \cr
 }
\author{ Thibaut Jombart \email{jombart@biomserv.univ-lyon1.fr} }
\examples{
## build an look at data
data(maples)
tre <- read.tree(text=maples$tre)
x <- phylo4d(tre, maples$tab)
omar <- par("mar")
par(mar=rep(.1,4))
s.phylo4d(x, cex.lab=.5, cex.sym=.6, ratio=.1) # note NAs in last trait ('x')

## function to replace NA
f1 <- function(vec){
if(any(is.na(vec))){
m <- mean(vec, na.rm=TRUE)
vec[is.na(vec)] <- m
}
return(vec)
}

## compute a PCA
dat <- apply(maples$tab,2,f1) # replace NAs
x.noNA <- phylo4d(tre, as.data.frame(dat))
if(require(ade4)){
ppca1 <- ppca(x.noNA, scannf=FALSE, method="Abouheif")
ppca1

## some graphics
par(mar=omar)
barplot(ppca1$eig, main="Eigenvalues") # screeplot

a <- ppca1$c1[,1] # loadings on PC 1
names(a) <- row.names(ppca1$c1)
dotchart(a) # plot of loadings
abline(v=median(a), lty=3)
highContrib <- a[a< quantile(a,0.1) | a>quantile(a,0.9)]

plot(ppca1) # ppca plot

datSel <- cbind.data.frame(dat[, names(highContrib)], ppca1$li)
temp <- phylo4d(tre, datSel)
s.phylo4d(temp) # plot of most structured traits

## phylogenetic autocorrelation tests for these traits
prox <- proxTips(tre, method="Abouheif")
gearymoran(prox, dat[, names(highContrib)]) # one test per trait
}
}
\keyword{multivariate}