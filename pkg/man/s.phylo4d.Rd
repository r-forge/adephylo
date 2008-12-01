\name{s.phylo4d}
\alias{s.phylo4d}
\title{scatterplot for a phylogeny and traits}
\description{
  This function displays scatterplots of traits provided for the tips of
  a phylogeny. Plotted objects must be valid \linkS4class{phylo4d}
  objects (implemented by the \code{phylobase} package). Current version
  allows plotting of a tree and one or more quantitative traits
  (possibly containing missing data, represented by an 'x').\cr
  
  The plot of phylogenies is performed by a call to \code{\link[pkg:ape]{plot.phylo}}
  from the \code{ape} package. Hence, many of the arguments of
  \code{\link[pkg:ape]{plot.phylo}} can be passed to \code{s.phylo4d},
  through the \dots argument, but their names must be complete.
}
\usage{
s.phylo4d(x, treetype=c("phylogram","cladogram"), symbol=c("circles", "squares"),
          center=TRUE, scale=TRUE, legend=TRUE, grid=TRUE, box=TRUE,
          show.tip.label=TRUE, show.node.label=TRUE, show.var.label=TRUE,
          ratio.tree=1/3, font=3,
         tip.label=x@tip.label, var.label=colnames(x@tip.data),
         cex.symbol=1, cex.label=1, cex.legend=1, coord.legend=NULL, \dots)
}
\arguments{
  \item{x}{a \linkS4class{phylo4d} object}
  \item{treetype}{the type of tree to be plotted ("phylogram" or "cladogram")}
  \item{symbol}{the type of symbol used to represent data ("circles" or "squares")}
  \item{center}{a logical stating whether variables should be centred
    (TRUE, default) or not (FALSE)}
  \item{scale}{a logical stating whether variables should be scaled
    (TRUE, default) or not (FALSE)}
  \item{legend}{a logical stating whether a legend should be added to
    the plot (TRUE) or not (FALSE, default)}
  \item{grid}{a logical stating whether a grid should be added to the plot
    (TRUE, default) or not (FALSE)}
  \item{box}{a logical stating whether a box should be added around the plot
    (TRUE, default) or not (FALSE)}
 \item{show.tip.label}{a logical stating whether tip labels should be
   printed (TRUE, default) or not (FALSE)}
 \item{show.node.label}{a logical stating whether node labels should be
   printed (TRUE, default) or not (FALSE)}
 \item{show.var.label}{a logical stating whether labels of variables should be
   printed (TRUE, default) or not (FALSE)}
 \item{ratio.tree}{the proportion of width of the figure occupied by the tree}
 \item{font}{an integer specifying the type of font for the labels: 1
   (plain text), 2 (bold), 3 (italic, default), or 4 (bold italic).}
 \item{tip.label}{a character vector giving the tip labels}
 \item{var.label}{a character vector giving the labels of variables}
 \item{cex.symbol}{a numeric giving the factor scaling the symbols}
 \item{cex.label}{a numeric giving the factor scaling all labels}
 \item{cex.legend}{a numeric giving the factor scaling the legend}
 \item{coord.legend}{an optional list with two components 'x' and 'y'
   indicating the lower-left position of the legend. Can be set to
   \code{locator(1) to position the legend interactively.}}
 \item{\dots}{further arguments to be passed to plot methods from
  \code{ape}. See \code{?\link[pkg:ape]{plot.phylo}}.}
}
\seealso{The \linkS4class{phylo4d} class for storing \code{phylogeny+data}.\cr

  \code{\link[pkg:ape]{plot.phylo}} from the \code{ape} package.\cr

  An alternative (deprecated) representation is available from
  \code{\link[pkg:ade4]{dotchart.phylog}}.
}
\details{
  The function \code{s.phylo4d} is based on former plot method for
  \linkS4class{phylo4d} objects from the \code{phylobase} package.
  It replaces the deprecated \code{ade4} functions
  \code{\link[pkg:ade4]{symbols.phylog}} and
  \code{\link[pkg:ade4]{table.phylog}}.
}
\author{Thibaut Jombart \email{jombart@biomserv.univ-lyon1.fr}}
\examples{
if(require(ape) & require(phylobase) & require(ade4)){

## simulated data
tr <- rtree(20)
dat <- data.frame(a = rnorm(20), b = scale(1:20), c=runif(20,-2,2) )
dat[3:6, 2] <- NA # introduce some NAs
obj <- phylo4d(tr, dat, use.tip=FALSE) # build a phylo4d object
s.phylo4d(obj) # default scatterplot
s.phylo4d(obj,cex.leg=.6, use.edge.length=FALSE) # customized
s.phylo4d(obj,treetype="clad", show.node=FALSE, cex.leg=.6,
use.edge.length=FALSE, edge.color="blue", edge.width=3) # more customized


## teleost fishes data
data(mjrochet)
temp <- read.tree(text=mjrochet$tre) # make a tree
mjr <- phylo4d(x=temp,tip.data=mjrochet$tab) # male a phylo4d object
s.phylo4d(mjr,cex.lab=.5,show.node=FALSE,symb="square")


## lizards data
data(lizards)
liz.tr <- read.tree(tex=lizards$hprA) # make a tree
liz <- phylo4d(liz.tr, lizards$traits) # make a phylo4d object
s.phylo4d(liz)


## plotting principal components
liz.pca1 <- dudi.pca(lizards$traits, scannf=FALSE, nf=2) # PCA of traits
myPC <- phylo4d(liz.tr, liz.pca1$li) # store PC in a phylo4d object
varlab <- paste("Principal \ncomponent", 1:2) # make labels for PCs
s.phylo4d(myPC, ratio=.8, var.lab=varlab) # plot the PCs
add.scatter.eig(liz.pca1$eig,2,1,2,posi="topleft", inset=c(0,.15))
title("Phylogeny and the principal components")
}
}
\keyword{hplot}
\keyword{multivariate}