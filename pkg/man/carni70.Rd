\encoding{latin1}
\name{carni70}
\alias{carni70}
\docType{data}
\title{Phylogeny and quantitative traits of carnivora}
\description{
This data set describes the phylogeny of 70 carnivora as reported by Diniz-Filho and Torres (2002). It also gives the geographic range size and body size corresponding to these 70 species.
}
\usage{data(carni70)}
\format{
\code{carni70} is a list containing the 2 following objects: 
\describe{
   \item{tre}{is a character string giving the phylogenetic tree in Newick format.
   Branch lengths are expressed as divergence times (millions of years)} 
   \item{tab}{is a data frame with 70 species and two traits: size (body size (kg)) ; range (geographic range size (km)).}    
}}
\source{
Diniz-Filho, J. A. F., and N. M. T�rres. (2002) Phylogenetic comparative methods and the 
geographic range size-body size relationship in new world terrestrial carnivora. \emph{Evolutionary Ecology}, \bold{16}, 351--367.
}
\details{
  This dataset replaces the former version in ade4.
}
\examples{
data(carni70)
tre <- read.tree(text=carni70$tre)
x <- phylo4d(tre, carni70$tab)
s.phylo4d(x)

par(mar=rep(.1,4))
s.phylo4d(x,cex.lab=.5, show.n=FALSE, ratio=.5)

## transform size in log
size <- log(carni70$tab)[,1]
names(size) <- row.names(carni70$tab)
tdata(x) <- cbind.data.frame(sizeLog=size, tdata(x))
s.phylo4d(x,cex.lab=.5, show.n=FALSE, ratio=.5)

}
\keyword{datasets}