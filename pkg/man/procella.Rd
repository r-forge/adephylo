\name{procella}
\alias{procella}
\docType{data}
\title{Phylogeny and quantitative traits of birds}
\description{
This data set describes the phylogeny of 19 birds as reported by Bried et al. (2002). It also gives 6 traits corresponding to these 19 species.
}
\usage{data(procella)}
\format{
\code{procella} is a list containing the 2 following objects:  
\describe{     
   \item{tre}{is a character string giving the phylogenetic tree in Newick format.} 
   \item{traits}{is a data frame with 19 species and 6 traits}    
}}
\details{
Variables of \code{procella$traits} are the following ones: \cr
- site.fid: a numeric vector that describes the percentage of site fidelity\cr
- mate.fid: a numeric vector that describes the percentage of mate fidelity\cr
- mass: an integer vector that describes the adult body weight (g)\cr
- ALE: a numeric vector that describes the adult life expectancy (years)\cr
- BF: a numeric vector that describes the breeding frequencies\cr
- col.size: an integer vector that describes the colony size (no nests monitored)
}
\references{
Bried, J., Pontier, D. and Jouventin, P. (2002) Mate fidelity in monogamus birds: a re-examination of the Procellariiformes. 
\emph{Animal Behaviour}, \bold{65}, 235--246. 

See a data description at \url{http://pbil.univ-lyon1.fr/R/pps/pps037.pdf} (in French).
}
\details{
  This dataset replaces the former version in ade4.
}
\examples{
## load data, make tree and phylo4d object
data(procella)
tre <- read.tree(text=procella$tre)
x <- phylo4d(tre, procella$traits)
par(mar=rep(.1,4))
s.phylo4d(x,cex.lab=.7)

## test phylogenetic autocorrelation in traits
if(require(ade4)){

f1 <- function(vec){ # to replace missing data
m <- mean(vec,na.rm=TRUE)
vec[is.na(vec)] <- m
return(vec)
}

traits <- data.frame(lapply(procella$traits, f1))
prox <- proxTips(tre)

myTests <- gearymoran(prox, traits)
plot(myTests)
}

}
\keyword{datasets}