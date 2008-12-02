\encoding{latin1}
\name{ungulates}
\alias{ungulates}
\docType{data}
\title{Phylogeny and quantitative traits of ungulates.}
\description{
This data set describes the phylogeny of 18 ungulates as reported by P\'elabon et al. (1995). It also gives 4 traits corresponding to these 18 species.
}
\usage{data(ungulates)}
\format{
\code{fission} is a list containing the 2 following objects : 
\describe{ 
   \item{tre}{is a character string giving the phylogenetic tree in Newick format.} 
   \item{tab}{is a data frame with 18 species and 4 traits}    
}}
\details{        
  Variables of \code{ungulates$tab} are the following ones : \cr
  
  - afbw: is a numeric vector that describes the adult female body weight (g) \cr
  - mnw: is a numeric vector that describes the male neonatal weight (g) \cr
  - fnw: is a numeric vector that describes the female neonatal weight (g) \cr
  - ls: is a numeric vector that describes the litter size \cr
}
\source{
Data were obtained from P\'elabon, C., Gaillard, J.M., Loison, A. and Portier, A. (1995) 
Is sex-biased maternal care limited by total maternal expenditure in polygynous ungulates? 
\emph{Behavioral Ecology and Sociobiology}, \bold{37}, 311--319.
}
\note{
  This dataset replaces the former version in ade4.
}
\examples{
## load data
data(ungulates)
tre <- read.tree(text=ungulates$tre)
plot(tre)

## look at two traits
afbw <- log(ungulates$tab[,1])
neonatw <- log((ungulates$tab[,2]+ungulates$tab[,3])/2)
names(afbw) <- tre$tip.label
names(neonatw) <- tre$tip.label
plot(afbw, neonatw) # relationship between traits
abline(lm(neonatw~afbw))
x <- phylo4d(tre, cbind.data.frame(afbw, neonatw)) # traits on the phylogeny


## orthogram(ung.x,ung.phy$Bscores,nrep=9999)
## ung.z <- residuals(lm(ung.y~ung.x))
## names(ung.z) <- names(ung.phy$leaves)
## dotchart.phylog(ung.phy,ung.z,cle=1.5,cno=1.5,cdot=1,ceti=0.75)
## orthogram(ung.z,ung.phy$Bscores,nrep=9999)
}
\keyword{datasets}
