\name{IDD.IRR}
\alias{IDD.IRR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Integration-driven Discovery and Integration-driven Revision Rates }
\description{
  Calculates the gain or the loss of differentially expressed genes due to meta-analysis
compared to individual studies.
}
\usage{
IDD.IRR(meta_de, ind_de)
}
\arguments{
  \item{meta_de}{ Vector of differentially expressed tags (or indices of these tags) with the meta-analysis}
  \item{ind_de}{ List of vectors storing differentially expressed tags (or indices of these tags) in each individual study }
}
\value{
  \item{DE }{Number of Differentially Expressed (DE) genes }
  \item{IDD }{Integration Driven Discoveries: number of genes that are declared DE in
the meta-analysis that were not identified in any of
the individual studies alone.}
  \item{Loss}{Number of genes that are declared DE in individual studies but not in meta-analysis.}
  \item{IDR}{Integration-driven Discovery Rate: proportion of genes that are identified as DE in the meta-analysis that were not
identified in any of the individual studies alone.}
  \item{IRR}{Integration-driven Revision Rate: percentage of genes
that are declared DE in individual studies but not in meta-analysis.}
}
\references{  Marot, G., Foulley, J.-L., Mayer, C.-D., Jaffrezic, F. (2009) Moderated effect size and p-value combinations for microarray meta-analyses. Bioinformatics. 25 (20): 2692-2699. }
\author{ Guillemette Marot }
\examples{
data(rawpval)
adjpval<-lapply(rawpval, FUN=function(x) p.adjust(x, method="BH"))
ind_smalladjp<-lapply(adjpval, FUN=function(x) which(x <= 0.05))
#indicators corresponding to the inverse normal p-value combination
invnormcomb <- invnorm(rawpval,nrep=c(8,8), BHth = 0.05) 
IDD.IRR(invnormcomb$DEindices,ind_smalladjp)
#indicators corresponding to the p-value combination with Fisher's method
fishcomb <- fishercomb(rawpval, BHth = 0.05)
IDD.IRR(fishcomb$DEindices,ind_smalladjp)
}
\keyword{ methods }
\keyword{ models }% __ONLY ONE__ keyword per line
