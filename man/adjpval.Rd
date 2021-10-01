\name{adjpval}
\alias{adjpval}
\docType{data}
\title{
Simulated adjusted p-values
}
\description{
The adjusted p-values provided here result from the following procedure:
 1) simulation of two RNA-seq experiments with four replicates in each condition via the \code{\link{sim.function}}, 2) analysis of differentially expressed tags using the DESeq package.
}
\usage{data(adjpval)}
\format{
List of length 2, where each list is a vector containing the adjusted p-values for 14,456 genes from
individual differential analyses (obtained using DESeq v1.8.3) of each of the simulated RNA-seq datasets. 
}
\details{
It is possible to reproduce these adjusted p-values using the procedure described in the package vignette.
}
\references{
A. Rau, G. Marot and F. Jaffrezic (2014). Differential meta-analysis of RNA-seq data. \emph{BMC Bioinformatics} \bold{15}:91
}
\examples{
data(adjpval)
## Maybe str(adjpval)
}
\keyword{datasets}
