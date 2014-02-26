\name{dispFuncs}
\alias{dispFuncs}
\docType{data}
\title{
Gamma regression parameters describing the mean-dispersion relationship for two real datasets.
}
\description{
Gamma regression parameters describing the mean-dispersion relationship for each of the 
two real datasets considered in the associated paper, as estimated using the DESeq
package version 1.8.3 (Anders and Huber, 2010).
}
\usage{data(dispFuncs)}
\format{
List of length 2, where each list is a vector containing the two estimated coefficients (\eqn{\alpha_0} and 
\eqn{\alpha_1}) for the gamma regression in each study (see details below).
}
\details{
The \code{dispFuncs} object contains the estimated coefficients from the parametric gamma regressions
describing the mean-dispersion relationship for the two real datasets considered in the associated paper.
The gamma regressions were estimated using the DESeq package version 1.8.3 (Anders and Huber, 2010).

Briefly, after estimating a per-gene mean expression and dispersion values, the DESeq package fits a curve 
through these estimates. These fitted values correspond to an estimation of the typical relationship between 
mean expression values \eqn{\mu} and dispersions \eqn{\alpha} within a given dataset. By default, this relationship 
is estimated using a gamma-family generalized linear model (GLM), where two coefficients \eqn{\alpha_0} and \eqn{\alpha_1} 
are found to parameterize the fit as \eqn{\alpha = \alpha_0 + \alpha_1 / \mu}.

For the first dataset (F078), the estimated mean-dispersion relationship is described using the following
gamma-family GLM:
\deqn{\alpha = 0.024 + 14.896 / \mu.}

For the second dataset (F088), the estimated mean-dispersion relationship is described using the following
gamma-family GLM:
\deqn{\alpha = 0.00557 + 1.54247 / \mu.}

These gamma-family GLMs describing the mean-dispersions relationship in each of the two datasets are used in this
package to simulate data using dispersion parameters that are as realistic as possible.
}
\references{
A. Rau, G. Marot and F. Jaffrezic (submitted, 2012). Differential meta-analysis of RNA-seq data.

S. Anders and W. Huber (2010). Differential expression analysis for sequence count data. 
\emph{Genome Biology}, 11:R106.
}
\examples{
data(dispFuncs)
}
\keyword{datasets}
\seealso{
\code{\link{sim.function}}
}
