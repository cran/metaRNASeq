\name{invnorm}
\alias{invnorm}
\title{
P-value combination using the inverse normal method
}
\description{
Combines one sided p-values using the inverse normal method.
}
\usage{
invnorm(indpval, nrep, BHth = 0.05)
}
\arguments{
  \item{indpval}{
List of vectors of one sided p-values to be combined.
}
  \item{nrep}{
Vector of numbers of replicates used in each study to calculate the previous one-sided p-values. 
}
  \item{BHth}{
Benjamini Hochberg threshold. By default, the False Discovery Rate is controlled at 5\%.
}
}
\details{
For each gene \emph{g}, let
\deqn{N_g = \sum_{s=1}^S \omega_s \Phi^{-1}(1-p_{gs}),}
where \eqn{p_{gs}} corresponds to the raw \emph{p}-value obtained for gene \emph{g} in a differential
analysis for study \emph{s} (assumed to be uniformly distributed under the null hypothesis), \eqn{\Phi} the
cumulative distribution function of the standard normal distribution, and \eqn{\omega_s} a set of weights.
We define the weights \eqn{\omega_s} as in Marot and Mayer (2009):
\deqn{\omega_s = \sqrt{\frac{\sum_c R_{cs}}{\sum_\ell \sum_c R_{c\ell}}},}
where \eqn{\sum_c R_{cs}} is the total number of biological replicates in study \emph{s}. This allows
studies with large numbers of biological replicates to be attributed a larger weight than smaller studies.

Under the null hypothesis, the test statistic \eqn{N_g} follows a N(0,1) distribution. A unilateral
test on the righthand tail of the distribution may then be performed, and classical procedures for the
correction of multiple testing, such as that of Benjamini and Hochberg (1995), may subsequently be applied to
the obtained \emph{p}-values to control the false discovery rate at a desired level \eqn{\alpha}.
}
\value{
\item{DEindices }{Indices of differentially expressed genes at the chosen Benjamini Hochberg threshold.}
\item{TestStatistic }{Vector with test statistics for differential expression in the meta-analysis.}
\item{rawpval }{Vector with raw p-values for differential expression in the meta-analysis.}
\item{adjpval }{Vector with adjusted p-values for differential expression in the meta-analysis.}
}
\references{ 
Y. Benjamini and Y. Hochberg (1995). Controlling the false discovery rate: a pratical and powerful approach
to multiple testing. \emph{JRSS B} (57): 289-300.

Hedges, L. and Olkin, I. (1985). Statistical Methods for Meta-Analysis. Academic Press. 

Marot, G. and Mayer, C.-D. (2009). Sequential analysis for microarray data based on sensitivity and meta-analysis.
\emph{SAGMB} 8(1): 1-33.

A. Rau, G. Marot and F. Jaffrezic (2012). Differential meta-analysis of RNA-seq data. (Manuscript in progress)}
\note{
This function resembles the function \code{directpvalcombi} in the \emph{metaMA} R package; there is, however, one 
important difference in the calculation of \emph{p}-values. In particular, for microarray data, it is typically 
advised to separate under- and over-expressed genes prior to the meta-analysis. In the case of RNA-seq data, 
differential analyses from individual studies typically make use of negative binomial models and exact tests, 
which lead to one-sided, rather than two-sided, p-values. As such, we suggest performing a meta-analysis over 
the full set of genes, followed by an a posteriori check, and if necessary filter, of genes with conflicting 
results (over vs. under expression) among studies.
}
\seealso{
\code{\link{metaRNASeq}}
}
\examples{
data(rawpval)
## 8 replicates simulated in each study
invnormcomb <- invnorm(rawpval,nrep=c(8,8), BHth = 0.05)       
DE <- ifelse(invnormcomb$adjpval<=0.05,1,0)
hist(invnormcomb$rawpval,nclass=100)

## A more detailed example is given in the vignette of the package:
## vignette("metaRNASeq")
}
\keyword{ methods }
\keyword{ models }
