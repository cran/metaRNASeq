\name{fishercomb}
\alias{fishercomb}
\title{
P-value combination using Fisher's method
}
\description{
Combines one sided p-values using Fisher's method.
}
\usage{
fishercomb(indpval, BHth = 0.05)
}
\arguments{
  \item{indpval}{
List of vectors of one sided p-values to be combined.
}
  \item{BHth}{
Benjamini Hochberg threshold. By default, the False Discovery Rate is controlled at 5\%.
}
}
\details{
The test statistic for each gene \emph{g} is defined as
\deqn{F_g = -2 \sum_{s=1}^S ln(p_{gs})}
where \eqn{p_{gs}} corresponds to the raw \emph{p}-value obtained for gene \emph{g} in a differential
analysis for study \emph{s} (assumed to be uniformly distributed under the null hypothesis). Under the
null hypothesis, the test statistic \eqn{F_g} follows a \eqn{\chi^2}{chi-squared} distribution with \emph{2S}
degrees of freedom. Classical procedures for the correction of multiple testing, such as that of Benjamini
and Hochberg (1995) may subsequently be applied to the obtained \emph{p}-values to control the false
discovery rate at a desired rate \eqn{\alpha}. 
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

M. Brown (1975). A method for combining non-independent, one-sided tests of significance. \emph{Biometrics} \bold{31}(4): 987-992.

A. Rau, G. Marot and F. Jaffrezic (2012). Differential meta-analysis of RNA-seq data. (Manuscript in progress)
}
\seealso{
\code{\link{metaRNASeq}}
}
\examples{
data(rawpval)
fishcomb <- fishercomb(rawpval, BHth = 0.05)
DE <- ifelse(fishcomb$adjpval<=0.05,1,0)
hist(fishcomb$rawpval,nclass=100)

## A more detailed example is given in the vignette of the package:
## vignette("metaRNASeq")
}
\keyword{ methods }
\keyword{ models }