\name{metaRNASeq-package}
\alias{metaRNASeq-package}
\alias{metaRNASeq}
\docType{package}
\title{
Meta-analysis for RNA-seq data.
}
\description{
This package implements two p-value combination techniques (inverse normal and Fisher methods). It also provides a vignette explaining how to perform a meta-analysis from two independent RNA-seq experiments.
}
\details{
\tabular{ll}{
Package: \tab metaRNASeq\cr
Type: \tab Package\cr
Version: \tab 0.4\cr
Date: \tab 2014-02-26\cr
License: \tab GPL \cr
}
}
\author{
Andrea Rau, Guillemette Marot, Florence Jaffrezic 

Maintainer: Guillemette Marot <guillemette.marot@inria.fr>

}
\references{
A. Rau, G. Marot and F. Jaffrezic (2013). Differential meta-analysis of RNA-seq data from multiple studies. arXiv:1306.3636
}
\keyword{ package }
\seealso{
\code{\link{invnorm}}
\code{\link{fishercomb}}
}
\examples{
# See examples in User's Guide:
vignette("metaRNASeq")
}
