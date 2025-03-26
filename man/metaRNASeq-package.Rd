\name{metaRNASeq-package}
\alias{metaRNASeq-package}
\alias{metaRNASeq}
\docType{package}
\title{
Meta-analysis for RNA-seq data.
}
\description{
Implementation of two p-value combination techniques (inverse normal and Fisher methods). A vignette is provided to explain how to perform a meta-analysis from two independent RNA-seq experiments.
}
\details{
\tabular{ll}{
Package: \tab metaRNASeq\cr
Type: \tab Package\cr
Version: \tab 1.0.2\cr
Date: \tab 2015-01-26\cr
License: \tab GPL \cr
}
}
\author{
Andrea Rau, Guillemette Marot, Florence Jaffrezic 

Maintainer: Guillemette Marot <guillemette.marot@inria.fr>

}
\references{
A. Rau, G. Marot and F. Jaffrezic (2014). Differential meta-analysis of RNA-seq data. \emph{BMC Bioinformatics} \bold{15}:91
}
\keyword{ package }
\seealso{
\code{\link{invnorm}}
\code{\link{fishercomb}}
}
\examples{
#An User's guide with detailed examples can be downloaded in interactive R sessions
if(interactive()){
	vignette("metaRNASeq")
	}
}
