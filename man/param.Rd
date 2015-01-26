\name{param}
\alias{param}
\docType{data}
\title{
Mean simulation parameters
}
\description{
Mean simulation parameters obtained from the analysis of a real dataset
}
\usage{data(param)}
\format{
  A data frame with 26408 observations on the following 3 variables.
  \describe{
    \item{\code{mucond1}}{a numeric vector with mean parameters for condition 1}
    \item{\code{mucond2}}{a numeric vector with mean parameters for condition 2}
    \item{\code{DE}}{a logical vector indicating which tags are differentially expressed (value 1)}
  }
}
\details{
Mean parameters provided in this package are empirical means (obtained after normalization for library size differences) of real data described in the following references.
}
\source{
Supplementary material of (Dillies et al., 2013) paper.
}
\references{
M.A. Dillies, A. Rau, J. Aubert, C. Hennequet-Antier, M. Jeanmougin, N. Servant, C. Keime, G. Marot, D. Castel, J. Estelle, G. Guernec, B. Jagla, L. Jouneau, D. Laloe, C. Le Gall, B. Schaeffer, S. Le Crom, M. Guedj and F. Jaffrezic, on behalf of the French StatOmique Consortium (2013) 
A comprehensive evaluation of normalization methods for Illumina high-throughput RNA sequencing data analysis. \emph{Briefings in Bioinformatics} \bold{14}(6):671-83 .

T. Strub, S. Giuliano, T. Ye, et al. (2011) Essential role of microphthalmia transcription factor for DNA replication, mitosis and genomic stability in melanoma. emph{Oncogene} \bold{30}:2319-32.
}
\examples{
data(param)
}
\keyword{datasets}
