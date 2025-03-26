\name{sim.function}
\alias{sim.function}
\alias{extractfromsim}
\title{
Simulation of multiple RNA-seq experiments
}
\description{
Simulate data arising from multiple independent RNA-seq experiments 
}
\usage{
sim.function(param, dispFuncs, nrep = 4, classes = NULL, inter.sd = 0.3)
}
\arguments{
  \item{param}{ Mean expression levels: \code{param} must be a data frame containing at least two columns named 
"mucond1" and "mucond2" and one row per gene.
}
  \item{dispFuncs}{
List of length equal to the number of studies to be simulated, containing the gamma regression parameters 
describing the mean-dispersion relationship for each one; these are the mean-dispersion functions linking mean 
and intra-study variability for each independent experiment. 
}
  \item{nrep}{
Number of replicates to be simulated in each condition in each study. Ignored if \code{classes} is filled.
}
  \item{classes}{
List of class memberships, one per study (maximum 20 studies). Each vector or factor of the list can only 
contain two levels which correspond to the two conditions studied. If NULL, \code{classes} is built as a 
list of two vectors with \code{nrep} labels 1 (for condition 1) and \code{nrep} labels 2 (for condition 2).
}
  \item{inter.sd}{
Inter-study variability. By default, \code{inter.sd} is set to 0.3, which corresponds to a moderate inter-study 
variability in the case where \code{\link{param}} and \code{\link{dispFuncs}} parameters are used to simulate data.
}
}
\value{
A matrix with simulated expression levels, one row per gene and one column per replicate. Names of studies are given 
in the column names of the matrix. 
}
\details{
Details about the simulation procedure are given in the following paper: 
}
\references{
A. Rau, G. Marot and F. Jaffrezic (2014). Differential meta-analysis of RNA-seq data. \emph{BMC Bioinformatics} \bold{15}:91
}
\note{
If the \code{\link{param}} data provided in this package are not used to simulate data, one should check that the 
per-condition means in \code{param} are reasonable. Note also that for genes to be simulated as non-differentially
expressed, the values of "mucond1" and "mucond2" in \code{param} should be equal.
}

\seealso{
\code{\link{metaRNASeq}}
}
\examples{
## Load simulation parameters
data(param)
data(dispFuncs)

## Simulate data
matsim <- sim.function(param = param, dispFuncs = dispFuncs)
sim.conds <- colnames(matsim)
rownames(matsim) <-paste("tag", 1:dim(matsim)[1],sep="")

# extract simulated data from one study
simstudy1 <- extractfromsim(matsim,"study1")
head(simstudy1$study)
simstudy1$pheno
}
\keyword{ methods }
