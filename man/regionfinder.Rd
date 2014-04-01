\name{regionfinder}
\alias{regionfinder}
\title{IDENTIFICATION OF PROBE DENSE REGIONS}
\description{
Identifies genomic regions with sufficient probe density, i.e. a high number of probes positioned closely together. 
}
\usage{
regionfinder(targetID,chr,position,pclass,r_excl,u_pclass,d=200,m.v,beta.v,n_min=4, debug=FALSE,gs)
}
\arguments{
  \item{targetID}{See \code{\link{DMRforPairs}} parameters}
  \item{chr}{See \code{\link{DMRforPairs}} parameters}
  \item{position}{See \code{\link{DMRforPairs}} parameters}
  \item{pclass}{See \code{\link{merge_classes}} output}
  \item{r_excl}{See \code{\link{merge_classes}} output ($no.pclass). If a probe is not annotated to one of the classes of interest in the recoding scheme, that probe is discarded.}
  \item{u_pclass}{See \code{\link{merge_classes}} output}
  \item{d}{See \code{\link{DMRforPairs}} parameters}
  \item{m.v}{See \code{\link{DMRforPairs}} parameters}
  \item{beta.v}{See \code{\link{DMRforPairs}} parameters}
  \item{n_min}{See \code{\link{DMRforPairs}} parameters}
  \item{debug}{See \code{\link{DMRforPairs}} parameters}
  \item{gs}{See \code{\link{DMRforPairs}} parameters}
}
\details{
Regions are identified (per unique functional class) based on a minimal number of probes that lay in close proximity to each other. This can lead to identical (start to end) regions in different classes due to overlapping annotations. regionfinder() merges such regions to prevent duplications. These regions can be identified by a combined class name (e.g. "tss;island").
}
\value{
List of objects containing
\item{$boundaries}{information about the regions identified (chr, start and end bp position, length of region, n of probes in region and information about the functional class(es) the region is associated to.}
\item{$perprobe}{m x p matrix with region IDs. Rows represent the valid probes (see below) and there is one column per functional class (see u_pclass parameter). Each probe can be associated to one region within each class. Region IDs are specific to a dataset and to a set of DMRforPairs parameters. They are not interchangeable between datasets/experiments and serve mainly as an identifier during exploration of the dataset.}
\item{$valid.m}{m x n matrix with M values for all samples. Only probes associated with at least one class of interest are included (=valid probes).}
\item{$valid.beta}{see valid.m. This matrix contains the beta values.}
}
\seealso{
\code{\link{DMRforPairs}}, \code{\link{merge_classes}}.
}
\author{Martin Rijlaarsdam, \email{m.a.rijlaarsdam@gmail.com}}
\examples{
#regionfinder() is an integrated part of the DMRforPairs() wrapper and is 
#not usually called by the user directly. Please see DMRforPairs() for 
#an example.
}
\keyword{analysis}
