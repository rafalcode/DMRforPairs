\name{testregion}
\alias{testregion}
\title{METHYLATION DIFFERENECE & TESTING FOR 1 REGIONID}
\description{
Serves as a wrapper for \code{\link{calc_stats}}. \code{\link{calc_stats}} is a generic function using row numbers of adjacent probes which are to be queried for their difference in methylation levels between samples. testregion() identifies the row numbers of probes included in a specific regionID x and then calls \code{\link{calc_stats}}. It returns the difference in methylation level for region x and (if appropriate) the results of formal tests that evaluate this difference.
}
\usage{
testregion(x,probes,m,b,n,dMth,do.format=FALSE)
}
\arguments{
  \item{x}{region ID, see output \code{\link{regionfinder}}.}
  \item{probes}{see output of \code{\link{regionfinder}} (regions$perprobe).}
  \item{m}{see output of \code{\link{regionfinder}} (regions$valid.m).}
  \item{b}{see output of \code{\link{regionfinder}} (regions$valid.beta).}
  \item{n}{total number of samples.}
  \item{dMth}{see min_dM parameter of \code{\link{DMRforPairs}}}
  \item{do.format}{If TRUE, the function returns a formatted data frame with statistics regarding the region. If FALSE, it returns a basic vector. This is more computationally efficient if testregion is used in a loop (like in \code{\link{DMRforPairs}}).}
}
\details{
This function is rarely called by the user directly.
}
\value{
Data frame or vector (depending on the setting of do.format) with statistics and test results concerning region x. For specification, please see \code{\link{calc_stats}}.
}
\seealso{
\code{\link{calc_stats}}
}
\author{Martin Rijlaarsdam, \email{m.a.rijlaarsdam@gmail.com}}
\examples{
#testregion() is an integrated part of the DMRforPairs() wrapper and is 
#not usually called by the user directly. Please see this function for 
#an example.
}
\keyword{analysis}
