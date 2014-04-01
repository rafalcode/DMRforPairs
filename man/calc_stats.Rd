\name{calc_stats}
\alias{calc_stats}
\title{CALCULATE DESCRIPTIVE AND TEST STATISTICS}
\description{
Calculate the (difference in) median beta and M values for the rows specified in probe_rows (descriptive statistics). probe_rows must specify a continuous range of adjacent probes. If a regions shows an high difference (>dM_min) between any pair of samples, the differences are formally tested (Mann Whitney U if n=2, Kruskal Wallis n>2). If n>2 and the omnibus test is significant , pairwise testing is performed to determine which differ significantly. 
}
\usage{
calc_stats(probe_rows,probes,m,b,n,dMth,do.format=FALSE)
}
\arguments{
  \item{probe_rows}{Vector of row numbers in regions$valid.m, regions$valid.beta or regions$perprobe identifying adjacent probes. See output \code{\link{regionfinder}}}
  \item{probes}{see output of \code{\link{regionfinder}} (regions$perprobe)}
  \item{m}{see output of \code{\link{regionfinder}} (regions$valid.m)}
  \item{b}{see output of \code{\link{regionfinder}} (regions$valid.beta)}
  \item{n}{total number of samples}
  \item{dMth}{see min_dM parameter of \code{\link{DMRforPairs}}}
  \item{do.format}{see \code{\link{testregion}}}
}
\details{
The p-values reported by this function are not yet adjusted for multiple testing. This function is rarely called by the user directly.
}
\value{
Data frame or vector with statistics and test results concerning the probes specified in probe_rows (please see description).
}
\seealso{
\code{\link{DMRforPairs}}, 
\code{\link{regionfinder}}, 
\code{\link{testregion}}.
}
\author{Martin Rijlaarsdam, \email{m.a.rijlaarsdam@gmail.com}}
\examples{
#calc_stats() is an integrated part of the DMRforPairs() wrapper and is 
#not usually called by the user directly. Please see DMRforPairs() for 
#an example.
}
\keyword{analysis}
