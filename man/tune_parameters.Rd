\name{tune_parameters}
\alias{tune_parameters}
\title{TUNE MIN_DISTANCE AND MIN_N PARAMETERS}
\description{
  This function evaluates the number of regions identified for various settings of min_distance and n_min (see \code{\link{DMRforPairs}}). It also summarizes the number of probes available inclusion using a specific setting of recode and the number of probes that is covered by the detected regions.
}
\usage{
tune_parameters(parameters,classes_gene,classes_island,targetID, chr, position,m.v,beta.v,recode=1,sep=";",gs,do.parallel=0)
}
\arguments{
  \item{parameters}{A data frame with two columns: [1] min_distance and [2] min_n. The number of regions / probes will be evaluated for each pair of parameters (per row).}
  \item{classes_gene}{Please see \code{\link{DMRforPairs}}}
  \item{classes_island}{Please see \code{\link{DMRforPairs}}}
  \item{targetID}{Please see \code{\link{DMRforPairs}}}
  \item{chr}{Please see \code{\link{DMRforPairs}}}
  \item{position}{Please see \code{\link{DMRforPairs}}}
  \item{m.v}{Please see \code{\link{DMRforPairs}}}
  \item{beta.v}{Please see \code{\link{DMRforPairs}}}
  \item{recode}{Please see \code{\link{DMRforPairs}}. (Default=1)}
  \item{sep}{Please see \code{\link{DMRforPairs}}. (Default=";")}
  \item{gs}{Please see \code{\link{DMRforPairs}}}
  \item{do.parallel}{Tuning of the parameters requires a significant amount of compute time and power. Parallelization if offered to do this more efficiently. 0=no parallelization, -1=use all available cores, n>1 use n cores (default=0)}
}
\details{
This function requires a significant amount of compute time/power. Parallelization is highly recommended. When parallelization is enabled, no progress indicator is displayed by R.
}
\value{
A data frame is returned with the folowing columns
  \item{min_distance}{The setting of min_distance. Also see \code{\link{DMRforPairs}}}
  \item{min_n}{The setting for min_n. Also see \code{\link{DMRforPairs}}}
  \item{n.regions}{The number of regions identified in your data using the before mentioned combination of settings and the setting of recode.}
  \item{n.valid.probes}{The number of probes that were assigned to one or more classes in the \code{\link{merge_classes}} function. This is constant for a specific data set as long as the recode parameter is the same.}
  \item{n.probes.included}{The number of probes included in the identified regions.}
}
\seealso{
  \code{\link{regionfinder}}, 
  \code{\link{DMRforPairs}}
}
\author{Martin Rijlaarsdam}
\examples{
data(DMRforPairs_data)
CL.methy=CL.methy[which(CL.methy$position<=1.07E+8 & 
                        CL.methy$position>=1.06E+8),]
#just two relevant values for the parameters 
#min_distance (100 & 300) and min_n (4 & 8)
parameters=expand.grid(min_distance = c(200), min_n = c(4,5))
results_parameters=	tune_parameters(parameters,
					classes_gene=CL.methy$class.gene.related, 
					classes_island=CL.methy$class.island.related, 
					targetID=CL.methy$targetID, 
					chr=CL.methy$chromosome, 
					position=CL.methy$position, 
					m.v=CL.methy[,c(7:8)], 
					beta.v=CL.methy[,c(11:12)],
					recode=1,
					gs=CL.methy$gene.symbol,
					do.parallel=0)
results_parameters  
}
\keyword{analysis}