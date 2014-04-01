\name{plot_annotate_probes}
\alias{plot_annotate_probes}
\title{PLOT METHYLATION+ANNOTATE A SET OF PROBES}
\description{
Function that visualizes the methylation status of a selection of adjacent probes and (if requested) identifies transcripts associated with that region. If n>2 plots are generated for all samples together as well as in a pairwise fashion.
}
\usage{
plot_annotate_probes(regions,title_x,probe_rows,margin=10000,ID=NA,clr=NA,annotate=TRUE,scores=NA,path)
}
\arguments{
  \item{regions}{See output of \code{\link{regionfinder}}.}
  \item{title_x}{String which will be used as the main title of the scatter plots of the beta and M values.}
  \item{probe_rows}{Row numbers of the adjacent probes of interest in the regions object. See output of \code{\link{regionfinder}}:$perprobe, $valid.probes, $valid.m, $valid.beta.}
  \item{margin}{see \code{\link{export_data}}}
  \item{ID}{String which will be used as the file name of the figures generated. (default = "CustomRegion")}
  \item{clr}{see \code{\link{export_data}}}
  \item{annotate}{see \code{\link{plot_annotate_region}}}
  \item{scores}{see \code{\link{plot_annotate_region}}}
  \item{path}{Path relative to the current working directory to save the files (=experiment.name in \code{\link{export_data}})}
}
\details{
A plot will be generated using only probes present in the DMRforPairs output (i.e. associated to one or more of the recoded classes, see \code{\link{merge_classes}}). This function is rarely called directly by the user. There are a number of wrappers available which lookup the probes (probe_rows) associated with a regionID from DMRforPairs, a custom genomic region or a gene of interest and then call plot_annotate_probes().
}
\value{
A number of figures (including pairwise plots if n>2) is generated visualizing methylation status in the samples at / near the probes investigated. Statistics and optional annotation information are reported back by the function.
}
\seealso{
\code{\link{export_data}}, 
\code{\link{plot_annotate_region}}, 
\code{\link{regionfinder}}, 
\code{\link{merge_classes}}.
}
\author{Martin Rijlaarsdam, \email{m.a.rijlaarsdam@gmail.com}}
\examples{
#This function is called by the other plot_annotate_...() wrappers and is 
#not usually called by the user directly. 
#Please see these functions for examples.  
}
\keyword{export_visualization}
