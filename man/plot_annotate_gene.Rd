\name{plot_annotate_gene}
\alias{plot_annotate_gene}
\title{PLOT METHYLATION VALUES AROUND A GENE}
\description{
Wrapper for \code{\link{plot_annotate_probes}}. Generates a plot of the methylation status of / around a specific gene.
}
\usage{
plot_annotate_gene(gs,regions,margin=10000,ID,clr=NA,annotate=TRUE,path)
}
\arguments{
  \item{gs}{Gene Symbol of interest}
  \item{regions}{Output of \code{\link{regionfinder}}}
  \item{margin}{See \code{\link{export_data}}}
  \item{ID}{Custom string which will be used as the file name of the figures generated}
  \item{clr}{See \code{\link{export_data}}}
  \item{annotate}{See \code{\link{plot_annotate_region}}}
  \item{path}{Path relative to the current working directory to save the files (=experiment.name in \code{\link{export_data}})}
}
\details{
Expects a subfolder with the name of your experiment to exist as created by export_data(). Writing permissions are required in the working directory to use the export and visualization functions of DMRforPairs. Internet access is required to use the annotation features.\cr
}
\value{
A number of figures  (including pairwise plots if n>2) is generated visualizing methylation status in the samples at / near the gene investigated. Statistics and optional annotation information are reported back by the function.
}
\seealso{
\code{\link{export_data}}, 
\code{\link{plot_annotate_probes}}, 
\code{\link{regionfinder}}, 
}
\author{Martin Rijlaarsdam, \email{m.a.rijlaarsdam@gmail.com}}
\examples{
#For the this step, please ensure that you have writing 
#permissions in your working directory. If any of the
#annotate parameters is set to TRUE, an internet connection
#is required.

#This example expects a folder "example_DMRforPairs" to exist and
#be writable in the current working directory.
data(DMRforPairs_data)
experiment="example_DMRforPairs"

plot_annotate_gene(gs="DUS4L",
                   regions=example_output_DMRforPairs$regions,
                   margin=10000,
                   ID="DUS4L",
                   annotate=FALSE,path=experiment)			   
}
\keyword{export_visualization}
