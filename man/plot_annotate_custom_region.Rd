\name{plot_annotate_custom_region}
\alias{plot_annotate_custom_region}
\title{PLOT METHYLATION OF A CUSTOM GENOMIC REGION}
\description{
Wrapper for \code{\link{plot_annotate_probes}}. Generates a plot of the methylation status of a custom genomic region.
}
\usage{
plot_annotate_custom_region(chr,st,ed,regions,margin=10000,ID="CustomRegion",clr=NA,annotate=TRUE,path)
}
\arguments{
  \item{chr}{Chromosome on which the custom region is located. Number/X/Y only (i.e. no 'chr' prefix.}
  \item{st}{Start position (bp) of the custom region}
  \item{ed}{End position (bp) of the custom region}
  \item{regions}{See output of \code{\link{regionfinder}}}
  \item{margin}{See \code{\link{export_data}}}
  \item{ID}{String which will be used as the file name of the figures generated. (default = "CustomRegion")}
  \item{clr}{See \code{\link{export_data}}}
  \item{annotate}{See \code{\link{plot_annotate_region}}}
  \item{path}{Path relative to the current working directory to store the files (=experiment.name in \code{\link{export_data}})}
}
\details{
A plot will be generated using only the probes present in the \code{\link{DMRforPairs}} output (i.e. associated to one or more of the recoded classes, see \code{\link{merge_classes}}).\cr
Expects a writable subfolder of the curent working directory (as created by export_data()). Internet access is required to use the annotation features.\cr
}
\value{
A number of figures (including pairwise plots if n>2) is generated, visualizing methylation status in the samples at / near the custom region investigated. Statistics and optional annotation information are reported back by the function.
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

#visualize a custom region on the genome
plot_annotate_custom_region(chr=7,
                            st=1.05E+8,
                            ed=1.06E+8,
                            example_output_DMRforPairs$regions,
                            margin=500,
                            ID="custom_region",
                            annotate=FALSE,path=experiment)

}
\keyword{export_visualization}