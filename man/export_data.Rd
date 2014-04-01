\name{export_data}
\alias{export_data}
\title{COMPLETE EXPORT AND VISUALIZATION OF DMRforPairs RESULTS}
\description{
Complete export of the DMRforPairs results. Export includes pdf's, thumbnails and TSV files describing the methylation status of all relevant (median delta M > dM_min) regions. M and beta values plotted against genomic position. Plots are generated for all samples and sample pairs (if n>2). Annotation information can be automatically looked up in the Ensembl database. If this is requested, integrated figures of the methylation pattern are generated including all transcripts annotated to the region of interest (+margin). This function also generates overviews (tables) in tsv and HTML format with information about the regions and their methylation status in the samples. This is done separately for all regions, regions with a relevant difference and regions with a significant difference.
}
\usage{
export_data (tested,regions,th=0.05,annotate.relevant=FALSE,annotate.significant=TRUE,FigsNotRelevant=FALSE, min_n=4,min_dM=1.4,min_distance=200,margin=10000,clr=NA,method="fdr", experiment.name,debug=FALSE)
}
\arguments{
  \item{tested}{output of \code{\link{testregion}}}
  \item{regions}{output of \code{\link{regionfinder}}}
  \item{th}{p-value threshold for a DMR to be included as significant. (default=0.05; after correction for multiple testing)}
  \item{annotate.relevant}{By default annotations are only looked up for regions showing significant differential methylation patterns. Setting this parameter to TRUE will query Ensembl for the annotation info of each of the relevant (but not significant) regions as well. Please note that the database queries are slow and Ensembl might not accept an unlimited number of requests. (default=FALSE)}
  \item{annotate.significant}{By default annotations are only looked up for regions showing significant differential methylation patterns. Setting this parameter to FALSE will disable this feature (relevant on machines without an internet connection)(default=TRUE)}
  \item{FigsNotRelevant}{If TRUE, export_data() also exports figures (no annotation) for regions that do not show a relevant median difference in M value. Note that this will result in an extremely high number of files in the figures directory. (default=FALSE)}
  \item{min_n}{see \code{\link{DMRforPairs}}}
  \item{min_dM}{see \code{\link{DMRforPairs}}}
  \item{min_distance}{see \code{\link{DMRforPairs}}}
  \item{margin}{Margin in bp outside of each region that is queried in Ensembl and subsequently plotted). (default=10kb)}
  \item{clr}{Defines the colors to be used when plotting the different samples. If NA, colors are picked from the rainbow pallet. (default= NA)}
  \item{method}{see \code{\link{DMRforPairs}}}
  \item{experiment.name}{String. Is used to create a folder to export the results.}
  \item{debug}{see \code{\link{DMRforPairs}}}
}
\details{
In the output is organized according to all, relevant and significant regions:
\enumerate{
\item{all: all regions identified by DMRforPairs}
\item{relevant: regions identified by DMRforPairs with a median difference in M values between any pair of samples of > min_dM}
\item{significant: relevant + Kruskall Wallis test shows a significant result (<th)}
}
Writing permissions  are required in the working directory when using the export and visualization functions of DMRforPairs. Internet access is required to use the annotation features.\cr
By default, DMRforPairs creates a folder (experiment.name) within the current working directory in which the output is stored. This is done within the \code{\link{export_data}} wrapper which should be run before any other (specific) analysis are performed (exploration phase, see \code{\link{DMRforPairs-package}}. A complete export generates a large number of files which are stored out of sight (experiment.name/figures). Figures and overviews of statistics are most conveniently accessed through the HTML tables.
}
\value{
The output is comparable to that of \code{\link{testregion}}, supplemented with external links (Ensembl, UCSC) and annotation info if requested.
}
\seealso{
\code{\link{DMRforPairs}}, 
\code{\link{regionfinder}}, 
\code{\link{testregion}},
\code{\link{plot_annotate_region}},
\code{\link{plot_annotate_probes}}.
}
\author{Martin Rijlaarsdam, \email{m.a.rijlaarsdam@gmail.com}}
\examples{
#For the this step, pleas ensure that you have writing 
#permissions in your working directory. If any of the
#annotate parameters is set to TRUE, an internet connection
#is required.
data(DMRforPairs_data)
experiment="example_DMRforPairs"

#Export the results from output using the default settings;
#annotation is disabled as to prevent problems on (build)
#systems without an internet connection. A folder 
#"example_DMRforPairs" containing the resuls is created in 
#your working directory.
tested_inclannot=export_data(tested=example_output_DMRforPairs$tested, 
regions=example_output_DMRforPairs$regions, 
experiment.name=experiment,
annotate.significant=FALSE,
annotate.relevant=FALSE)
}
\keyword{export_visualization}
