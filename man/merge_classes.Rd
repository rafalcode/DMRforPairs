\name{merge_classes}
\alias{merge_classes}
\title{MERGE ILLUMINA PROBE CLASSES}
\description{
Illumina assigns quite specific functional classes to the probes. 11 classes are defined (Bibikova et al. 2009 & 2011):
\itemize{
\item{Relation to gene: \emph{Body, 5'UTR, 3'UTR, 1stExon, TSS1500, TSS200}}
\item{Relation to CpG island: \emph{Island, N_Shelf, N_Shore, S_Shelf, S_Shore}}
}
As DMRforPairs looks at probes in close proximity to each other within each class, very specific annotations might result in too few probes per region per class and a subsequent drop in the number of identified regions and/or statistical power. This function therefore allows grouping and / or selection of classes of interest.
}
\usage{
merge_classes(refgene_class,island_class,recode=1,sep=";")
}
\arguments{
  \item{refgene_class}{see description and classes_gene parameter in \code{\link{DMRforPairs}}}
  \item{island_class}{see description and classes_island parameter in \code{\link{DMRforPairs}}}
  \item{recode}{recoding scheme to use for the functional classes. Can be custom scheme (data frame) or build in scheme (0, 1 or 2). See details.}
  \item{sep}{Separator used in the second column of the recode parameter. Use ";" or do not specify if using the build in schemes.}
}
\details{
The recode parameter can be set to use one of the build in recoding schemes:
\itemize{
\item{0: analyze all 11 classes annotated by Illumina separately}
\item{1: group Body, 5'UTR, 3'UTR into one category "gene" and TSS1500 and TSS200 in another ("tss"). All island associated classes are merged in one class "island"}
\item{2: all probes without subdivision into classes (also included probes associated with no class).}
}
The recode parameter can also be set to a custom recoding scheme (data.frame). For example: \cr\cr
\emph{data.frame(c("gene","tss"),c("Body;5\'UTR;3\'UTR;1stExon","TSS1500;TSS200"))}. \cr\cr
In this scheme the classes are merged into two categories: TSS or other gene region. Probes solely associated with CpG island-related classes are discarded. Probes not annotated to any of the 11 classes are always discarded by DMRforPairs, except when option 2 is used which collects all probes into one class (i.e. ignores classes).\cr\cr
If classes are unknown an m x 1 character vector with "unknown.gene" and "unknown.island" for all m rows (probes) can be used for the refgene_class and island_class parameters respectively (set "recode" to 2!).
}
\value{
List of objects containing:
\item{$pclass}{original classes per probe (gene and island classes are merged) (m x 1 data frame)}
\item{$pclass_merged}{classes after recoding (m x 1 data frame)}
\item{$no.pclass}{row indexes of probes with no annotation to any of the classes specified in the recoding scheme.}
\item{$u_pclass}{unique list of the classes of interest after recoding (i.e. the first column of the recode data fram)}
}
\references{
\itemize{
\item{Bibikova, M., et al., \emph{High density DNA methylation array with single CpG site resolution.} Genomics, 2011. 98(4): p. 288-95.}
\item{Bibikova, M., et al., \emph{Genome-wide DNA methylation profiling using Infinium(R) assay.} Epigenomics, 2009. 1(1): p. 177-200.}
}
}
\author{Martin Rijlaarsdam, \email{m.a.rijlaarsdam@gmail.com}}
\examples{
#merge_classes() is an integrated part of the DMRforPairs() wrapper and is 
#not usually called by the user directly. Please see DMRforPairs() for 
#an example.
}
\keyword{analysis}
