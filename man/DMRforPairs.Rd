\name{DMRforPairs}
\alias{DMRforPairs}
\title{WRAPPER FOR THE DMRforPairs ANALSYSIS}
\description{
  Wrapper for the DMRforPairs analysis. Includes recoding of the probe classes, identification of sufficiently probe-dense regions and analysis (calculation of statistics and testing).
}
\usage{
DMRforPairs(classes_gene,classes_island,targetID, chr, position,m.v,beta.v,min_n=4,min_distance=200,min_dM=1.4,recode=1,sep=";",method="fdr",debug.v=FALSE,gs,do.parallel=0)
}
\arguments{
  \item{classes_gene}{m x 1 data frame with the relation to gene for each probe. Illumina annotates this as Body, 5'UTR, 3'UTR, 1stExon, TSS1500 or TSS200. When the lumiMethyR function of the lumi package is used to import the data from GenomeStudio's final report, this information can be extracted via the fData(function), e.g. fData(data) $UCSC_REFGENE_GROUP where data is the MethyLumiM object resulting from the import. If these classes are unknown an m x 1 character vector with "unknown.gene" can be supplied and recode should be set to 2.}
  \item{classes_island}{m x 1 data frame with the relation to CpG island. Illumina annotates this as Island, N_Shelf, N_Shore, S_Shelf or S_Shore. Information can be accessed via (fData(data)$RELATION_TO_UCSC_CPG_ISLAND. If these classes are unknown an m x 1 character vector with "unknown.island" can be supplied and recode should be set to 2.}
  \item{targetID}{m x 1 data frame containing the identifier of each of the m probes, i.e. cgxxxxxxxx. Information can be accessed via (fData(data)$TargetID}
  \item{chr}{m x 1 data frame containing the chromosome each probe is annotated to. Information can be accessed via (fData(data)$CHR}
  \item{position}{m x 1 data frame containing the genomic position each probe is annotated to. Information can be accessed via (fData(data)$MAPINFO}
  \item{m.v}{m x n matrix or data frame containing the M-values for each probe for each sample. These can be directly extracted from the MethyLumiM using the exprs(estimateM()) functions form the lumi package. Alternatively any m x n matrix can be used (for example output from an external normalization algorithm that does not use the MethyLumiM format or methylation values from another platform.}
  \item{beta.v}{See m.v. Use and exprs(estimateBeta(data)) to extract beta values from the MethyLumiM object.}
  \item{min_n}{Minimal number of probes required to consider a group of subsequent probes for inclusion as a region (=potential DMR). Default (minimum) is 4 as the Mann Whitney U test requires a minimum of 7 samples to ever reach a value < 0.05 (2x4=8).}
  \item{min_distance}{Maximal distance between adjacent probes to accept when including them in one region. The default value of 200bp is based on the findings in Marimba et al. 2013 en Eckhardt et al 2006 regarding rapid drop of co-methylation of adjacent probes when these are further apart.}
  \item{min_dM}{Minimal median difference between M-values in a region to consider the region for formal testing. Default value of 1.4 is based on the findings in Du et al 2010.}
  \item{recode}{Recoding scheme to group or discard probes annotated to certain functional genomic regions (see also classes_gene and classes_island parameters and the \code{\link{merge_classes}} function). (default=1)}
  \item{sep}{Separator used in the second column of the recode parameter. Use ";" or do not specify if using the built-in schemes. (default=";")}
  \item{method}{Method to use for correction for multiple testing. See p.adjust() function in R for possible settings. Default is 'fdr' implicating correction according to Benjamini-Hochberg.}
  \item{debug.v}{For development. If TRUE, only the first chromosome is analyzed. (default=FALSE)}
  \item{gs}{m x 1 data frame with associated gene symbols. From a MethyLumiM object this information can be extracted via the fData(function), e.g. fData(data) $UCSC_REFGENE_NAME}
  \item{do.parallel}{Use parallel processes to compute statistics per region. 0=no parallelization, -1=use all available cores, n>1 use n cores (default=0)}
}
\details{
This wrapper subsequently:
\enumerate{
\item{Recodes probe classes according to a custom or build in scheme. See recode parameter and \code{\link{merge_classes}}.}
\item{Identifies regions with sufficient probe density (i.e. number of probes and proximity) over all genomic regions at which probes are annotated in the dataset. See \code{\link{regionfinder}}}
\item{Calculates relevant statistics (e.g. median (difference in) M and beta values). If the median difference is sufficiently large (>=min_dM), the function performs formal testing of the difference. If indicated, pairwise testing is performed as well. See \code{\link{testregion}} and \code{\link{calc_stats}}.}
}
Export & visualization are not directly included in the wrapper, as these require an internet connection and the analysis itself might take place on a compute server that is not connected to the internet. See \code{\link{export_data}} for more information about this. This wrapper does not store the results, but returns a list of data structures to be processed further or saved by the user.
}
\value{
Returns a list of objects resulting from the different steps of the DMRforPairs  algorithm. Please see the description of the associated functions for more information.
	\item{$classes}{includes the results from \code{\link{merge_classes}}}
	\item{$regions}{includes the results from \code{\link{regionfinder}}}
	\item{$tested}{includes the results from \code{\link{testregion}}}
}
\seealso{
  \code{\link{merge_classes}}, 
  \code{\link{regionfinder}}, 
  \code{\link{testregion}}, 
  \code{\link{calc_stats}}.
}
\references{
\itemize{
\item{Marabita, F., et al., \emph{An evaluation of analysis pipelines for DNA methylation profiling using the Illumina HumanMethylation450 BeadChip platform.} Epigenetics, 2013. 8(3): p. 333-46.}
\item{Eckhardt, F., et al., \emph{DNA methylation profiling of human chromosomes 6, 20 and 22.} Nat Genet, 2006. 38(12): p. 1378-85.}
\item{Du, P., et al., \emph{Comparison of Beta-value and M-value methods for quantifying methylation levels by microarray analysis.} BMC Bioinformatics, 2010. 11: p. 587.}
}
}
\author{Martin Rijlaarsdam, \email{m.a.rijlaarsdam@gmail.com}}
\examples{
#load vignette data & select a region of 4 MB on chr 7 for speed
data(DMRforPairs_data)
CL.methy=CL.methy[which(CL.methy$position<=1.07E+8 & 
                        CL.methy$position>=1.06E+8),]

output=DMRforPairs(classes_gene=CL.methy$class.gene.related, 
                classes_island=CL.methy$class.island.related, 
                targetID=CL.methy$targetID, 
                chr=CL.methy$chromosome, 
                position=CL.methy$position, 
                m.v=CL.methy[,c(7:8)], 
                beta.v=CL.methy[,c(11:12)],gs=CL.methy$gene.symbol,
				do.parallel=0)

#primary output of merge_classes()
head(output$classes$pclass) #orginal probe classes
head(output$classes$pclass_recoded) #recoded probe classes
head(output$classes$no.pclass) #row numbers of probes without a recoded class
head(output$classes$u_pclass) #classes used for recoding

#primary output of regionfinder()
#output listing the potential regions of interest identified
head(output$regions$boundaries)
#probes with associated class after recoding (=valid)
head(output$regions$valid.probes) 
#... and associated m values for all samples
head(output$regions$valid.m) 
#... and associated beta values for al samples
head(output$regions$valid.beta) 
#matrix of valid probes (rows) and recoded probe classes (columns) with 
#either NA if not included in any potential region of interest or the 
#ID of the region the probe is assigned to.
head(output$regions$perprobe,10) 

#primary output of testregion() and calc_stats()
#these results are similar to output$regions$boundaries 
#but are supplemented with descriptive statistics 
#and formal test results per region.
head(output$tested)
}
\keyword{analysis}
