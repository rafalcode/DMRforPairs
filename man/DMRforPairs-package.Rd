\name{DMRforPairs-package}
\alias{DMRforPairs-package}
\title{DMRforPairs package}
\docType{package}
\description{
Please see the DESCRIPTION file.
}
\note{
DMRforPairs is independent of the platform used to generate the data and the pipeline used to pre-process it. It can be used with any dataset as long as data frames are provided that have one matching row for each of the m probes.  The examples provided in the vignette and associated publication focus mainly on data aquired using the Infinium HumanMethylation450 BeadChip (450K) platform. Methods to import an preprocess 450K data have been described extensively (see references below). Normalization and preprocessing are readily available in various pipelines and are not part of DMRforPairs. When assassing regionf containing multiple probes, special care should be taken to normalize for probe type bias (for review, see Wilhelm-Benartzi et al 2013). In the following lists some general pointers are given on how to access the data for DMRforPairs using frequenly used and publically avaliable 450K pipelines. An example on how to process data from the CHARM platform is given in the supplementary data of the original DMRforPairs publication (Rijlaarsdam et al 2014, submitted).\cr
\itemize{
	\item{lumi: DMRforPairs is easily applied using 450K data imported by the lumi package (MethyLumiM object). The lumiMethyR() function can be used to read in a final report file as exported from GenomeStudio. The import includes intensity data and annotation information as exported by GenomeStudio. The lumi package offers extensive quality control and pre-processing options.  The estimateBeta() and estimateM() functions should be used to extract beta and M values. The fData() method can be used to access annotation info. This manual focusses on the lumi package and supplies additional pointers on how to extract the necessary information using this package in the documentation of \code{\link{DMRforPairs}}. For more information, please see the \href{http://www.bioconductor.org/packages/release/bioc/html/lumi.html}{LUMI Bioconductor page}}
	\item{IMA: the IMA package also uses the final report exported by GenomeStudio as the source for methylation data and annotation information (IMA.methy450R() function). Pre-processing is done via the IMA.methy450PP() function. beta values are readily available after import, while M values can be computed by calculating the logit2 of beta (M=log2(beta/1-beta)). For more information, please see the \href{http://ima.r-forge.r-project.org/}{IMA homepage}}
	\item{minfi: minfi reads in the intensity information directly from the IDAT files, subsequently using the read.450k.sheet() and read.450k.exp() functions. Minfi also offers a range of QC and pre-processing options. the getBeta() and getM() functions are used to access the methylation quantities needed by DMRforPairs. Since bioconductor 2.13, this package depends on the annotation package IlluminaHumanMethylation450kanno.ilmn12.hg19. For more information, please see the \href{http://www.bioconductor.org/packages/2.12/bioc/html/minfi.html}{minfi} and \href{http://www.bioconductor.org/packages/release/data/annotation/html/IlluminaHumanMethylation450kanno.ilmn12.hg19.html}{IlluminaHumanMethylation450kanno.ilmn12.hg19} Bioconductor pages.}
}}
\note{
Writing permissions  are required in the working directory to use the export and visualization functions of DMRforPairs. Internet access is required to use the annotation features.\cr
The figure below shows the relation between the functions described in this manual.\cr\cr
\figure{FunctionSchedule.png}
}
\references{
\itemize{
\item{Bibikova, M., et al., \emph{High density DNA methylation array with single CpG site resolution.} Genomics, 2011. 98(4): p. 288-95.}
\item{Bibikova, M., et al., \emph{Genome-wide DNA methylation profiling using Infinium(R) assay.} Epigenomics, 2009. 1(1): p. 177-200.}
\item{Dedeurwaerder, S., et al., \emph{Evaluation of the Infinium Methylation 450K technology.} Epigenomics, 2011. 3(6): p. 771-84.}
\item{Du, P., W.A. Kibbe, and S.M. Lin, \emph{lumi: a pipeline for processing Illumina microarray.} Bioinformatics, 2008. 24(13): p. 1547-8.}
\item{Du, P., et al., \emph{Comparison of Beta-value and M-value methods for quantifying methylation levels by microarray analysis.} BMC Bioinformatics, 2010. 11: p. 587.}
\item{Maksimovic, J., L. Gordon, and A. Oshlack, \emph{SWAN: Subset-quantile within array normalization for illumina infinium HumanMethylation450 BeadChips.} Genome Biol, 2012. 13(6): p. R44.}
\item{Marabita, F., et al., \emph{An evaluation of analysis pipelines for DNA methylation profiling using the Illumina HumanMethylation450 BeadChip platform.} Epigenetics, 2013. 8(3): p. 333-46.}
\item{Pidsley, R., et al., \emph{A data-driven approach to preprocessing Illumina 450K methylation array data.} BMC Genomics, 2013. 14(1): p. 293.}
\item{Teschendorff, A.E., et al., \emph{A beta-mixture quantile normalization method for correcting probe design bias in Illumina Infinium 450 k DNA methylation data.} Bioinformatics, 2013. 29(2): p. 189-96.}
\item{Touleimat, N. and J. Tost, \emph{Complete pipeline for Infinium((R)) Human Methylation 450K BeadChip data processing using subset quantile normalization for accurate DNA methylation estimation.} Epigenomics, 2012. 4(3): p. 325-41.}
\item{Wang, D., et al., \emph{IMA: an R package for high-throughput analysis of Illumina's 450K Infinium methylation data.} Bioinformatics, 2012. 28(5): p. 729-30.}
\item{Wilhelm-Benartzi, C.S., et al., \emph{Review of processing and analysis methods for DNA methylation array data.} Br J Cancer, 2013. 109(6): p. 1394-402.}
\item{Eckhardt, F., et al., \emph{DNA methylation profiling of human chromosomes 6, 20 and 22.} Nat Genet, 2006. 38(12): p. 1378-85.}
\item{Rijlaarsdam, M.A., et al., \emph{DMRforPairs: identifying Differentially Methylated Regions between unique samples using array based methylation profiles.} Submitted}
}
}
\author{Martin Rijlaarsdam, \email{m.a.rijlaarsdam@gmail.com}}
