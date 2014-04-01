\name{CL.methy}
\docType{data}
\alias{CL.methy}
\title{METHYLATION DATA OF 4 CELL LINES (CHR 7)}
\description{
  This dataset provides the average methylation values on chromosomes 7 (33-35MG, 105-108MB) of two commercially available EBV transfected lymphoblastoid cell lines from healthy individuals (\href{http://ccr.coriell.org/Sections/Search/Sample_Detail.aspx?Ref=GM17105}{NA17105} (African American male) and \href{http://ccr.coriell.org/Sections/Search/Sample_Detail.aspx?Ref=GM17018}{NA17018} (Chinese female), Coriell Institute for Medical Research (NJ, USA)). The dataset also contains the breast cancer cell line MCF7 (Soule et al. 1973) and HPV negative squamous-cell vulva carcinoma cell line A431 (Giard et al. 1973 & Hietanen et al. 1995). Raw data is publically avaliable via the \href{http://support.illumina.com/downloads/genomestudio_software_20111.ilmn}{Illumina website} and discussed in Bibikova et al. 2009 & 2011. Import and pre-processing was carried out using the LUMI package (Du et al. 2008) in R following the optimized "lumi: QN+BMIQ" pipeline (Marabita et al. 2013). This includes exclusion of poorly performing probes (p<0.01, n=713), color adjustment, quantile normalization and correction for probe type bias (Infinium I vs II) using the BMIQ algorithm (Teschendorf 2013).
}
\usage{CL.methy}
\format{A data frame of 14 variables containing 711 observations. Variables: targetID (probe name), chromosome, position (bp), class.gene.related (functional class related to gene), class.island.related (functional class related to CpG island), gene.symbol, A431.M, MCF7.M, NA17105.M, NA17018.M, A431.beta, MCF7.beta, NA17105.beta, NA17018.beta}
\source{Bibikova et al 2009 & 2011, \href{http://support.illumina.com/downloads/genomestudio_software_20111.ilmn}{download from Illumina website}.}
\references{
\itemize{
\item{Soule, H.D., et al., \emph{A human cell line from a pleural effusion derived from a breast carcinoma.} J Natl Cancer Inst, 1973. 51(5): p. 1409-16.}
\item{Giard, D.J., et al., \emph{In vitro cultivation of human tumors: establishment of cell lines derived from a series of solid tumors.} J Natl Cancer Inst, 1973. 51(5): p. 1417-23.}
\item{Hietanen, S., et al., \emph{Human papillomavirus in vulvar and vaginal carcinoma cell lines.} Br J Cancer, 1995. 72(1): p. 134-9.}
\item{Bibikova, M., et al., \emph{High density DNA methylation array with single CpG site resolution.} Genomics, 2011. 98(4): p. 288-95.}
\item{Bibikova, M., et al., \emph{Genome-wide DNA methylation profiling using Infinium(R) assay.} Epigenomics, 2009. 1(1): p. 177-200.}
\item{Du, P., W.A. Kibbe, and S.M. Lin, \emph{lumi: a pipeline for processing Illumina microarray.} Bioinformatics, 2008. 24(13): p. 1547-8.}
\item{Marabita, F., et al., \emph{An evaluation of analysis pipelines for DNA methylation profiling using the Illumina HumanMethylation450 BeadChip platform.} Epigenetics, 2013. 8(3): p. 333-46.}
\item{Teschendorff, A.E., et al., \emph{A beta-mixture quantile normalization method for correcting probe design bias in Illumina Infinium 450k DNA methylation data.} Bioinformatics, 2013. 29(2): p. 189-96.}
}
}
\keyword{datasets}