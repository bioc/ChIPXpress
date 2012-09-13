\name{Oct4ESC_ChIPgenes}
\alias{Oct4ESC_ChIPgenes}
\docType{data}
\title{
Predicted Oct4 bound genes in embryonic stem cells (ESC) obtained from analyzing ChIP-seq data
}
\description{
Vector of Oct4 bound genes predicted by analyzing ESC ChIP-seq data from GSE11724.
}
\usage{data(Oct4ESC_ChIPgenes)}
\format{
 A data frame with 5158 observations on the following 21 variables.
  \describe{
    \item{\code{Rank}}{a numeric vector}
    \item{\code{Chr}}{a character vector}
    \item{\code{Start}}{a numeric vector}
    \item{\code{End}}{a numeric vector}
    \item{\code{Strand}}{a character vector}
    \item{\code{Annotation}}{a character vector}
    \item{\code{Gene}}{a character vector}
    \item{\code{EntrezID}}{a numeric vector}
    \item{\code{peak_length}}{a numeric vector}
    \item{\code{FDR}}{a numeric vector}
    \item{\code{left_peakboundary}}{a numeric vector}
    \item{\code{right_peakboundary}}{a numeric vector}
    \item{\code{peak_summit}}{a numeric vector}
    \item{\code{bound_center}}{a numeric vector}
    \item{\code{bound_width}}{a numeric vector}
    \item{\code{maxT}}{a numeric vector}
    \item{\code{maxT_pos}}{a numeric vector}
    \item{\code{max_log2FC}}{a numeric vector}
    \item{\code{maxFC_pos}}{a numeric vector}
    \item{\code{minuslog10_minPoisP}}{a numeric vector}
    \item{\code{minPoisP_pos}}{a numeric vector}
  }
}
\details{
To obtain the TF bound gene predictions, the ChIP-seq data is processed using CisGenome with the default parameters. Only peaks significant at a FDR of 0.10 are retained and annotated by assigning peaks to genes if the peak falls within 10kbp upstream or 5kbp downstream of the gene transcription start site. Only the highest ranking peak for each gene is retained in the data frame for input. To be clear, the Rank in the data frame corresponds to the original peak ranking by CisGenome. The ChIPx ranking is simply the order of genes in the data frame. 
}
\source{
www.ncbi.nlm.nih.gov/geo/
}
\references{
Marson A. et al. (2008) Connecting microRNA genes to the core trancsriptional regulatory circuitry of embryonic stem cells. \emph{Cell} \bold{134}, 521-533.

Barrett T., et al. (2007) NCBI GEO: mining tens of millions of expression profiles - database and tools update. \emph{Nucl. Acids Res.} \bold{35}, D760-D765.
}
\examples{
data(Oct4ESC_ChIPgenes)
}
\keyword{datasets,Oct4}
