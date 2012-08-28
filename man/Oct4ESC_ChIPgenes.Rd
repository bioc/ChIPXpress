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
  The format is:
 chr [1:5158] "235345" "22059" "51869" "21429" "109113" "64297" ...
}
\details{
To obtain the TF bound gene predictions, the ChIP-seq data is processed using CisGenome with the default parameters. Only peaks significant at a FDR of 0.10 are retained and annotated by assigning peaks to genes if the peak falls within 10kbp upstream or 5kbp downstream of the gene transcription start site. 
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
