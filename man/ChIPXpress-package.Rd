\name{ChIPXpress-package}
\alias{ChIPXpress-package}
\docType{package}
\title{
ChIPXpress ranks TF target genes given predicted TF targets from ChIP-chip and ChIP-seq data using publicly available gene expression data
}
\description{
Takes as input predicted TF bound targets from ChIP-chip or ChIP-seq data and ranks them according to the most likely to be a TF target gene. ChIPXpress rankings are more consistent, robust, and accurate than standard ChIP-chip or ChIP-seq rankings, which prioritize genes only on the strength of the observed peaks. ChIPXpress is able to accomplish this because it searches for TF bound genes that are highly correlated in expression with the TF across a database of highly diverse gene expresison profiles collected from different diseases, tissues, and cell types.
}
\details{
\tabular{ll}{
Package: \tab ChIPXpress\cr
Type: \tab Package\cr
Version: \tab 0.99.0\cr
Date: \tab 2012-07-31\cr
License: \tab GPL-2\cr
}
}
\author{
Author: George Wu
Maintainer: George Wu <gewu@jhsph.edu>
}
\references{
Wu G. and Ji H. (2012) ChIPXpress: enhanced ChIP-seq and ChIP-chip target gene identification using publicly available gene expression data. \emph{In preparation}.
}
\keyword{ package, ChIPXpress }
