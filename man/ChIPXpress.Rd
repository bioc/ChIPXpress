\name{ChIPXpress}
\alias{ChIPXpress}
\title{
ChIPXpress Ranking
}
\description{
Ranks list of predicted TF bound genes from ChIPx data using a database of publicly available gene expression profiles.
}
\usage{
ChIPXpress(TFID, ChIP, DB)
}
\arguments{
  \item{TFID}{
    A character value corresponding to the Entrez GeneID of the TF-of-interest.
}
  \item{ChIP}{
    A vector of Entrez GeneIDs corresponding to the predicted TF bound genes from ChIPx data. Entrez GeneIDs should be reported as character values.
}
  \item{DB}{a
    A big.matrix of expression values that is used as the database of publicly availiable gene expression profiles. This database can be either built by the user directly using the functions 'build database' and 'clean database', or loaded from a pre-built database in package \emph{ChIPXpressData}. \emph{ChIPXpressData} currently contains two pre-built databases avaliable for analysis purposes - GPL1261 for mouse ChIPx data and GPL570 for human ChIPx data. See example code below on how to load a big.matrix formatted gene expression database.  
}
}
\details{
ChIPXpress works by ranking the TF bound genes predicted from ChIPx data by the absolute correlation between each gene and the TF in the database of gene expression profiles. Genes that are more highly correlated, either positively or negatively, are ranked as the most likely to be an actual TF target gene. Note, the absolute correlation is calculated only using the samples in which the TF expression is above 0, i.e. the mean TF expression across all samples, in order to improve prediction performance. See reference below for more information on the rationale behind the ranking algorithm of ChIPXpress.
}
\value{
    Returns a list with two vectors: the first vector contains the absolute correlations of each predicted TF bound gene with the TF in ranked order, where the names of the vector correspond to the Entrez GeneID of each gene, and the second vector contains the Entrez GeneIDs of the predicted TF bound genes not found in the database.
}
\references{
Wu G. and Ji H. (2012) ChIPXpress: enhanced ChIP-seq and ChIP-chip target gene identification using publicly available gene expression data. \emph{In preparation}.
}
\author{
George Wu
}
\examples{
## Example analyses of real Oct4 bound genes predicted
## from ChIP-seq data in ESC using pre-built GPL1261 
## database

## Load predicted Oct4-bound genes from ChIPx data
data(Oct4ESC_ChIPgenes)

## Load example GPL1261 Database
library(bigmemory)
path <- system.file("extdata",package="ChIPXpressData")
DB_GPL1261 <- attach.big.matrix("DB_GPL1261.bigmemory.desc",path=path)
## To load the human GPL570 data, replace 'DB_GPL1261' with 'DB_GPL570'.

## Run ChIPXpress ("18999" is Entrez GeneID of Oct4)
out <- ChIPXpress(TFID="18999",ChIP=Oct4ESC_ChIPgenes,DB=DB_GPL1261)
head(out[[1]]) ## GeoChIP Rankings
head(out[[2]]) ## Missing genes not found

}
\keyword{ ChIPXpress }
