\name{buildDatabase}
\alias{buildDatabase}
\title{
Builds a new database of gene expression profiles from a specified platform or sample files in NCBI GEO.
}
\description{
Takes as input a specified platform (GPL ID) or vector of sample files (GSM IDs). All files for the given platform or the vector of given GSM files are then downloaded, processed by fRMA, and stored in a large matrix of expression values in big.matrix format. 

The user will still need to convert the rowIDs from probeIDs to Entrez GeneIDs and then  run 'cleanDatabase' on the resulting matrix in order to modify the databse to the format required for ChIPXpress analyses. 
}
\usage{
buildDatabase(GPL_id, GSMfiles = NULL, SaveDir = NULL, LoadPrevious=FALSE)
}
\arguments{
  \item{GPL_id}{
    A character value of format 'GPLXXX' indicating the platform from which a database of gene expression profiles will be built.
}
  \item{GSMfiles}{
    A vector of character values of format 'GSMXXX' indicating the samples from which a database of gene expression profiles will be built. The GSMfiles must be from the same platform, and if the GSMfiles are provided, GPL_id does not need to be specified.
}
  \item{SaveDir}{
    Path of an EMPTY temporary directory for the built database and cel files to be downloaded to. The user will be required to create the directory beforehand and make sure it is empty. All cel files will be removed after being processed.
}
  \item{LoadPrevious}{
    If LoadPrevious is set to TRUE, then buildDatabase will continue from the last unprocessed sample file. Use this when the buildDatabase is abruptly stopped and unable to finish. Be sure to input the other arguements exactly as they were inputted previously.
}
}
\details{
    Use this function only if you would like to build your own database of gene expression profiles. For the common user interested in ranking TF bound genes from mouse or human ChIPx data, it is easier and faster to load a pre-built database. To do, see the example code in the man page of the ChIPXpress function.

The overall process of creating a database proceeds in three steps:
(1) Run buildDatabase
(2) Convert rownames from probeIDs into EntrezGeneIDs
(3) Run cleanDatabase
The user will be required step 2 themselves. This means annotating the rows by downloading and installing the appropriate annotation package and replacing the probeID rownames with the appropriate gene IDs. Feel free to do this with any package of method that is most convenient. Note, it is possible to annotate the rows using a different ID format than the recommended Entrez GeneID. If you choose to do this, then be sure later, that when running the ChIPXpress analyses that the inputted list of TF bound genes is also of the SAME ID format.

buildDatabase uses GEOquery to download the files and then processes them using frma. Thus, the user will need to have the required frmavecs package installed apriori, so frma can process the raw gene expression files. See the help files from the frma package for more information.

The process of downloading and processing a large collection of gene expression profiles can take an extremely long time. For example, if we wanted to build a database from all of the GPL1261 samples in NCBI GEO, this will take 2 weeks since this would require processing 29086 samples! Be sure to specify an empty directory for the files to be downloaded to, since the function will remove all .gz files after each iteration in order to save hard drive space.

If for any reason the database building process is abruptly stopped, the user can continue from their previous point by running the function with the SAME inputs, save directory, and setting LoadPrevious=TRUE.
}
\value{
    A big.matrix of gene expression values with rows corresponding to the probes on the platform and columns corresponding to the number of samples processed. This is essentially the output from frma() concatenated by column and then formatted as a big.matrix. The database will still need to be annotated by converting the probeIDs into Entrez GeneID format, and then inputted into the cleanDatabase function to obtain a finished database suitable for input into the ChIPXpress function.
}
\references{
McCall M.N., Bolstad B.M., and Irizarry R.A. (2010) Frozen robust multiarray analysis (fRMA). Biostatistics 11, 242-253.

Barrett T., et al. (2007) NCBI GEO: mining tens of millions of expression profiles - database and tools update. Nucl. Acids Res. 35, D760-D765.

Wu G. and Ji H. (2012) ChIPXpress: enhanced ChIP-seq and ChIP-chip target gene identification using publicly available gene expression data. \emph{In preparation}.
}
\author{
George Wu
}
\examples{

\dontrun{
## In this example, we construct a database from 
## samples stored in NCBI GEO on the GPL1261 platform.
## The output would then be used as the database 
## to input into the ChIPXpress function.

## Load required libraries
library(bigmemory)
library(biganalytics)
library(GEOquery)
library(affy)
library(frma)

library(mouse4302frmavecs) 
## GPL1261 corresponds to the mouse430 2.0 array (required by frma)
## If your database is from a different platform, you will
## need to download the corresponding frmavecs package.

SaveDir <- tempdir()
## Make sure the save directory is empty and created.

DB <- buildDatabase(GPL_id='GPL1261',SaveDir=SaveDir)
## This will take up to 2 weeks to finish since
## it processes all GPL1261 samples currently
## stored in NCBI GEO (over 29000 samples).

## Alternatively, the user can also specify by GSM ids.
GSM_ids <- c("GSM24056","GSM24058","GSM24060","GSM24061",
             "GSM94856","GSM94857","GSM94858","GSM94859")
DB <- buildDatabase(GSMfiles=GSM_ids,SaveDir=SaveDir)
## This will only take approximately 5 minutes
## since only 8 GPL1261 samples are specified.
## Note, the database will need to include at least
## 2 samples in which the TF is above the mean TF expression
## to calculate a correlation estimate as required by the
## ChIPXpress algorithm. The more samples the better!!

## Annotates database by converting rowIDs from probeIDs into
## Entrez GeneIDs. Feel free to use any annotation method
## that is most convenient.
library(mouse4302.db)
EntrezID <- mget(as.character(rownames(DB)),mouse4302ENTREZID)
rownames(DB) <- as.character(EntrezID)

## Clean Database - cleans database into format required by ChIPXpress
cleanDB <- cleanDatabase(DB,SaveFile="newDB_GPL1261.bigmemory",
                         SavePath=SaveDir)

head(cleanDB) ## Final database ready for direct input into ChIPXpress

## To load the saved database in new R sessions, simply attach using
## the saved description file of the big.matrix as follows:
cleanDB <- attach.big.matrix("newDB_GPL1261.bigmemory.desc",path=SaveDir)
head(cleanDB)

}

}
\keyword{ buildDatabase }
