\name{cleanDatabase}
\alias{cleanDatabase}
\title{
Cleans the annotated output from buildDatabase to the required format for ChIPXpress analyses
}
\description{
Given the output matrix of 'buildDatabase' after the rows have been properly annotated (preferably into Entrez GeneIDs), this function will clean the database such that each row is normalized and in 1-to-1 correspondance with a single Entrez GeneID or alternative gene ID format.
}
\usage{
cleanDatabase(DB, SaveFile="newDB.bigmemory", SavePath=".")
}
\arguments{
  \item{DB}{
    A big.matrix of expression values from 'buildDatabase' after the rows of the matrix are converted into Entrez GeneID format. Alternative gene ID formats for the row names is also suitable.
}
  \item{SaveFile}{
    Specifies the name of the big.matrix and big.matrix description file. The default file name is "newdb.bigmemory". This will be used to load the database for future R sessions.
}
  \item{SavePath}{
    Specifies the directory in which the final big.matrix and big.matrix descrption file is saved. The default is to save into the current working directory.
}
}
\details{
    This function is to be used after buildDatabase and row annotation has already been completed. Specifically, the entire process of creating a database proceeds in the following three steps:
(1) Run buildDatabase
(2) Convert rownames from probeIDs into EntrezGeneIDs
(3) Run cleanDatabase
The user will be required to complete step 2 themselves. This means annotating the rows by downloading and installing the appropriate annotation package and replacing the probeID rownames with the appropriate gene IDs. Feel free to do this with any package or method that is most convenient. Note, it is possible to annotate the rows using a different ID format than the recommended Entrez GeneID. If you choose to do this, then be sure later, that when running the ChIPXpress analyses that the inputted list of TF bound genes is also of the SAME ID format.

cleanDatabase will first find all rownames that corresponds to multiple rows and then retain only the row with expression values with the highest variance. This is to ensure each row corresponds to a single gene ID (or Entrez GeneID). Then the function will normalize by rows - subtracting by the mean and dividing by the standard deviation.

}
\value{
    A big.matrix of normalized gene expression values. Each row will uniquely annotated to a single Entrez GeneID or alternative gene ID format. The big.matrix is ready for input as the database of gene expression profiles utilized by the ChIPXpress function.
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
\keyword{ cleanDatabase }
