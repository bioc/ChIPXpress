buildDatabase <- function(GPL_id,GSMfiles=NULL,SaveDir=NULL,LoadPrevious=FALSE){
    ## Obtain sample names
    if(is.null(GSMfiles)){
        GSMfilenames <- Meta(getGEO(GPL_id))$sample_id
    } else {
        GSMfilenames <- GSMfiles
    }

    ## Download and Process
    if(is.null(SaveDir)){
        stop("ERROR: Please enter a directory to save output.")
    } else {
        j=1
        options(bigmemory.allow.dimnames=TRUE)
        ## Check if previous file exists
        if(file.exists(paste(SaveDir,"tmpmatrixGEO.bigmemory.desc",sep="/")) &&
           file.exists(paste(SaveDir,"tmpmatrixGEO.bigmemory",sep="/")) &&
           LoadPrevious){
            
            ## Remove old .gz files
            files <- list.files()
            files <- files[grep(".gz",files)]
            system(paste("rm",paste(files,collapse=" ")))

            message("Loading Previous File...")
            bigdb <- attach.big.matrix("tmpmatrixGEO.bigmemory.desc",path=SaveDir)
            colcnts <- colsum(bigdb)
            j=which(colcnts==0)[1]
        }

        for(i in j:length(GSMfilenames)){
            message("Downloading File ",i,": ",GSMfilenames[i],"...")
            getGEOSuppFiles(GSMfilenames[i],makeDirectory=FALSE,baseDir=SaveDir)

            message("Processing File ",i,": ",GSMfilenames[i],"...")
            object <- ReadAffy(celfile.path=SaveDir) ## Read CEL files
            object <- exprs(frma(object)) ## Process with frma

            ## Store experssion Output
            if(i==1)   bigdb <- filebacked.big.matrix(nrow=nrow(object),
               ncol=length(GSMfilenames), type="double",
               backingfile="tmpmatrixGEO.bigmemory",
               dimnames=list(rownames(object), GSMfilenames),
               backingpath=SaveDir, descriptorfile="tmpmatrixGEO.bigmemory.desc")
            bigdb[,i] <- object[,i]

            ## Remove downloaded files
            files <- list.files(SaveDir)
            files <- files[grep(".gz",files)]
            files <- paste(SaveDir,files,sep="/")
            system(paste("rm",paste(files,collapse=" ")))
        }
        return(bigdb)
    }
}

