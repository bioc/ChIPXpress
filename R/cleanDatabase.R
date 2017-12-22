cleanDatabase <- function(DB,SaveFile="newDB.bigmemory",SavePath="."){
    ## Remove NAs
    DB <- DB[rownames(DB)!="NA",]

    ## Retain high variance probe
    MultiIDs <- table(rownames(DB))
    rowdex <- length(MultiIDs)
    MultiIDs <- names(MultiIDs)[MultiIDs > 1]
    multilen <- length(MultiIDs)

    cleanDB <- filebacked.big.matrix(nrow=rowdex,
               ncol=ncol(DB), type="double",
               backingfile=SaveFile,
               dimnames=list(1:rowdex, colnames(DB)),
               backingpath=SavePath, 
               descriptorfile=paste(SaveFile,"desc",sep="."))

    if(multilen > 0){
        index <- rownames(DB) %in% MultiIDs
        cleanDB[(multilen+1):rowdex,] <- DB[!index,]
        rownames(cleanDB)[(multilen+1):rowdex] <- rownames(DB)[!index]

        for(i in 1:length(MultiIDs)){
            DBother <- DB[rownames(DB) %in% MultiIDs[i],]
            DBvar <- apply(DBother,1,var)
            cleanDB[i,] <- DBother[which(DBvar==max(DBvar))[1],]
        }
        rownames(cleanDB)[1:multilen] <- MultiIDs
    } else {
        cleanDB[1:nrow(DB),] <- DB[1:nrow(DB),]
        rownames(cleanDB) <- rownames(DB)
    }
    ### Write normalization code
    if(rowdex <= 1000) {
        cleanDB[1:rowdex,] <- t(scale(t(cleanDB[1:rowdex,])))
    } else {
      tmpdex <- c(seq(0,rowdex,by=1000),rowdex)
      for(k in 1:(length(tmpdex)-1)){
        cleanDB[(tmpdex[k]+1):tmpdex[k+1],] <- t(scale(t(cleanDB[(tmpdex[k]+1):
                                                         tmpdex[k+1],])))
      }
    }

    return(cleanDB)
}
