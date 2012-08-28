ChIPXpress <-
function(TFID,ChIP,DB){
  if(sum(TFID %in% rownames(DB))==0) {
    stop("ERROR: TF EntrezID cannot be found in database")
  } else {
    index <- DB[as.character(TFID),] > 0
    if(sum(index) < 2) {
        stop("ERROR: Need more samples in compendium! 
Need at least 2 samples in which TF is above the mean TF expression")
    }
    TFmat <- DB[as.character(TFID),index]
    DBmat <- DB[,index]
    Out <- apply(DBmat,1,function(tmpmat) cor(tmpmat,TFmat))
    Out <- Out[names(Out) %in% ChIP]
    Out <- sort(Out,decreasing=T)
    Miss.Genes <- ChIP[!(ChIP %in% rownames(DB))]

    message("Success!")
    message(length(Out)," out of ",length(ChIP)," (",
            100*round(length(Out)/length(ChIP),3),"%) genes processed")
    
    return(list(ChIPXpress.ranking = Out,
                MissingEntrezGeneIDs = Miss.Genes))
  }
}
