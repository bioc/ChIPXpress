ChIPXpress <-
function(TFID,ChIP,DB,w=0.1,c=0,warn=FALSE,DBmu=NULL,DBvar=NULL,mu.co=0.1,var.co=0.1){
  if(sum(TFID %in% rownames(DB))==0) {
    stop("ERROR: TF EntrezID cannot be found in database")
  } else {
    if(warn) {
        if(DBvar[as.character(TFID)] < var.co) 
            warning("Variance of TF is below low variance cutoff")
        if(|DBmu[as.character(TFID)]| < mu.co)
            warning("Mean TF expression is below low expression cutoff")
    }
    index <- DB[as.character(TFID),] > c
    if(sum(index) < 2) {
        stop("ERROR: Need more samples in compendium! 
Need at least 2 samples in which TF is above TF expression cutoff")
    }
    TFmat <- DB[as.character(TFID),index]
    DBmat <- DB[,index]
    Out <- apply(DBmat,1,function(tmpmat) cor(tmpmat,TFmat))
    Out <- Out[names(Out) %in% ChIP]
    Out <- sort(Out,decreasing=T)
    Miss.Genes <- ChIP[!(ChIP %in% rownames(DB))]

    ## Weight Ranking Combination
    x <- unique(ChIP); y <- unique(names(Out))
    xc <- x[x %in% y]; yc <- y[y %in% x]
    Score <- rep(0, length(xc))
    for(i in 1:length(xc)){
        Score[i] <- w*i+(1-w)*which(yc == xc[i])
    }
    names(Score) <- xc
    Score <- sort(Score)

    message("Success!")
    message(length(Score)," out of ",length(ChIP)," (",
            100*round(length(Score)/length(ChIP),3),"%) genes processed")
    
    return(list(ChIPXpress.ranking = Score,
                MissingEntrezGeneIDs = Miss.Genes))
  }
}
