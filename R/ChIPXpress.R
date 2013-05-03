ChIPXpress <-
function(TFID,ChIP,DB,w=0.1,c=0,warn=FALSE,DBmu=NULL,DBvar=NULL){
  if(sum(TFID %in% rownames(DB))==0) {
    stop("ERROR: TF EntrezID cannot be found in database")
  } else {
    if(warn) {
        DBcv <- sqrt(DBvar)/DBmu
        print(paste("Variance of TF:",round(DBvar[as.character(TFID)],3)))
        print(paste("Mean of TF:",round(DBmu[as.character(TFID)],3)))
        print(paste("CV of TF:",round(DBcv[as.character(TFID)],3)))

        if(DBmu[as.character(TFID)] < quantile(DBmu,0.25)) 
            print("TF mean below 25% quantile.")
        if(DBvar[as.character(TFID)] < quantile(DBvar,0.25)) 
            print("TF variance below 25% quantile.")
        if(DBcv[as.character(TFID)] < quantile(DBcv,0.25)) 
            print("TF coefficient of variation below 25% quantile.")

        if(DBmu[as.character(TFID)] < quantile(DBmu,0.05)) 
            print("TF mean below 5% quantile.")
        if(DBvar[as.character(TFID)] < quantile(DBvar,0.05)) 
            print("TF variance below 5% quantile.")
        if(DBcv[as.character(TFID)] < quantile(DBcv,0.05)) 
            print("TF coefficient of variation below 5% quantile.")
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
