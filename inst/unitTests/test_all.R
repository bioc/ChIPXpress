test_all <- function()
{
    GSM_ids <- c("GSM24056","GSM24058","GSM24060",
                 "GSM94856","GSM94857","GSM94858")
    t <- tempdir()
    DB <- buildDatabase(GSMfiles=GSM_ids,SaveDir=t)
    library(mouse4302.db)
    EntrezID <- mget(as.character(rownames(DB)),mouse4302ENTREZID)
    rownames(DB) <- as.character(EntrezID)
    cleanDB <- cleanDatabase(DB,SaveFile="newDB_GPL1261.bigmemory",
                             SavePath=t)
    data(Oct4ESC_ChIPgenes)
    out <- ChIPXpress(TFID="18999",ChIP=Oct4ESC_ChIPgenes,DB=cleanDB)
    checkEquals(class(out), "list")
    checkEquals(length(out), 2)
    
}
