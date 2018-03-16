CalcModifiedPCscores <- function(object, loadings, n.genes=30){
  loadings
  ncp <- ncol(loadings) - 1
  
  pc.sigs <- as.data.frame(do.call(cbind, lapply(dims, function(dim) {
    setorderv(loadings,dim,-1)
    genes.pos <- loadings[,gene][1:num.sig.genes]
    setorderv(loadings,dim,1)
    genes.neg <- loadings[,gene][num.sig.genes:1]
    
    # generate signature value of top PC +/- correlated genes 
    # normalize every gene
    # Renames colums to PCx.pos, PCx.neg, and PCx for the xth PC
    cast.sig <- as.data.frame(cbind(rowSums(data.scaled[, genes.pos]), rowSums(data.scaled[, genes.neg])))
    cast.sig$pc.score <- cast.sig[, 1] - cast.sig[, 2]
    colnames(cast.sig) <- c(paste0(dim,'.pos'), paste0(dim,'.neg'), dim)
    
    return(as.matrix(cast.sig))
  }
  )))
}