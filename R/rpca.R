
#' Perform robust PCA (Hubert PCA)
#'
#' This function performs the variant of robust PCA developed by Mia Hubert et al. on a
#' matrix of counts
#'
#' @param data Log-transformed gene expression counts
#' @return A list containing the output of PcaHubert, the (default) PC scores, and the PC loadings
do.robpca <- function(data, ncp=10,...){
    pca <- rrcov::PcaHubert(data, k = ncp,kmax=ncp)

    scores <- as.data.frame(rrcov::getScores(pca))
    scores$cell.name <- rownames(scores)

    rpca.loadings <- as.data.frame(rrcov::getLoadings(pca))

    eigenvalues = rrcov::getEigenvalues(pca)
    list(scores, rpca.loadings, eigenvalues)
}


CalcModifiedPCscores <- function(object, loadings, num.sig.genes=30){
  ncp <- ncol(loadings)

  dims <- colnames(loadings)
  pc.sigs <- as.data.frame(do.call(cbind, lapply(dims, function(dim) {

    dim.order <- order(loadings[, dim], decreasing = T)
    genes.pos <- rownames(loadings)[dim.order[1:num.sig.genes]]
    #print(genes.pos)

    dim.order <- order(loadings[, dim], decreasing = F)
    genes.neg <- rownames(loadings)[dim.order[num.sig.genes:1]]

    # generate signature value of top PC +/- correlated genes
    # normalize every gene
    # Renames colums to PCx.pos, PCx.neg, and PCx for the xth PC
    data.pos <- t(as.matrix(object@data[genes.pos, ]))
    data.neg <- t(as.matrix(object@data[genes.neg, ]))

    data.pos <- data.pos / apply(data.pos, 2, max)
    data.neg <- data.neg / apply(data.neg, 2, max)

    cast.sig <- as.data.frame(cbind(rowSums(data.pos),
                                    rowSums(data.neg)))

    cast.sig$pc.score <- cast.sig[, 1] - cast.sig[, 2]
    colnames(cast.sig) <- c(paste0(dim,'.pos'), paste0(dim,'.neg'), dim)

    return(as.matrix(cast.sig))
  }
  )))

  pc.sigs <- pc.sigs[object@cell.names, ]
  pc.sigs <- as.matrix(pc.sigs[, grepl("PC[0-9]+$", colnames(pc.sigs))])


  return(pc.sigs)
}

#' @export
RunRobPCA <- function(object, npcs=10, pc.genes=NULL){

  print('Running rPCA')
  if(is.null(pc.genes)){
    pc.genes = object@var.genes
  } else{
    pc.genes = rownames(object@data)
  }

  mat.for.pca <- t(as.matrix(object@data[pc.genes, ]))
  output <- do.robpca(mat.for.pca, ncp = npcs)

  print('Calculating modified PCs')
  pc.sigs <- CalcModifiedPCscores(object = object, loadings = output[[2]], num.sig.genes = 30)


  reduction.name = 'rpca'

  gene.loadings = as.matrix(output[[2]])
  cell.embeddings = pc.sigs
  reduction.key = 'PC'
  sdev=numeric(0)

  pca.obj <- new(Class = "dim.reduction", gene.loadings = gene.loadings,
                 cell.embeddings = cell.embeddings, sdev = sdev, key = reduction.key)
  eval(expr = parse(text = paste0("object@dr$", reduction.name,
                                  "<- pca.obj")))
  parameters.to.store <- as.list(environment(), all = TRUE)[names(formals("RunPCA"))]
  object <- Seurat:::SetCalcParams(object = object, calculation = "RunPCA",
                          ... = parameters.to.store)
  if (is.null(object@calc.params$RunPCA$pc.genes)) {
    object@calc.params$RunPCA$pc.genes <- pc.genes
  }

#
  print('Adding modified PCs to Seurat object')
  object <- SetDimReduction(object, 'rpca', 'cell.embeddings', new.data = pc.sigs)

  print('Adding rPCA loadings to Seurat object')
  object <- SetDimReduction(object, 'rpca', 'gene.loadings', new.data = as.matrix(output[[2]]))

  object <- SetDimReduction(object, 'rpca', 'key', new.data = 'PC')

  return(object)
}

#' @export
RPCAPlot <- function(object){
  object <- SetDimReduction(object, 'pca','cell.embeddings',
                            new.data = GetDimReduction(object=object, reduction.type = 'rpca',slot = 'cell.embeddings'))
  PCAPlot(object)
}
