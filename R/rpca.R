
#' Perform robust PCA (Hubert PCA)
#'
#' This function performs the variant of robust PCA developed by Mia Hubert et al. on a
#' matrix of counts
#'
#' @param data Log-transformed gene expression counts, rows = samples; columns = features (genes)
#' @return A list containing the output of PcaHubert, the (default) PC scores, and the PC loadings
do.robpca <- function(data, ncp=10,...){
    pca <- rrcov::PcaHubert(data, k = ncp)

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

#' Perform robust PCA (Hubert PCA) with modified PC scores on a Seurat object
#'
#' \code{scRobustPCA} performs the variant of robust PCA developed by Mia Hubert
#' et al. on the gene expression matrix "data" in a Seurat object. Set
#' reduction.type='rpca' in other Seurat functions to use the rPCA results, for
#' example to calculate clusters with FindClusters.
#'
#' @param object Seurat object
#' @param npcs Number of principal components to calculate
#' @param pc.genes Genes as input to PC. If \code{pc.genes==NULL}, the var.genes
#'   slot is used. If var.genes is empty and \code{pc.genes==NULL}, then all
#'   genes are used.
#' @param use.modified.pcscores If \code{FALSE}, then the raw pc score output
#'   from robust PCA is used. If \code{TRUE}, then pc scores are replaced with
#'   the sum of the top 30 genes by positive loading minus the sum of the top 30
#'   genes by negative loading. Each gene is scaled to a max of 1 and min of 0.
#' @return A Seurat object with the 'rpca' field filled.
#'
#' @examples
#'
#' @export
RunRobPCA <- function(object, npcs=10, pc.genes=NULL, use.modified.pcscores=TRUE){

  print('Running rPCA')
  if(is.null(pc.genes)){
    pc.genes = object@var.genes
  } else{
    print("Variable genes not found. Using all genes in dataset.")
    pc.genes = rownames(object@data)
  }

  mat.for.pca <- t(as.matrix(object@data[pc.genes, ]))

  output <- do.robpca(mat.for.pca, ncp = npcs)

  # print(rownames(output[[1]]))
  print(colnames(output[[1]]))
  print(colnames(output[[2]]))
  # print(rownames(output[[2]]))
  if(use.modified.pcscores){
    print('Calculating modified PCs')
    pc.sigs <- CalcModifiedPCscores(object = object, loadings = output[[2]], num.sig.genes = 30)
  } else{
    pc.sigs <- as.matrix(output[[1]])
  }


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
  object <- Seurat::SetDimReduction(object, 'rpca', 'cell.embeddings', new.data = pc.sigs)

  print('Adding rPCA loadings to Seurat object')
  object <- Seurat::SetDimReduction(object, 'rpca', 'gene.loadings', new.data = as.matrix(output[[2]]))

  object <- Seurat::SetDimReduction(object, 'rpca', 'key', new.data = 'PC')

  return(object)
}

#' Plot the results of rPCA with cells colored by their identity class.
#'
#' @param object Seurat object
#' @param ... Additional parameters to DimPlot; for example, which dimensions to
#'   plot.
#' @export
RPCAPlot <- function(object,...){
  object <- Seurat::SetDimReduction(object, 'pca','cell.embeddings',
                            new.data = Seurat::GetDimReduction(object=object, reduction.type = 'rpca',slot = 'cell.embeddings'))
  Seurat::PCAPlot(object = object)
}
