library(devtools)
devtools::document('~/Dropbox/singlecell-pipeline/scRobustPCA/')
devtools::install('~/Dropbox/singlecell-pipeline/scRobustPCA')
# rm(list=ls())
library(scRobustPCA)
# load('~/Dropbox/singlecell-pipeline/maca/endothelial/data/MACA_HeartCap.Rdata')

tiss.cap <- FindVariableGenes(tiss.cap, do.plot = F)

npcs=5
tiss.cap <- RunRobPCA(tiss.cap, npcs = npcs)
tiss.cap <- RunTSNE(tiss.cap, reduction.use='rpca', dims.use = 1:npcs)
tiss.cap <- FindClusters(tiss.cap, reduction.type = 'rpca', dims.use = 1:npcs, force.recalc = T)
TSNEPlot(tiss.cap)
TSNEPlot(SetAllIdent(tiss.cap, id = 'subtissue'))
pairs(GetDimReduction(tiss.cap, reduction.type = 'rpca'))

npcs=5
tiss.cap <- RunPCA(tiss.cap, pcs.compute = npcs, do.print = F)
tiss.cap <- RunTSNE(tiss.cap, reduction.use='pca', dims.use = 2:npcs)
tiss.cap <- FindClusters(tiss.cap, reduction.type = 'pca', dims.use = 2:npcs, force.recalc = T)
TSNEPlot(tiss.cap)
TSNEPlot(SetAllIdent(tiss.cap, id = 'subtissue'))
pairs(GetDimReduction(tiss.cap, reduction.type = 'pca'))


sessionInfo()
