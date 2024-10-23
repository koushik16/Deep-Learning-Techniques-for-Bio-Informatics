library(devtools)

devtools::install_github("aertslab/AUCell")
devtools::install_github("aertslab/RcisTarget")

devtools::install_github("aertslab/cisTopic")


source("https://bioconductor.org/biocLite.R")
biocLite(c('Rsubread', 'umap', 'Rtsne', 'ComplexHeatmap', 'fastcluster', 'data.table', 'rGREAT', 'ChIPseeker', 'TxDb.Hsapiens.UCSC.hg19.knownGene', 'org.Hs.eg.db'))

cisTopicObject <- readmtx('matrix')

pathToBams <- 'dataz/'
bamFiles <- paste(pathToBams, list.files(pathToBams), sep='')
cisTopicObject <- createcisTopicObjectFromBAM(bamFiles, regions, project.name='hw')

# If you want to rename cells
cell.names <- cisTopicObject@cell.names
new.cell.names <- sapply(strsplit(cell.names, split = ".", fixed=TRUE), "[", 3)
cisTopicObject <- renameCells(cisTopicObject, new.cell.names)

data(counts_mel) 
cisTopicObject <- createcisTopicObject(counts_mel, project.name='hw')
rm(counts_mel)

data(cellData_mel)
cisTopicObject <- addCellMetadata(cisTopicObject, cell.data = cellData_mel)
rm(cellData_mel)

cisTopicObject <- runWarpLDAModels(cisTopicObject, topic=c(2:15, 20, 25, 35, 40, 45, 50), seed=123, nCores=17, addModels=FALSE)

# For Warp LDA
cisTopicObject@calc.params[['runWarpLDAModels']]$seed <- seed
cisTopicObject@calc.params[['runWarpLDAModels']]$iterations <- iterations
cisTopicObject@calc.params[['runWarpLDAModels']]$alpha <- alpha
cisTopicObject@calc.params[['runWarpLDAModels']]$alphaByTopic <- alphaByTopic
cisTopicObject@calc.params[['runWarpLDAModels']]$beta <- beta

# For CGS
cisTopicObject@calc.params[['runCGSModels']]$seed <- seed
cisTopicObject@calc.params[['runCGSModels']]$burnin <- burnin
cisTopicObject@calc.params[['runCGSModels']]$iterations <- iterations
cisTopicObject@calc.params[['runCGSModels']]$alpha <- alpha
cisTopicObject@calc.params[['runCGSModels']]$alphaByTopic <- alphaByTopic
cisTopicObject@calc.params[['runCGSModels']]$beta <- beta

par(mfrow=c(3,3))
cisTopicObject <- selectModel(cisTopicObject, type='maximum')
cisTopicObject <- selectModel(cisTopicObject, type='perplexity')
cisTopicObject <- selectModel(cisTopicObject, type='derivative')

cisTopicObject <- runUmap(cisTopicObject, target='cell')

par(mfrow=c(1,2))
plotFeatures(cisTopicObject, method='Umap', target='cell', topic_contr=NULL, colorBy=c('cellLine', 'LineType'), cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE, col.low='darkgreen', col.mid='yellow', col.high='brown1', intervals=20)

par(mfrow=c(2,5))
plotFeatures(cisTopicObject, method='Umap', target='cell', topic_contr='Probability', colorBy=NULL, cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE)

cellTopicHeatmap(cisTopicObject, method='Probability', colorBy=c('LineType', 'cellLine'))

pred.matrix <- predictiveDistribution(cisTopicObject)

# Obtain signatures
path_to_signatures <- '/data/peaks'
ChIP_Seq_signatures <- paste(path_to_signatures, list.files(path_to_signatures), sep='')
labels  <- c('MITF', 'SOX10', 'TFAP2A')
cisTopicObject <- getSignaturesRegions(cisTopicObject, ChIP_Seq_signatures, labels=labels, minOverlap = 0.4)

# Compute cell rankings
library(AUCell)
aucellRankings <- AUCell_buildRankings(pred.matrix, plot=FALSE, verbose=FALSE)

# Check signature enrichment in cells
cisTopicObject <- signatureCellEnrichment(cisTopicObject, aucellRankings, selected.signatures='all', aucMaxRank = 0.1*nrow(aucellRankings), plot=FALSE)

# Plot
par(mfrow=c(2,2))
plotFeatures(cisTopicObject, method='Umap', target='cell', topic_contr=NULL, colorBy=c('LineType', 'MITF', 'SOX10', 'TFAP2A'), cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE, intervals=20)

cisTopicObject <- getRegionsScores(cisTopicObject, method='NormTop', scale=TRUE)
