# Most of this code is taken from the tutorial at https://hms-dbmi.github.io/scw/WGCNA.html
# Plotting additions and interface with Seurat are from me
#' Title
#'
#' @param seuratObj
#' @param min.module.size
#' @param filter.mito.ribo.genes
#'
#' @import dplyr
#' @importFrom Seurat FetchData WhichCells SetAllIdent SubsetData FindAllMarkers
#' @importFrom WGCNA pickSoftThreshold adjacency TOMsimilarityFromExpr labels2colors plotDendroAndColors TOMplot bicor
#' @importFrom flashClust flashClust
#' @importFrom dynamicTreeCut cutreeDynamic
#' @importFrom stats as.dist
#'
#' @return
#' @export
#'
#' @examples
#'
#' @param seuratObj
#' @param min.module.size

seuratClusterWGCNA <- function(seuratObj,
                               min.module.size = 50,
                               filter.mito.ribo.genes = FALSE) {

    options(stringsAsFactors = FALSE)
  WGCNA::allowWGCNAThreads()

  seuratObj.data <- as.matrix(seuratObj@data) %>%
    dplyr::na_if(y = 0) %>%
    t()
  names(seuratObj.data) <- colnames(seuratObj.data)
  gsg <- goodSamplesGenes(seuratObj.data,
                          minFraction = 0.25,
                          verbose = 6)
  datExpr <- seuratObj.data[gsg$goodSamples, gsg$goodGenes]

  if (isTRUE(filter.mito.ribo.genes)) {
    datExpr <- datExpr[, grep(pattern = "^MT|RP[LS]|MALAT",
                              x = colnames(datExpr),
                              value = TRUE,
                              invert = TRUE
                              )
                       ]
  }

  # Choose a set of soft-thresholding powers
  powers <- c(c(1:10), seq(from = 10, to=30, by=2))

  # Call the network topology analysis function
  sft <- pickSoftThreshold(data = datExpr,
                           powerVector = powers,
                           removeFirst = TRUE,
                           verbose = 5,
                           networkType = "signed")

  # Plot the results:
  # Scale-free topology fit index as a function of the soft-thresholding power
  # Red line corresponds to using an R^2 cut-off
  topology.fit.index <- ggplot(data = sft$fitIndices,
                               mapping = aes(x = Power,
                                             y = -sign(slope) * SFT.R.sq)) +
    geom_text(mapping = aes(label = Power,
                            color = "Red")) +
    labs(x = "Soft Threshold (power)",
         y = "Scale Free Topology Model Fit, signed R^2",
         title = "Scale independence") +
    geom_hline(aes(yintercept = 0.8, color = "Red")) +
    theme(legend.position = "none")

  # Mean connectivity as a function of the soft-thresholding power
  mean.connect <- ggplot(data = sft$fitIndices,
                         mapping = aes(x = Power,
                                       y = mean.k.)) +
    geom_text(mapping = aes(label = Power,
                            color = "Red")) +
    labs(x = "Soft Threshold (power)",
         y = "Mean Connectivity",
         title = "Mean Connectivity") +
    theme(legend.position = "none")

  if (max(sft$fitIndices$SFT.R.sq) >= 0.8) {
    sft.values <- sft$fitIndices %>%
      dplyr::filter(SFT.R.sq >= 0.8)
    softPower <- sft$fitIndices[which(sft$fitIndices$SFT.R.sq ==
                                        min(sft.values$SFT.R.sq)) - 1, ]$Power
  } else {
    softPower <- which(sft$fitIndices$SFT.R.sq == max(sft$fitIndices$SFT.R.sq))
  }


  adj <- adjacency(datExpr,
                   type = "signed",
                   power = softPower,
                   corFnc = "bicor")

  # Turn adjacency into topological overlap
  TOM <- TOMsimilarity(adjMat = adj,
                       TOMType = "signed",
                       verbose = 6,
                       TOMDenom = "mean")

  rownames(TOM) <- colnames(TOM) <- SubGeneNames <- colnames(datExpr)
  dissTOM <- 1 - TOM

  #hierarchical clustering of the genes based on the TOM dissimilarity measure
  geneTree <- flashClust(d = as.dist(dissTOM))

  # Module identification using dynamic tree cut
  minModuleSize <- min.module.size

  dynamicMods <- cutreeDynamic(dendro = geneTree,
                               method = "tree",
                               minClusterSize = minModuleSize,
                               distM = dissTOM,
                               deepSplit = 2,
                               pamRespectsDendro = FALSE
  )

  dynamicColors <- labels2colors(dynamicMods)

  pDC <- plotDendroAndColors(geneTree,
                             dynamicColors,
                             "Dynamic Tree Cut",
                             dendroLabels = FALSE,
                             hang = 0.03,
                             addGuide = TRUE,
                             guideHang = 0.05,
                             main = "Gene dendrogram and module colors"
  )

  #set the diagonal of the dissimilarity to NA
  diag(dissTOM) <- NA

  # Visualize the TOM plot.
  # Raise the dissimilarity matrix to a power
  # to bring out the module structure
  tmp <- TOMplot(dissTOM ^ 4, geneTree, as.character(dynamicColors))

  module_colors <- unique(dynamicColors)

  modules <- lapply(module_colors, function(x) {
    SubGeneNames[which(dynamicColors == x)]
  })

  names(modules) <- module_colors

  plots <- plot_grid(topology.fit.index,
                     mean.connect,
                     pDC,
                     tmp,
                     ncol = 2
  )

  return(list("plots" = plots,
              "modules" = modules,
              "adjacency.matrix" = adj,
              "exprMatrix" = datExpr
  )
  )
}
