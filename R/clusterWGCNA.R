# Most of this code is taken from the tutorial at https://hms-dbmi.github.io/scw/WGCNA.html
# Plotting additions and interface with Seurat are from me
#' Title
#'
#' @param seuratObj
#' @param min.module.size
#' @param filter.mito.ribo.genes
#' @param assay.use
#' @param slot.use
#' @param merge.similar.modules
#' @param merge.similarity.threshold
#'
#' @import dplyr
#' @importFrom Seurat FetchData WhichCells SetAllIdent SubsetData FindAllMarkers
#' @importFrom WGCNA pickSoftThreshold adjacency TOMsimilarityFromExpr labels2colors plotDendroAndColors TOMplot bicor
#' @importFrom flashClust flashClust
#' @importFrom dynamicTreeCut cutreeDynamic
#' @importFrom stats as.dist
#' @importFrom glue glue
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
                               minFraction = 0.25,
                               filter.mito.ribo.genes = FALSE,
                               assay.use = NULL,
                               slot.use = 'data',
                               merge.similar.modules = FALSE,
                               merge.similarity.threshold = 0.25) {

  options(stringsAsFactors = FALSE)
  WGCNA::allowWGCNAThreads()

  if(!is.null(assay.use)){
    seuratObj.data <- GetAssayData(seuratObj,
                                   assay.type = assay.use,
                                   slot = slot.use) %>%
      dplyr::na_if(y = 0) %>%
      t()
  } else {
    seuratObj.data <- as.matrix(seuratObj@data) %>%
      dplyr::na_if(y = 0) %>%
      t()
  }

  names(seuratObj.data) <- colnames(seuratObj.data)
  gsg <- goodSamplesGenes(seuratObj.data,
                          minFraction = minFraction,
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

  print(glue("Using {softPower} for softPower"))

  print(glue("Calculating adjacency matrix..."))
  adj <- adjacency(datExpr,
                   type = "signed",
                   power = softPower,
                   corFnc = "bicor")

  print(glue("Turning adjacency into topological overlap..."))
  TOM <- TOMsimilarity(adjMat = adj,
                       TOMType = "signed",
                       verbose = 6,
                       TOMDenom = "mean")

  rownames(TOM) <- colnames(TOM) <- SubGeneNames <- colnames(datExpr)
  dissTOM <- 1 - TOM

  print(glue("Hierarchical clustering of the genes based on the TOM dissimilarity measure..."))
  geneTree <- flashClust(d = as.dist(dissTOM))

  print(glue("Module identification using dynamic tree cut..."))
  dynamicMods <- cutreeDynamic(dendro = geneTree,
                               method = "hybrid",
                               minClusterSize = min.module.size,
                               distM = dissTOM,
                               deepSplit = 3,
                               pamRespectsDendro = FALSE
  )

  dynamicColors <- labels2colors(dynamicMods)

  # pdc <- plotDendroAndColors(geneTree,
  #                            dynamicColors,
  #                            "Dynamic Tree Cut",
  #                            dendroLabels = FALSE,
  #                            hang = 0.03,
  #                            addGuide = TRUE,
  #                            guideHang = 0.05,
  #                            main = "Gene dendrogram and module colors"
  #                            )
  
  #set the diagonal of the dissimilarity to NA
  # diag(dissTOM) <- NA

  # Visualize the TOM plot.
  # Raise the dissimilarity matrix to a power
  # to bring out the module structure
  # tmp <- TOMplot(dissTOM ^ 4, geneTree, as.character(dynamicColors))
  
  print(glue("Assembling modules"))
  module_colors <- unique(dynamicColors)

  modules <- lapply(module_colors, function(x) {
    SubGeneNames[which(dynamicColors == x)]
  })

  names(modules) <- module_colors

  print(glue("Calculating module eigengenes..."))
  MEList <- moduleEigengenes(expr = datExpr,
                             colors = dynamicColors,
                             softPower = softPower)
  MEs <- MEList$eigengenes

  if (isTRUE(merge.similar.modules)){
    print(glue("Prior to merging, found {length(modules)} modules."))
    print(glue("Calculate dissimilarity of module eigengenes"))
    MEDiss <- 1 - WGCNA::bicor(MEs,
                               use = 'pairwise.complete.obs')
    # Cluster module eigengenes
    METree <- hclust(as.dist(MEDiss),
                    method = "complete")

    print(glue("Merging similar modules..."))
    merge <- mergeCloseModules(exprData = datExpr,
                               colors = dynamicColors,
                               corFnc = "bicor",
                               cutHeight = merge.similarity.threshold,,
                               verbose = 3)

    moduleColors = merge$colors
    colorOrder = c("grey", standardColors(50));
    moduleLabels = match(moduleColors, colorOrder)-1;
    MEs <- merge$newMEs

    mergedModules <- lapply(unique(moduleColors), function(x){
      colnames(datExpr)[which(moduleColors == x)]
    })

    names(mergedModules) <- unique(moduleColors)

    mergedModules <- lapply(names(mergedModules), function(x) {
      stringr::str_remove(string = mergedModules[[x]], pattern = "DCA_")
    })

    names(mergedModules) <- unique(moduleColors)
    modules <- mergedModules
  } else {
    print(glue("Found {length(modules)} modules."))
  }

  return(list("modules" = modules,
              "moduleEigengenes" = MEs,
              "adjacency.matrix" = adj,
              "exprMatrix" = datExpr
              )
         )
}
