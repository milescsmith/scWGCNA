# This code is adapted from the tutorial at
# https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/
# Plotting additions and interface with Seurat are from me
#' @title scWGCNA
#'
#' @description Perform Weighted Gene Co-expression Network Analysis on scRNAseq
#' data in a Seurat object
#'
#' @param object Processed scRNAseq object
#' @param min.module.size Minimum number of genes needed to form an expression
#' module. Default: 50.
#' @param filter_mito_ribo_genes Should mitochondrial and ribosomal genes be
#' removed? Default: FALSE
#' @param minFraction When determining genes to keep, what is the minimum
#' fraction of samples a gene must be detected in before being remove?
#' Default: 0.25
#' @param assay_use Instead of using expression data from the obj@data slot,
#' use data stored in the obj@assay slot
#' @param slot_use Assay data slot to use (i.e. "raw.data", "data", or
#' "scale.data").  Default: 'data'
#' @param merge_similar_modules Should similar modules be merged? Default: FALSE
#' @param merge_similarity_threshold Dissimilarity (i.e., 1-correlation) cutoff
#' used to determine if modules should be merged. Default: 0.25
#'
#' @importFrom dplyr filter na_if
#' @importFrom ggplot2 ggplot geom_text geom_hline theme labs geom_text aes
#' @importFrom cowplot plot_grid
#' @importFrom WGCNA pickSoftThreshold adjacency TOMsimilarityFromExpr TOMsimilarity
#' labels2colors plotDendroAndColors TOMplot bicor allowWGCNAThreads goodSamplesGenes
#' mergeCloseModules moduleEigengenes
#' @importFrom flashClust flashClust hclust
#' @importFrom dynamicTreeCut cutreeDynamic
#' @importFrom stats as.dist
#' @importFrom glue glue
#' @importFrom stringr str_remove
#'
#' @return
#' @export
#'
#' @examples
#'
#' @param object
#' @param min.module.size
scWGCNA <- function(object,
                    min.module.size = 50,
                    minFraction = 0.25,
                    filter_mito_ribo_genes = FALSE,
                    assay_use = NULL,
                    slot_use = NULL,
                    merge_similar_modules = FALSE,
                    merge_similarity_threshold = 0.25) {

  object_data <- GatherData(object = object,
                         assay = assay_use,
                         slot_use = slot_use) %>%
    na_if(y = 0) %>%
    t()

  names(object_data) <- colnames(object_data)
  gsg <- goodSamplesGenes(object_data,
                          minFraction = minFraction,
                          verbose = 6
  )
  datExpr <- object_data[gsg$goodSamples, gsg$goodGenes]

  if (isTRUE(filter_mito_ribo_genes)) {
    datExpr <- datExpr[, grep(
      pattern = "^MT|RP[LS]|MALAT",
      x = colnames(datExpr),
      value = TRUE,
      invert = TRUE
    )]
  }

  # Choose a set of soft-thresholding powers
  powers <- c(c(1:10), seq(from = 10, to = 30, by = 2))

  # Call the network topology analysis function
  sft <- pickSoftThreshold(
    data = datExpr,
    powerVector = powers,
    removeFirst = TRUE,
    verbose = 5,
    networkType = "signed"
  )

  # Plot the results:
  # Scale-free topology fit index as a function of the soft-thresholding power
  # Red line corresponds to using an R^2 cut-off
  topology.fit.index <- ggplot(
    data = sft$fitIndices,
    mapping = aes(
      x = Power,
      y = -sign(slope) * SFT.R.sq
    )
  ) +
    geom_text(mapping = aes(
      label = Power,
      color = "Red"
    )) +
    labs(
      x = "Soft Threshold (power)",
      y = "Scale Free Topology Model Fit, signed R^2",
      title = "Scale independence"
    ) +
    geom_hline(aes(yintercept = 0.8, color = "Red")) +
    theme(legend.position = "none")

  # Mean connectivity as a function of the soft-thresholding power
  mean.connect <- ggplot(
    data = sft$fitIndices,
    mapping = aes(
      x = Power,
      y = mean.k.
    )
  ) +
    geom_text(mapping = aes(
      label = Power,
      color = "Red"
    )) +
    labs(
      x = "Soft Threshold (power)",
      y = "Mean Connectivity",
      title = "Mean Connectivity"
    ) +
    theme(legend.position = "none")

  plot_grid(topology.fit.index, mean.connect)

  if (max(sft$fitIndices$SFT.R.sq) >= 0.8) {
    sft.values <- sft$fitIndices %>%
      filter(SFT.R.sq >= 0.8)
    softPower <- sft$fitIndices[which(sft$fitIndices$SFT.R.sq ==
                                        min(sft.values$SFT.R.sq)) - 1, ]$Power
  } else {
    softPower <- which(sft$fitIndices$SFT.R.sq == max(sft$fitIndices$SFT.R.sq))
  }

  message(glue("Using {softPower} for softPower"))

  message(glue("Calculating adjacency matrix..."))
  adj <- adjacency(datExpr,
                   type = "signed",
                   power = softPower,
                   corFnc = "bicor"
  )

  message(glue("Turning adjacency into topological overlap..."))
  TOM <- TOMsimilarity(
    adjMat = adj,
    TOMType = "signed",
    verbose = 6,
    TOMDenom = "mean"
  )

  rownames(TOM) <- colnames(TOM) <- SubGeneNames <- colnames(datExpr)
  dissTOM <- 1 - TOM

  message(glue("Hierarchical clustering of the genes based on the
             TOM dissimilarity measure..."))
  geneTree <- flashClust(d = as.dist(dissTOM))

  message(glue("Module identification using dynamic tree cut..."))
  dynamicMods <- cutreeDynamic(
    dendro = geneTree,
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

  # set the diagonal of the dissimilarity to NA
  # diag(dissTOM) <- NA

  # Visualize the TOM plot.
  # Raise the dissimilarity matrix to a power
  # to bring out the module structure
  # tmp <- TOMplot(dissTOM ^ 4, geneTree, as.character(dynamicColors))

  message(glue("Assembling modules"))
  module_colors <- unique(dynamicColors)

  modules <- lapply(module_colors, function(x) {
    SubGeneNames[which(dynamicColors == x)]
  })

  names(modules) <- module_colors

  message(glue("Calculating module eigengenes..."))
  MEList <- moduleEigengenes(
    expr = datExpr,
    colors = dynamicColors,
    softPower = softPower
  )
  MEs <- MEList$eigengenes

  if (isTRUE(merge_similar_modules)) {
    print(glue("Prior to merging, found {length(modules)} modules."))
    print(glue("Calculate dissimilarity of module eigengenes"))
    MEDiss <- 1 - bicor(MEs,
                        use = "pairwise.complete.obs"
    )
    # Cluster module eigengenes
    METree <- hclust(as.dist(MEDiss),
                     method = "complete"
    )

    message(glue("Merging similar modules..."))
    merge <- mergeCloseModules(
      exprData = datExpr,
      colors = dynamicColors,
      corFnc = "bicor",
      cutHeight = merge_similarity_threshold,
      verbose = 3
    )

    moduleColors <- merge$colors
    MEs <- merge$newMEs

    mergedModules <- lapply(unique(moduleColors), function(x) {
      colnames(datExpr)[which(moduleColors == x)]
    })

    names(mergedModules) <- unique(moduleColors)

    mergedModules <- lapply(names(mergedModules), function(x) {
      str_remove(string = mergedModules[[x]], pattern = "DCA_")
    })

    names(mergedModules) <- unique(moduleColors)
    modules <- mergedModules
  } else {
    message(glue("Found {length(modules)} modules."))
  }

  return(list(
    "modules" = modules,
    "moduleEigengenes" = MEs,
    "adjacency.matrix" = adj,
    "exprMatrix" = datExpr
  ))
}
