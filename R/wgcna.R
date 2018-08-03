# Most of this code is taken from the tutorial at https://hms-dbmi.github.io/scw/WGCNA.html
# Plotting additions and interface with Seurat are from me
#' Title
#'
#' @param seuratObj
#' @param group.by
#' @param subset.ident
#' @param cluster.use
#' @param markers
#' @param do.return
#'
#' @import dplyr
#' @importFrom Seurat FetchData WhichCells SetAllIdent SubsetData FindAllMarkers
#' @importFrom WGCNA pickSoftThreshold adjacency TOMsimilarityFromExpr labels2colors plotDendroAndColors TOMplot
#' @importFrom flashClust flashClust
#' @importFrom dynamicTreeCut cutreeDynamic
#' @importFrom stats as.dist cor
#'
#' @return
#' @export
#'
#' @examples
clusterWGCNA <- function(seuratObj,
                         group.by = NULL,
                         subset.ident = NULL,
                         cluster.use = NULL,
                         markers = NULL,
                         minModuleSize = 5,
                         clustering.method = "average",
                         do.return = TRUE) {
  options(stringsAsFactors = FALSE)
  allowWGCNAThreads()
  enableWGCNAThreads()

  if (!is.null(subset.ident) & !is.null(cluster.use)) {
    obj.mat <- FetchData(
      seuratObj,
      cells.use = WhichCells(
        seuratObj,
        subset.name = subset.ident,
        accept.value = cluster.use
      ),
      vars.all = rownames(seuratObj@data)
    ) %>%
      as.data.frame()
  } else {
    obj.mat <- FetchData(seuratObj,
                         vars.all = rownames(seuratObj@data))
  }

  if (!is.null(group.by)) {
    seuratObj <- SetAllIdent(
      object = SubsetData(
        seuratObj,
        subset.name = subset.ident,
        accept.value = cluster.use
      ),
      id = group.by
    )
  }

  if (is.null(markers)) {
    markers <- FindAllMarkers(seuratObj, min.pct = 0.3)
  }

  obj.mat <- obj.mat[, unique(markers$gene)]
  if (dim(obj.mat)[[2]] == 0) {
    stop("There is no data to analyze.")
  }

  powers <- c(1:10,
              seq(from = 12,
                  to = 40,
                  by = 2))

  sft <- pickSoftThreshold(
    obj.mat,
    dataIsExpr = TRUE,
    powerVector = powers,
    corFnc = cor,
    corOptions = list(use = 'p'),
    networkType = "signed"
  )

  # Plot the results
  # Scale-free topology fit index as a function of the soft-thresholding power
  # Red line corresponds to using an R^2 cut-off
  #   topology.fit.index <- ggplot(data = sft$fitIndices,
  #                                mapping = aes(x = Power,
  #                                              y = -sign(slope)*SFT.R.sq)
  #   ) +
  #     geom_text(mapping = aes(label = Power,
  #                             color = 'Red')) +
  #     labs(x = "Soft Threshold (power)",
  #          y = "Scale Free Topology Model Fit, signed R^2",
  #          title = "Scale independence") +
  #     geom_hline(aes(yintercept = 0.8, color = 'Red')) +
  #     theme(legend.position = "none")

  #   # Mean connectivity as a function of the soft-thresholding power
  #   mean.connect <- ggplot(data = sft$fitIndices,
  #                          mapping = aes(x = Power,
  #                                        y = mean.k.)
  #   ) +
  #     geom_text(mapping = aes(label = Power,
  #                             color = 'Red')) +
  #     labs(x = "Soft Threshold (power)",
  #          y = "Mean Connectivity",
  #          title = "Mean Connectivity") +
  #     theme(legend.position = "none")
  #   plot_grid(topology.fit.index, mean.connect, ncol = 2)

  # Set the soft-thresholding power
  if (max(sft$fitIndices$SFT.R.sq) >= 0.8) {
    sft.values <- sft$fitIndices %>% dplyr::filter(SFT.R.sq >= 0.8)
    softPower <-
      sft$fitIndices[which(sft$fitIndices$SFT.R.sq == min(sft.values$SFT.R.sq)) - 1, ]$Power
  } else {
    softPower <-
      which(sft$fitIndices$SFT.R.sq == max(sft$fitIndices$SFT.R.sq))
  }

  #calclute the adjacency matrix
  adj <- adjacency(obj.mat,
                   type = "signed",
                   power = softPower)

  #turn adjacency matrix into a topological overlap matrix (TOM) to minimize the effects of noise and spurious associations
  TOM <- TOMsimilarityFromExpr(obj.mat,
                               networkType = "signed",
                               TOMType = "signed",
                               power = softPower)

  rownames(TOM) <-
    colnames(TOM) <- SubGeneNames <- colnames(obj.mat)
  dissTOM = 1 - TOM

  #hierarchical clustering of the genes based on the TOM dissimilarity measure
  geneTree = flashClust(as.dist(dissTOM), method = clustering.method)

  # Module identification using dynamic tree cut

  dynamicMods <- cutreeDynamic(dendro = geneTree,
                               method = "tree",
                               minClusterSize = minModuleSize)

  dynamicColors <- labels2colors(dynamicMods)

  plotDendroAndColors(
    geneTree,
    dynamicColors,
    "Dynamic Tree Cut",
    dendroLabels = FALSE,
    hang = 0.03,
    addGuide = TRUE,
    guideHang = 0.05,
    main = "Gene dendrogram and module colors"
  )

  #set the diagonal of the dissimilarity to NA
  diag(dissTOM) = NA

  #Visualize the Tom plot. Raise the dissimilarity matrix to a power  to bring out the module structure
  TOMplot(dissTOM ^ 4, geneTree, as.character(dynamicColors))

  module_colors <- setdiff(unique(dynamicColors), "grey")
  modules <-
    lapply(module_colors, function(x) {
      SubGeneNames[which(dynamicColors == x)]
    })
  names(modules) <- module_colors

  # module.order <- unlist(tapply(1:ncol(obj.mat),as.factor(dynamicColors),I))
  # m < -t(t(obj.mat[,module.order])/apply(obj.mat[,module.order],2,max))
  #heatmaply(t(m),zlim=c(0,1),col=gray.colors(100),Rowv=NA,Colv=NA,scale="none",RowSideColors=dynamicColors[module.order])

  if (do.return) {
    return(modules)
  }
}

#' Title
#'
#' @param seuratObj
#' @param modules
#' @param ...
#'
#' @import ggplot2
#' @import magrittr
#' @importFrom seuratBubblePlot bubbleplot
#' @importFrom gridExtra grid.arrange
#'
#' @return
#' @export
#'
#' @examples
plotModules <- function(seuratObj, modules, ...) {
  modplots = lapply((1:length(modules)),
                    function(x) {
                      bp <- bubbleplot(
                        seuratObj = seuratObj,
                        genes.plot = modules[[x]],
                        do.return = TRUE,
                        ...
                      ) +
                        labs(title = names(modules)[[x]]) +
                        xlab(NULL) +
                        ylab(NULL) +
                        theme(legend.position = "none")
                    })

  grid.arrange(grobs = modplots)
}
