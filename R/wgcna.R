library(WGCNA)
library(flashClust)
library(Seurat)
library(tidyverse)

cluster.wcgna <- function(seuratObj, group.by, subset.ident, cluster, markers = NULL){
  options(stringsAsFactors = FALSE)
  allowWGCNAThreads()
  enableWGCNAThreads()

  obj.mat <- FetchData(seuratObj, 
                       cells.use = WhichCells(seuratObj, 
                                              subset.name = subset.ident, 
                                              accept.value = cluster), 
                       vars.all = rownames(seuratObj@data)
                       ) %>% 
    as.data.frame()
  
  seuratObj <- SetAllIdent(object = SubsetData(seuratObj, 
                                               subset.name = subset.ident, 
                                               accept.value = cluster), 
                           id = group.by)
  obj.markers <- FindAllMarkers(seuratObj, min.pct = 0.3)

  obj.mat <- obj.mat[,unique(obj.markers$gene)]
  
  powers <- c(c(1:10), seq(from = 12, to=40, by=2))
  sft <- pickSoftThreshold(obj.mat,
                           dataIsExpr = TRUE,
                           powerVector = powers,
                           corFnc = cor,
                           corOptions = list(use = 'p'),
                           networkType = "signed")
  
  # Plot the results
  
  # Scale-free topology fit index as a function of the soft-thresholding power
  # Red line corresponds to using an R^2 cut-off
  topology.fit.index <- ggplot(data = sft$fitIndices, 
                               mapping = aes(x = Power, 
                                             y = -sign(slope)*SFT.R.sq)
                               ) + 
    geom_text(mapping = aes(label = Power, 
                            color = 'Red')) +
    labs(x = "Soft Threshold (power)",
         y = "Scale Free Topology Model Fit, signed R^2",
         title = "Scale independence") + 
    geom_hline(aes(yintercept = 0.8, color = 'Red')) +
    theme(legend.position = "none")
  
  # Mean connectivity as a function of the soft-thresholding power
  mean.connect <- ggplot(data = sft$fitIndices, 
                         mapping = aes(x = Power, 
                                       y = mean.k.)
                         ) + 
    geom_text(mapping = aes(label = Power, 
                            color = 'Red')) +
    labs(x = "Soft Threshold (power)",
         y = "Mean Connectivity",
         title = "Mean Connectivity") + 
    theme(legend.position = "none")
  plot_grid(topology.fit.index, mean.connect, ncol = 2)

  # Set the soft-thresholding power
  if (max(sft$fitIndices$SFT.R.sq) >= 0.8){
    sft.values <- sft$fitIndices %>% filter(SFT.R.sq >= 0.8)
    softPower <- sft$fitIndices[which(sft$fitIndices$SFT.R.sq == min(sft.values$SFT.R.sq)) - 1,]$Power
  } else {
    softPower <- which(sft$fitIndices$SFT.R.sq == max(sft$fitIndices$SFT.R.sq))  
  }
  
  #calclute the adjacency matrix
  adj <- adjacency(obj.mat, type="signed", power=softPower)
  
  #turn adjacency matrix into a topological overlap matrix (TOM) to minimize the effects of noise and spurious associations
  TOM <- TOMsimilarityFromExpr(obj.mat,
                               networkType = "signed", 
                               TOMType = "signed", 
                               power = softPower)
  
  rownames(TOM) <- colnames(TOM) <- SubGeneNames <- colnames(obj.mat)
  dissTOM=1-TOM
  
  #hierarchical clustering of the genes based on the TOM dissimilarity measure
  geneTree = flashClust(as.dist(dissTOM), method="average")
  
  #plot the resulting clustering tree (dendrogram)
  plot(geneTree, xlab="", sub="", cex=0.3)
  
  # Set the minimum module size
  minModuleSize <- 5
  
  # Module identification using dynamic tree cut
  
  dynamicMods <- cutreeDynamic(dendro = geneTree,  method="tree", minClusterSize = minModuleSize)
  
  #the following command gives the module labels and the size of each module.
  dynamicColors = labels2colors(dynamicMods)
  
  #Plot the module assignment under the dendrogram note. 
  # The grey color is reserved for unassigned genes
  table(dynamicColors)
  dynamicColors <- labels2colors(dynamicMods)
  
  plotDendroAndColors(geneTree, 
                      dynamicColors, 
                      "Dynamic Tree Cut", 
                      dendroLabels = FALSE, 
                      hang = 0.03, 
                      addGuide = TRUE, 
                      guideHang = 0.05, 
                      main = "Gene dendrogram and module colors")
  
  #set the diagonal of the dissimilarity to NA 
  diag(dissTOM) = NA
  
  #Visualize the Tom plot. Raise the dissimilarity matrix to a power  to bring out the module structure
  TOMplot(dissTOM^4, geneTree, as.character(dynamicColors))
  
  module_colors <- setdiff(unique(dynamicColors), "grey")
  modules <- lapply(module_colors, function(x){SubGeneNames[which(dynamicColors==x)]})
  names(modules) <- module_colors
  
  # module.order <- unlist(tapply(1:ncol(obj.mat),as.factor(dynamicColors),I))
  # m < -t(t(obj.mat[,module.order])/apply(obj.mat[,module.order],2,max))
  #heatmaply(t(m),zlim=c(0,1),col=gray.colors(100),Rowv=NA,Colv=NA,scale="none",RowSideColors=dynamicColors[module.order]) 
  return(modules)
}
