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

  purrr::reduce(modplots, `+`)
}
