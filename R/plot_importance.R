#' Importance Plot: Visualize relative importance of variables
#' 
#' Importance is currently based on the PRISM filter model. For elastic net (filter_glmnet).
#' variables with non-zero coefficients are shown. For random forest variable importance 
#' (filter_ranger), variables are sorted by their p-values, and "top_n" will show only the
#' "top_n" most importance variables (based on p-values).
#'
#' @param object PRISM object
#' @param top_n Show top_n variables only, default=NULL (show all)
#' @param ... Additional arguments (currently ignored).
#' @return Plot (ggplot2) object
#' @export
#' @examples
#' library(StratifiedMedicine)
#' ## Continuous ##
#' dat_ctns = generate_subgrp_data(family="gaussian")
#' Y = dat_ctns$Y
#' X = dat_ctns$X
#' A = dat_ctns$A
#'
#' mod1 = filter_train(Y=Y, A=A, X=X)
#' plot_importance(mod1)

plot_importance <- function(object, top_n=NULL, ...) {
  
  if (class(object)=="filter_train") {
    filter.mod <- object$mod
  }
  if (class(object)=="PRISM") {
    filter.mod <- object$filter.mod
  }
  if (object$filter=="filter_glmnet" | object$filter=="glmnet") { 
    plt.out <- plot_vimp_glmnet(filter.mod)
  }
  if (object$filter=="filter_ranger" | object$filter=="ranger") {
    plt.out <- plot_vimp_ranger(filter.mod, top_n=top_n)
  }
  return(plt.out)
}