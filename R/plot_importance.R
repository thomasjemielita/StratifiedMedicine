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

plot_importance <- function(object, top_n=NULL, ...) {
  
  if (object$filter=="filter_glmnet") { 
    plt.out <- plot_vimp_glmnet(object$filter.mod)
  }
  if (object$filter=="filter_ranger") {
    plt.out <- plot_vimp_ranger(object$filter.mod, top_n=top_n)
  }
  return(plt.out)
}