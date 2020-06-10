#' Patient-Level Estimate Plot (plot_ple): Visualize distribution of estimates
#' 

#' Plots based on Patient-level estimate (see \code{ple_train}) model results. Options 
#' include "waterfall" and "density". Target controls which column of "mu_train" 
#' (from ple_train object) is shown on the plot.
#'
#' @param object \code{ple_train} object
#' @param target Which patient-level estimate to visualize. Default=NULL, which uses 
#' the estimated treatment difference. 
#' @param type TYpe of plot. Default="waterfall"; type="density" shows density plot.
#' @param ... Additional arguments (currently ignored).
#' @return Plot (ggplot2) object
#' @export
#' @examples
#' \donttest{
#' library(StratifiedMedicine)
#' ## Continuous ##
#' dat_ctns = generate_subgrp_data(family="gaussian")
#' Y = dat_ctns$Y
#' X = dat_ctns$X
#' A = dat_ctns$A
#'
#'
#' mod1 = ple_train(Y=Y, A=A, X=X, Xtest=X, ple="ranger", meta="X-learner")
#' plot_ple(mod1)
#' }
#' 
### PLE Plots (waterfall, density) ###
plot_ple <- function(object, target=NULL, type="waterfall", ...) {
  
  if (is.null(object$mu_train)){
    stop("Check ple model fit, no training estimates available.")
  }
  mu_hat <- object$mu_train
  if (is.null(mu_hat$Subgrps)) {
    mu_hat$Subgrps <- rep("Overall", dim(mu_hat)[1])
  }
  if (is.null(target)) {
    ple_name <- colnames(mu_hat)[grepl("diff", colnames(mu_hat))]
    ple_name <- ple_name[1]
    mu_hat$PLE <- mu_hat[[ple_name]]
  }
  if (!is.null(target)) {
    ple_name <- target
    mu_hat$PLE <- mu_hat[[target]]
  }
  if (!is.null(mu_hat$estimand)) {
    ## TO DO ##
  }
  ple.label <- ple_name ## TO DO ##
  y.label <- paste("Estimates:", ple.label)
  x.label <- y.label
  mu_hat$id = 1:nrow(mu_hat)
  if (type=="waterfall") {
    res = ggplot2::ggplot(mu_hat, ggplot2::aes(x=reorder(id, PLE), y=PLE, fill=Subgrps)) +
      ggplot2::geom_bar(stat="identity") +
      ggplot2::ggtitle( paste("Waterfall Plot: Patient-level Estimates,", ple.label)) +
      ggplot2::ylab(y.label) + 
      ggplot2::xlab("")+
      ggplot2::theme_bw() +
      ggplot2::theme(axis.line.x = ggplot2::element_blank(), 
                     axis.text.x = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank(),
                     axis.title.y = ggplot2::element_text(face="bold",angle=90))
  }
  if (type=="density") {
    x.label <- y.label
    res = ggplot2::ggplot(mu_hat, ggplot2::aes(PLE, fill=Subgrps)) + 
      ggplot2::geom_density(alpha=0.30) +
      ggplot2::xlab(x.label) +
      ggplot2::ggtitle(paste("Density Plot: Patient-Level Estimates,", ple.label))+
      ggplot2::theme(plot.title=ggplot2::element_text(size=16,face="bold"),
                     axis.text.y=ggplot2::element_blank(),
                     axis.ticks.y=ggplot2::element_blank(),
                     axis.text.x=ggplot2::element_text(face="bold"),
                     axis.title=ggplot2::element_text(size=12,face="bold"))+
      ggplot2::theme_bw()
  }
  return(res)
}