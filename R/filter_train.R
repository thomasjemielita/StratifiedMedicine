#' filter_train: Identify variables of interest
#'
#' Wrapper function to train a filter model to determine variables associated with 
#' the outcome and/or treatment.. Options include elastic net (glmnet) and random 
#' forest based variable importance (ranger). Used directly in PRISM. 
#'
#' @inheritParams PRISM
#' @param hyper Hyper-parameters for the filter model (must be list). Default is NULL.
#' See details below. 
#' @param ... Any additional parameters, not currently passed through.
#'
#'
#' @return Trained filter model and vector of variable names that pass the filter. 
#'
#'  \itemize{
#'   \item mod - trained model
#'   \item filter.vars - Variables that remain after filtering (could be all)
#' }
#' 
#' @details filter_train currently fits elastic net or random forest to
#' find a reduced set of variables which are likely associated with the outcome (Y) 
#' and/or treatment (A). Current options include:
#' 
#' 1. \strong{glmnet}: Wrapper function for the function "glmnet" from the glmnet package. Here, 
#' variables with estimated elastic net coefficients of 0 are filtered. Uses LM/GLM/cox 
#' elastic net for family="gaussian","binomial", "survival" respectively. Default is to 
#' regress Y~ENET(X) with hyper-parameters:
#' 
#' hyper = list(lambda="lambda.min", family="gaussian",interaction=FALSE))
#' 
#' If interaction=TRUE, then Y~ENET(X,A,X*A), and variables with estimated coefficients of 
#' zero in both the main effects (X) and treatment-interactions (X*A) are filtered. This 
#' aims to find variables that are prognostic and/or predictive. 
#' 
#' 
#' 2. \strong{ranger}: Wrapper function for the function "ranger" (ranger R package) to calculate
#' random forest based variable importance (VI) p-values. Here, for the test of VI>0,
#' variables are filtered if their one-sided p-value>=0.10. P-values are obtained
#' through subsampling based T-statistics (T=VI_j/SE(VE_j)) for feature j through the 
#' delete-d jackknife), as described in Ishwaran and Lu 2017. Used for continuous, binary, 
#' or survival outcomes. Default hyper-parameters are:
#' 
#' hyper=list(b=0.66, K=200, DF2=FALSE, FDR=FALSE, pval.thres=0.10)
#' 
#' where b=(\% of total data to sample; default=66\%), K=# of subsamples, FDR (FDR based
#' multiplicity correction for p-values), pval.thres=0.10 (adjust to change 
#' filtering threshold). DF2 fits Y~ranger(X, A, XA) and calculates the 
#' VI_2DF = VI_X+VI_XA, which is the variable importance of the main effect + the 
#' interaction effect (joint test). Var(VI_2DF) = Var(VI_X)+Var(VI_AX)+2cov(VI_X, VI_AX) 
#' where each component is calculated using the subsampling approach described above.
#' 
#' 
#' @examples
#' \donttest{
#' library(StratifiedMedicine)
#' ## Continuous ##
#' dat_ctns = generate_subgrp_data(family="gaussian")
#' Y = dat_ctns$Y
#' X = dat_ctns$X
#' A = dat_ctns$A
#'
#' # Fit ple_ranger directly (treatment-specific ranger models) #
#' mod1 = filter_train(Y, A, X, filter="filter_glmnet")
#' mod1$filter.vars
#' 
#' mod2 = filter_train(Y, A, X, filter="filter_glmnet", hyper=list(interaction=TRUE))
#' mod2$filter.vars
#' 
#' mod3 = filter_train(Y, A, X, filter="filter_ranger")
#' mod3$filter.vars
#' 
#' }
#'
#'
#' @export
#' @references Friedman, J., Hastie, T. and Tibshirani, R. (2008) Regularization Paths for
#'  Generalized Linear Models via Coordinate Descent,
#'  \url{https://web.stanford.edu/~hastie/Papers/glmnet.pdf} Journal of Statistical 
#'  Software, Vol. 33(1), 1-22 Feb 2010 Vol. 33(1), 1-22 Feb 2010.
#' @references Wright, M. N. & Ziegler, A. (2017). ranger: A fast implementation of 
#' random forests for high dimensional data in C++ and R. J Stat Softw 77:1-17. 
#' \url{https://doi.org/10.18637/jss.v077.i01}.
#' @references Ishwaran, H. Lu, M. (2017). Standard errors and confidence intervals 
#' for variable importance in random forest regression, classification, and survival.
#' Statistics in Medicine 2017. 
#' @seealso \code{\link{PRISM}}
#'
filter_train = function(Y, A, X, family="gaussian", filter="glmnet", hyper=NULL, ...) {
  
  if (filter %in% c("glmnet", "ranger")) {
    filter <- paste("filter", filter, sep="_")
  }
  if (is.Surv(Y) & family!="survival") {
    family <- "survival"
  }
  fit <- do.call(filter, append(list(Y=Y, A=A, X=X, family=family), hyper))
  
  res <- list(mod = fit$mod, filter.vars=fit$filter.vars, filter=filter)
  
  class(res) <- "filter_train"
  return(res)
}
