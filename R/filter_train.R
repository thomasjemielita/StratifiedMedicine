#' Filter: Train Filter Model
#'
#' Wrapper function to train a filter model. Options include elastic net (glmnet) and random forest based 
#' variable importance (ranger). Used directly in PRISM. 
#'
#' @param Y The outcome variable. Must be numeric or survival (ex; Surv(time,cens) )
#' @param A Treatment variable. (a=1,...A)
#' @param X Covariate space.
#' @param family Outcome type ("gaussian", "binomial", "survival"). Default is "gaussian".
#' @param filter Filter function. Potentially reduces covariate space, 
#' (Y, A, X) ==> (Y, A, Xstar).
#' @param hyper Hyper-parameters for the filter model (must be list). Default is NULL.
#' @param ... Any additional parameters, not currently passed through.
#'
#'
#' @return Trained filter model and vector of variable names that pass the filter. 
#'
#'  \itemize{
#'   \item mod - trained model
#'   \item filter.vars - Variables that remain after filtering (could be all)
#' }
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
#' mod2 = filter_train(Y, A, X, filter="filter_ranger")
#' mod2$filter.vars
#' 
#' }
#'
#'
#' @export
#' @seealso \code{\link{PRISM}}
#'
filter_train = function(Y, A, X, family="gaussian", filter, hyper=NULL, ...) {
  
  fit <- do.call(filter, append(list(Y=Y, A=A, X=X, family=family), hyper))
  
  res <- list(mod = fit$mod, filter.vars=fit$filter.vars)
  
  class(res) <- "filter_train"
  return(res)
}
