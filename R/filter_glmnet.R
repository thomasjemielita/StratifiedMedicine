#' Filter: Elastic Net (glmnet)
#'
#' Filter variables through elastic net (Zou and Hastie 2005). Variables with estimated
#' coefficients of zero (depends on lambda choice; default is lambda.min) are filtered.
#' Usable for continuous, binary, and survival outcomes.
#'
#' @param Y The outcome variable. Must be numeric or survival (ex; Surv(time,cens) )
#' @param A Treatment variable. (a=1,...A)
#' @param X Covariate matrix. Must be numeric.
#' @param lambda Lambda for elastic net model (default="lambda.min"). Other options include
#' "lambda.1se" and fixed values
#' @param family Outcome type ("gaussian", "binomial", "survival"), default is "gaussian"
#' @param ... Any additional parameters, not currently passed through.
#'
#' @return Filter model and variables that remain after filtering.
#'  \itemize{
#'   \item mod - Filtering model
#'   \item filter.vars - Variables that remain after filtering (could be all)
#' }
#' @export
#' @examples
#'
#' library(StratifiedMedicine)
#'
#' ## Continuous ##
#' dat_ctns = generate_subgrp_data(family="gaussian")
#' Y = dat_ctns$Y
#' X = dat_ctns$X
#' A = dat_ctns$A
#'
#' mod1 = filter_glmnet(Y, A, X)
#' mod2 = filter_glmnet(Y, A, X, lambda = "lambda.min") # same as default
#' mod3 = filter_glmnet(Y, A, X, lambda = "lambda.1se")
#' mod1$filter.vars
#' mod2$filter.vars
#' mod3$filter.vars

##### Elastic net (glmnet): Y~X ######
filter_glmnet = function(Y, A, X, lambda="lambda.min", family="gaussian", ...){

  ## Model matrix X matrix #
  X = model.matrix(~. -1, data = X )

  ##### Elastic Net on estimated ITEs #####
  set.seed(6134)
  if (family=="survival") { family = "cox"  }
  mod <- cv.glmnet(x = X, y = Y, nlambda = 100, alpha=0.5, family=family)

  ### Extract filtered variable based on lambda ###
  VI <- coef(mod, s = lambda)[,1]
  VI = VI[-1]
  filter.vars = names(VI[VI!=0])
  return( list(mod=mod, filter.vars=filter.vars) )
}
