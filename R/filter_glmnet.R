#' Filter: Elastic Net (glmnet)
#'
#' Filter variables through elastic net (Zou and Hastie 2005). Default is to regress
#' Y~X (search for prognostic variables). Variables with estimated coefficients of zero
#' (depends on lambda choice; default is lambda.min) are filtered. Usable for continuous,
#' binary, and survival outcomes.
#'
#' @param Y The outcome variable. Must be numeric or survival (ex; Surv(time,cens) )
#' @param A Treatment variable. (a=1,...A)
#' @param X Covariate space.
#' @param lambda Lambda for elastic net model (default="lambda.min"). Other options include
#' "lambda.1se" and fixed values
#' @param family Outcome type ("gaussian", "binomial", "survival"), default is "gaussian"
#' @param interaction Regress Y~X+A+A*X (interaction between covariates and treatment)?
#' Default is FALSE. If TRUE, variables with zero coefficients (both X and X*A terms)
#' are filtered.
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
#' # Default: Regress Y~X (search for prognostic factors) #
#' mod1 = filter_glmnet(Y, A, X)
#' mod2 = filter_glmnet(Y, A, X, lambda = "lambda.min") # same as default
#' mod3 = filter_glmnet(Y, A, X, lambda = "lambda.1se")
#' mod1$filter.vars
#' mod2$filter.vars
#' mod3$filter.vars
#'
#' # Interaction=TRUE; Regress Y~X+A+X*A (search for prognostic and/or predictive) #
#' mod4 = filter_glmnet(Y, A, X, interaction=TRUE)
#' mod4$filter.vars

##### Elastic net (glmnet): Y~X ######
filter_glmnet = function(Y, A, X, lambda="lambda.min", family="gaussian",
                         interaction=FALSE,...){

  ## Model Matrix #
  X.mat = X
  colnames(X.mat) = paste(colnames(X.mat), "_REMOVE_", sep="")
  X.mat = model.matrix(~., data = X.mat )
  X.mat = X.mat[, colnames(X.mat) != "(Intercept)"]
  W = X.mat

  if (interaction){
    X_inter = X.mat*A
    colnames(X_inter) = paste(colnames(X.mat), "_A", sep="")
    W = cbind(X.mat, A, X_inter)
  }

  ##### Elastic Net on estimated ITEs #####
  if (family=="survival") { family = "cox" }
  mod <- cv.glmnet(x = W, y = Y, nlambda = 100, alpha=0.5, family=family)

  ### Extract filtered variable based on lambda ###
  VI <- coef(mod, s = lambda)[,1]
  VI = VI[ names(VI) != "(Intercept)" ]
  # Extract variables that pass the filter ##
  filter.vars = names(VI[VI!=0])
  filter.vars = unique( gsub("_REMOVE_.*","",filter.vars) )

  return( list(mod=mod, filter.vars=filter.vars) )
}
