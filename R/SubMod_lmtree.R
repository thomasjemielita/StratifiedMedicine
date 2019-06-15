#' Subgroup Identification: Model-based partitioning (lmtree)
#'
#' Uses the lmtree (model-based partitioning, OLS) algorithm to identify subgroups.
#' Usable for continuous and binary outcomes.
#'
#' @param Y The outcome variable. Must be numeric or survival (ex; Surv(time,cens) )
#' @param A Treatment variable. (a=1,...A)
#' @param X Covariate matrix. Must be numeric.
#' @param Xtest Test set
#' @param mu_train Patient-level estimates (See PLE_models)
#' @param minsize Minimum number of observations in a tree node.
#' Default = floor( dim(train)[1]*0.05  )
#' @param maxdepth Maximum depth of any node in the tree (default=4)
#' @param ... Any additional parameters, not currently passed through.
#'
#' @import partykit
#'
#' @return lmtree model, predictions, identified subgroups, and subgroup rules/defitions.
#'  \itemize{
#'   \item mod - lmtree model object
#'   \item Subgrps.train - Identified subgroups (training set)
#'   \item Subgrps.test - Identified subgroups (test set)
#'   \item pred.train - Predictions (training set)
#'   \item pred.test - Predictions (test set)
#'   \item Rules - Subgroups rules/definitions
#' }
#'
#' @export
#' @examples
#' library(StratifiedMedicine)
#'
#' ## Continuous ##
#' dat_ctns = generate_subgrp_data(family="gaussian")
#' Y = dat_ctns$Y
#' X = dat_ctns$X
#' A = dat_ctns$A
#' train = data.frame(Y, A, X)
#' # Outcome/treatment must be labeled as Y/A #
#'
#' res_lmtree1 = SubMod_lmtree(Y, A, X, Xtest=X)
#' res_lmtree2 = SubMod_lmtree(Y, A, X, Xtest=X, maxdepth=2, minsize=100)
#' plot(res_lmtree1$mod)
#' plot(res_lmtree2$mod)
#'
#' @seealso \code{\link{PRISM}}, \code{\link{lmtree}}
#'
#### lmtree (MOB) ###
SubMod_lmtree = function(Y, A, X, Xtest, mu_train, minsize = floor( dim(X)[1]*0.05  ),
                         maxdepth = 4, ...){

  ## Fit Model ##
  mod <- lmtree(Y~A | ., data = X, maxdepth = maxdepth, minsize=minsize)
  ##  Predict Subgroups for Train/Test ##
  Subgrps.train = as.numeric( predict(mod, type="node") )
  Subgrps.test = as.numeric( predict(mod, type="node", newdata = Xtest) )
  Rules = list_rules(mod)
  if (length(unique(Subgrps.train))==1){
    Rules = data.frame(Subgrps = unique(Subgrps.train), Rules = "All" )
  }
  if (length(unique(Subgrps.train))>1){
    Rules = data.frame(Subgrps = as.numeric(names(Rules)), Rules = as.character(Rules) )
  }
  ## Predict E(Y|X=x, A=1)-E(Y|X=x,A=0) ##
  pred.train = predict( mod, newdata = data.frame(A=1, X) ) -
    predict( mod, newdata = data.frame(A=0, X) )
  pred.test =  predict( mod, newdata = data.frame(A=1, Xtest) ) -
    predict( mod, newdata = data.frame(A=0, Xtest) )
  ## Return Results ##
  return(  list(mod=mod, Subgrps.train=Subgrps.train, Subgrps.test=Subgrps.test,
                pred.train=pred.train, pred.test=pred.test, Rules=Rules) )
}
