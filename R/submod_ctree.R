#' Subgroup Identification: Conditional Inference Trees (ctree)
#'
#' Uses the ctree (conditional inference trees) algorithm to identify subgroups
#' (Hothorn, Hornik, Zeileis 2006). Usable for continuous, binary, or survival outcomes.
#' Option to use the observed outcome or PLEs for subgroup identification.
#'
#' @param Y The outcome variable. Must be numeric or survival (ex; Surv(time,cens) )
#' @param A Treatment variable. (a=1,...A)
#' @param X Covariate matrix. Must be numeric.
#' @param Xtest Test set
#' @param mu_train Patient-level estimates (See PLE_models)
#' @param minbucket Minimum number of observations in a tree node.
#' Default = floor( dim(train)[1]*0.05  )
#' @param maxdepth Maximum depth of any node in the tree (default=4)
#' @param outcome_PLE If TRUE, use PLE as outcome (mu_train must contain PLEs).
#' @param family Outcome type ("gaussian", "binomial", "survival), default is "gaussian"
#' @param ... Any additional parameters, not currently passed through.
#'
#' @import partykit
#'
#' @return Trained ctree model.
#'  \itemize{
#'   \item mod - ctree model object
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
#'
#' res_ctree1 = submod_ctree(Y, A, X, Xtest=X, family="gaussian")
#' res_ctree2 = submod_ctree(Y, A, X, Xtest=X, family="gaussian", maxdepth=2, minsize=100)
#' plot(res_ctree1$mod)
#' plot(res_ctree2$mod)
#'
#'
#### CTREE ###
submod_ctree = function(Y, A, X, Xtest, mu_train, minbucket = floor( dim(X)[1]*0.05  ),
                        maxdepth = 4, outcome_PLE=FALSE, family="gaussian", ...){

  ## Use PLE as outcome? ##
  if (outcome_PLE==TRUE){
    Y = mu_train$PLE
    family = "gaussian"
  }
  ## Fit Model ##
  mod <- ctree(Y ~ ., data = X,
               control = ctree_control(minbucket=minbucket, maxdepth=maxdepth))
  res = list(mod=mod, family=family)
  class(res) = "submod_ctree"
  ## Return Results ##
  return(  res )
}

#' Predict submod: CTREE
#'
#' Predict subgroups and obtain subgroup-specific estimates, E(Y|X) or PLE(X), for a
#' trained ctree model (depends on if outcome_PLE argument)
#'
#' @param object Trained ctree model.
#' @param newdata Data-set to make predictions at.
#' @param ... Any additional parameters, not currently passed through.
#'
#' @import partykit
#'
#' @return Identified subgroups with subgroup-specific predictions of E(Y|X) or PLE(X).
#' \itemize{
#'   \item Subgrps - Identified subgroups
#'   \item pred - Predictions, E(Y|X) or PLE(X) by subgroup.
#'}
#' @examples
#' library(StratifiedMedicine)
#'
#' ## Continuous ##
#' dat_ctns = generate_subgrp_data(family="gaussian")
#' Y = dat_ctns$Y
#' X = dat_ctns$X
#' A = dat_ctns$A
#'
#' res_ctree1 = submod_ctree(Y, A, X, Xtest=X)
#' # Predict subgroups / estimates #
#' out = predict(res_ctree1, newdata=X)
#'
#' @method predict submod_ctree
#' @export
#'
predict.submod_ctree = function(object, newdata, ...){

  # Extract mod/family #
  mod = object$mod
  family = object$family
  ##  Predict Subgroups ##
  Subgrps = as.numeric( predict(mod, type="node", newdata = newdata) )
  ## Response Predictions ##
  if (family=="gaussian"){ type.fam = "response"   } # E(Y|X)
  if (family=="binomial"){ type.fam = "prob"   } # probability
  if (family=="survival"){ type.fam = "response" } # median survival
  pred = predict( mod, newdata = newdata, type = type.fam )
  ## Return Results ##
  return(  list(Subgrps=Subgrps, pred=pred) )
}
