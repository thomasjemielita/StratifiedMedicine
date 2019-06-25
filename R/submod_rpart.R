#' Subgroup Identification: CART (rpart)
#'
#' Uses the CART algorithm (rpart) to identify subgroups. Usable for continuous and binary
#' outcomes. Option to use the observed outcome or PLEs for subgroup identification.
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
#' Else use observed outcome Y
#' @param family Outcome type ("gaussian", "binomial"), default is "gaussian"
#' @param ... Any additional parameters, not currently passed through.
#'
#' @return Trained rpart (CART).
#'  \itemize{
#'   \item mod - rpart model as partykit object
#' }
#'
#' @export
#' @examples
#'
#' \donttest{
#' library(StratifiedMedicine)
#'
#' ## Continuous ##
#' dat_ctns = generate_subgrp_data(family="gaussian")
#' Y = dat_ctns$Y
#' X = dat_ctns$X
#' A = dat_ctns$A
#'
#' require(rpart)
#' res_rpart1 = submod_rpart(Y, A, X, Xtest=X)
#' res_rpart2 = submod_rpart(Y, A, X, Xtest=X, maxdepth=2, minbucket=100)
#' plot(res_rpart1$mod)
#' plot(res_rpart2$mod)
#' }
#'
#'
## CART(rpart) ###
submod_rpart = function(Y, A, X, Xtest, mu_train, minbucket = floor( dim(X)[1]*0.05  ),
                       maxdepth = 4, outcome_PLE=FALSE, family="gaussian", ...){

  if (!requireNamespace("rpart", quietly = TRUE)) {
    stop("Package rpart needed for submod_rpart. Please install.")
  }

  ## Use PLE as outcome? #
  if (outcome_PLE==TRUE){
    Y = mu_train$PLE
  }
  ## Fit Model ##
  mod <- rpart::rpart(Y ~ ., data = X,
               control = rpart::rpart.control(minbucket=minbucket, maxdepth=maxdepth))
  mod = as.party(mod)
  res = list(mod=mod, family=family)
  class(res) = "submod_rpart"
  ## Return Results ##
  return(  res )
}

#' Predict submod: rpart
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
#'
#' \donttest{
#' library(StratifiedMedicine)
#'
#' ## Continuous ##
#' dat_ctns = generate_subgrp_data(family="gaussian")
#' Y = dat_ctns$Y
#' X = dat_ctns$X
#' A = dat_ctns$A
#'
#' res_rpart = submod_rpart(Y, A, X, Xtest=X)
#' # Predict subgroups / estimates #
#' out = predict(res_rpart, newdata=X)
#' }
#'
#' @method predict submod_rpart
#' @export
#'
predict.submod_rpart = function(object, newdata, ...){

  # Extract mod/family #
  mod = object$mod
  family = object$family
  ##  Predict Subgroups ##
  Subgrps = as.numeric( predict(mod, type="node", newdata = newdata) )
  ## Response Predictions ##
  if (family=="gaussian"){ type.fam = "response"   } # E(Y|X)
  if (family=="binomial"){ type.fam = "prob"   } # probability
  pred = predict( mod, newdata = newdata, type = type.fam )
  ## Return Results ##
  return(  list(Subgrps=Subgrps, pred=pred) )
}
