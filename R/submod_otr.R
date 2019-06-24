#' Subgroup Identification: Optimal Treatment Regime (through CTREE)
#'
#' For continuous, binary, or survival outcomes, regress I(PLE>thres)~X with
#' weights=abs(PLE) in CTREE.
#'
#' @param Y The outcome variable. Must be numeric or survival (ex; Surv(time,cens) )
#' @param A Treatment variable. (a=1,...A)
#' @param X Covariate matrix. Must be numeric.
#' @param Xtest Test set
#' @param mu_train Patient-level estimates (See PLE_models)
#' @param minbucket Minimum number of observations in a tree node.
#' Default = floor( dim(train)[1]*0.05  )
#' @param maxdepth Maximum depth of any node in the tree (default=4)
#' @param thres Threshold for I(PLE>thres) (default=0)
#' @param ... Any additional parameters, not currently passed through.
#'
#' @import partykit
#'
#' @return Trained ctree (optimal treatment regime) model.
#'  \itemize{
#'   \item mod - tree (OTR) model object
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
#' \donttest{
#' ## Estimate PLEs (through Ranger) ##
#' res.ple = ple_model(Y, A, X, Xtest=X, family="gaussian", ple="ple_ranger")
#'
#' ## Fit OTR Subgroup Model ##
#' res_otr = submod_otr(Y, A, X, Xtest=X, mu_train = res.ple$mu_train)
#' plot(res_otr$mod)
#' }
#'
#' @seealso \code{\link{PRISM}}, \code{\link{ctree}}

#### OTR: I(PLE>thres) ~ X, weights = abs(PLE) ###
submod_otr = function(Y, A, X, Xtest, mu_train, minbucket = floor( dim(X)[1]*0.05  ),
                      maxdepth = 4, thres=0, ...){
  ## Set up data ##
  ind_PLE = ifelse(mu_train$PLE>thres, 1, 0)
  w_PLE = abs(mu_train$PLE)
  hold = data.frame(ind_PLE, X)
  mod <- suppressWarnings( ctree(ind_PLE ~ ., data = hold, weights = w_PLE,
                                 control = ctree_control(minbucket=minbucket,
                                                         maxdepth=maxdepth)) )

  res = list(mod=mod)
  class(res) = "submod_otr"
  ## Return Results ##
  return(  res )
}

#' Predict submod: OTR CTREE
#'
#' Predict subgroups and obtain subgroup-specific estimates, P(PLE>thres), for a
#' trained ctree OTR model.
#'
#' @param object Trained ctree model.
#' @param newdata Data-set to make predictions at.
#' @param ... Any additional parameters, not currently passed through.
#'
#' @import partykit
#'
#' @return Identified subgroups with subgroup-specific predictions of P(PLE>thres).
#' \itemize{
#'   \item Subgrps - Identified subgroups
#'   \item pred - Predictions, P(PLE>thres) by identified subgroup.
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
#' \donttest{
#' ## Estimate PLEs (through Ranger) ##
#' res.ple = ple_model(Y, A, X, Xtest=X, family="gaussian", ple="ple_ranger")
#'
#' ## Fit OTR Subgroup Model ##
#' res_otr = submod_otr(Y, A, X, Xtest=X, mu_train = res.ple$mu_train)
#' out = predict(res_otr, newdata=X)
#' plot(res_otr$mod)
#' }
#'
#' @method predict submod_otr
#' @export
#'
predict.submod_otr = function(object, newdata, ...){

  # Extract mod/family #
  mod = object$mod
  family = object$family
  ##  Predict Subgroups ##
  Subgrps = as.numeric( predict(mod, type="node", newdata = newdata) )
  ## Predict P(PLE>thres|X) ##
  pred = as.numeric( predict(mod, newdata = newdata) )
  ## Return Results ##
  return(  list(Subgrps=Subgrps, pred=pred) )
}

