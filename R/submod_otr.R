#' Subgroup Identification: Optimal Treatment Regime (through CTREE)
#'
#' For continuous, binary, or survival outcomes, regress I(PLE>thres)~X with weights=abs(PLE) in CTREE
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
#' @return CTREE (OTR) model, predictions, and identified subgroups.
#'  \itemize{
#'   \item mod - CTREE (OTR) model object
#'   \item Subgrps.train - Identified subgroups (training set)
#'   \item Subgrps.test - Identified subgroups (test set)
#'   \item pred.train - Predictions (training set)
#'   \item pred.test - Predictions (test set)
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
#' mod_ple = ple_ranger(Y, A, X, Xtest=X)
#'
#' ## Fit OTR Subgroup Model ##
#' res_otr = submod_otr(Y, A, X, Xtest=X, mu_train = mod_ple$mu_train)
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
  ##  Predict Subgroups for Train/Test ##
  Subgrps.train = as.numeric( predict(mod, type="node") )
  Subgrps.test = as.numeric( predict(mod, type="node", newdata = Xtest) )
  ## Predict E(Y|X=x, A=1)-E(Y|X=x,A=0) ##
  pred.train = as.numeric( predict(mod) )
  pred.test = as.numeric( predict(mod, newdata = Xtest) )
  ## Return Results ##
  return(  list(mod=mod, Subgrps.train=Subgrps.train, Subgrps.test=Subgrps.test,
                pred.train=pred.train, pred.test=pred.test) )
}
