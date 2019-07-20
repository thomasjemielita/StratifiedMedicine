#' Subgroup Identification: Model-based partitioning (lmtree)
#'
#' Uses the lmtree (model-based partitioning, OLS) algorithm to identify subgroups
#' (Zeileis, Hothorn, Hornik 2008). Usable for continuous and binary outcomes.
#'
#' @param Y The outcome variable. Must be numeric or survival (ex; Surv(time,cens) )
#' @param A Treatment variable. (a=1,...A)
#' @param X Covariate space.
#' @param Xtest Test set
#' @param mu_train Patient-level estimates (See PLE_models)
#' @param minsize Minimum number of observations in a tree node.
#' Default = floor( dim(train)[1]*0.05  )
#' @param maxdepth Maximum depth of any node in the tree (default=4)
#' @param ... Any additional parameters, not currently passed through.
#'
#' @import partykit
#'
#' @return Trained lmtree model.
#'  \itemize{
#'   \item mod - lmtree model object
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
#' res_lmtree1 = submod_lmtree(Y, A, X, Xtest=X)
#' res_lmtree2 = submod_lmtree(Y, A, X, Xtest=X, maxdepth=2, minsize=100)
#' plot(res_lmtree1$mod)
#' plot(res_lmtree2$mod)
#'
#'
#### lmtree (MOB) ###
submod_lmtree = function(Y, A, X, Xtest, mu_train, minsize = floor( dim(X)[1]*0.05  ),
                         maxdepth = 4, ...){

  ## Fit Model ##
  mod <- lmtree(Y~A | ., data = X, maxdepth = maxdepth, minsize=minsize)

  res = list(mod=mod)
  class(res) = "submod_lmtree"
  ## Return Results ##
  return(  res )
}

#' Predict submod: lmtree
#'
#' Predict subgroups and obtain subgroup-specific estimates of E(Y|A=1)-E(Y|A=0).
#'
#' @param object Trained lmtree model.
#' @param newdata Data-set to make predictions at (Default=NULL, predictions correspond
#' to training data).
#' @param ... Any additional parameters, not currently passed through.
#'
#' @import partykit
#'
#' @return Identified subgroups with subgroup-specific predictions of E(Y|A=1)-E(Y|A=0).
#' \itemize{
#'   \item Subgrps - Identified subgroups
#'   \item pred - Predictions, E(Y|A=1)-E(Y|A=0) by subgroup.
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
#' res_lmtree1 = submod_lmtree(Y, A, X, Xtest=X)
#' # Predict subgroups / estimates #
#' out = predict(res_lmtree1)
#'
#'
#' @method predict submod_lmtree
#' @export
#'
predict.submod_lmtree = function(object, newdata=NULL, ...){

  # Extract model #
  mod = object$mod
  ##  Predict Subgroups ##
  Subgrps = as.numeric( predict(mod, type="node", newdata = newdata) )
  ### Predict E(Y|A=1,S=s)-E(Y|A=0,S=s) (based on lmtree) ##
  hold.dat = data.frame(Subgrps = Subgrps, pred = NA)
  for (s in unique(Subgrps)){
    hold = summary(mod)[[as.character(s)]]
    hold.dat$pred[hold.dat$Subgrps==s] = coefficients(hold)[2,1]
  }
  pred = hold.dat$pred

  ## Return Results ##
  return(  list(Subgrps=Subgrps, pred=pred) )
}
