#' Subgroup Identification: Conditional Inference Trees (CTREE)
#'
#' Uses the CTREE (conditional inference trees) algorithm to identify subgroups.
#' Usable for continuous and binary outcomes.
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
#' @return CTREE model, predictions, and identified subgroups.
#'  \itemize{
#'   \item mod - CTREE model object
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
#' res_ctree1 = submod_ctree(Y, A, X, Xtest=X, family="gaussian")
#' res_ctree2 = submod_ctree(Y, A, X, Xtest=X, family="gaussian", maxdepth=2, minsize=100)
#' plot(res_ctree1$mod)
#' plot(res_ctree2$mod)
#'
#' # Survival #
#' dat_surv = generate_subgrp_data(family="survival")
#' Y = dat_surv$Y
#' X = dat_surv$X
#' A = dat_surv$A
#' res_ctree3 = submod_ctree(Y, A, X, Xtest=X, family="survival")
#' plot(res_ctree3$mod)
#'
#' @seealso \code{\link{PRISM}}, \code{\link{ctree}}
#'
#### CTREE ###
submod_ctree = function(Y, A, X, Xtest, mu_train, minbucket = floor( dim(X)[1]*0.05  ),
                        maxdepth = 4, outcome_PLE=FALSE, family="gaussian", ...){

  ## Use PLE as outcome? ##
  if (outcome_PLE==TRUE){
    Y = mu_train$PLE
  }
  ## Fit Model ##
  mod <- ctree(Y ~ ., data = X,
               control = ctree_control(minbucket=minbucket, maxdepth=maxdepth))
  ##  Predict Subgroups for Train/Test ##
  Subgrps.train = as.numeric( predict(mod, type="node") )
  Subgrps.test = as.numeric( predict(mod, type="node", newdata = Xtest) )
  ## Response Predictions ##
  if (family=="gaussian"){ type.fam = "response"   } # E(Y|X)
  if (family=="binomial"){ type.fam = "prob"   } # probability
  if (family=="survival"){ type.fam = "response" } # median survival
  pred.train = predict( mod, newdata = data.frame(X), type = type.fam )
  pred.test = predict( mod, newdata = data.frame(Xtest), type = type.fam )
  ## Return Results ##
  return(  list(mod=mod, Subgrps.train=Subgrps.train, Subgrps.test=Subgrps.test,
                pred.train=pred.train, pred.test=pred.test) )
}
