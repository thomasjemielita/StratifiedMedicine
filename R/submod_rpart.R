#' Subgroup Identification: CART (rpart)
#'
#' Uses the CART algorithm (rpart) to identify subgroups. Usable for continuous and binary outcomes.
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
#' Else use observed outcome Y
#' @param family Outcome type ("gaussian", "binomial"), default is "gaussian"
#' @param ... Any additional parameters, not currently passed through.
#'
#' @import rpart
#'
#' @return rpart model, predictions, and identified subgroups.
#'  \itemize{
#'   \item mod - rpart model as partykit object
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
#' res_rpart1 = submod_rpart(Y, A, X, Xtest=X, family="gaussian")
#' res_rpart2 = submod_rpart(Y, A, X, Xtest=X, maxdepth=2, minbucket=100, family="gaussian")
#' plot(res_rpart1$mod)
#' plot(res_rpart2$mod)
#'
#' @seealso \code{\link{PRISM}}, \code{\link{rpart}}
#'
## CART(rpart) ###
submod_rpart = function(Y, A, X, Xtest, mu_train, minbucket = floor( dim(X)[1]*0.05  ),
                       maxdepth = 4, outcome_PLE=FALSE, family="gaussian", ...){

  ## Use PLE as outcome? #
  if (outcome_PLE==TRUE){
    Y = mu_train$PLE
  }
  ## Fit Model ##
  mod <- rpart(Y ~ ., data = X,
               control = rpart.control(minbucket=minbucket, maxdepth=maxdepth))
  mod = as.party(mod)
  ##  Predict Subgroups for Train/Test ##
  Subgrps.train = as.numeric( predict(mod, type="node") )
  Subgrps.test = as.numeric( predict(mod, type="node", newdata = Xtest) )
  ## Response Predictions ##
  if (family=="gaussian"){ type.fam = "response"   } # E(Y|X)
  if (family=="binomial"){ type.fam = "prob"   } # probability
  pred.train = predict( mod, newdata = data.frame(X), type = type.fam )
  pred.test = predict( mod, newdata = data.frame(Xtest), type = type.fam )
  ## Return Results ##
  return(  list(mod=mod, Subgrps.train=Subgrps.train, Subgrps.test=Subgrps.test,
                pred.train=pred.train, pred.test=pred.test) )
}
