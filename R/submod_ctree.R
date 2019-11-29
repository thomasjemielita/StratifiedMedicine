#' Subgroup Identification: Conditional Inference Trees (ctree)
#'
#' Uses the ctree (conditional inference trees) algorithm to identify subgroups
#' (Hothorn, Hornik, Zeileis 2006). Usable for continuous, binary, or survival outcomes.
#' Option to use the observed outcome or PLEs for subgroup identification.
#'
#' @param Y The outcome variable. Must be numeric or survival (ex; Surv(time,cens) )
#' @param A Treatment variable. (a=1,...A)
#' @param X Covariate space.
#' @param Xtest Test set
#' @param mu_train Patient-level estimates (See PLE_models)
#' @param alpha Significance level for variable selection (default=0.05)
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
submod_ctree = function(Y, A, X, Xtest, mu_train, alpha=0.05,
                        minbucket = floor( dim(X)[1]*0.10  ),
                        maxdepth = 4, outcome_PLE=FALSE, family="gaussian", ...){

  ## Use PLE as outcome? ##
  if (outcome_PLE==TRUE){
    Y <- mu_train$PLE
    family <- "gaussian"
  }
  # Fit Model #
  mod <- ctree(Y ~ ., data = X,
               control = ctree_control(alpha=alpha, 
                                       minbucket=minbucket, maxdepth=maxdepth))
  
  # Prediction Function #
  pred.fun <- function(mod, X=NULL, type="subgrp"){
    Subgrps <- NULL
    pred <- NULL
    Subgrps <- as.numeric( predict(mod, type="node", newdata = X) )
    if (type=="all"){
      pred <- data.frame(Subgrps=Subgrps,
                        mu=as.numeric( predict(mod, newdata = X, type="response")) )
    }
    return( list(Subgrps=Subgrps, pred=pred) )
  }
  res <- list(mod=mod, pred.fun=pred.fun)
  class(res) = "submod_ctree"
  ## Return Results ##
  return(  res )
}