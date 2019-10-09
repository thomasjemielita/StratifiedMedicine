#' Subgroup Identification: Optimal Treatment Regime (through ctree)
#'
#' For continuous, binary, or survival outcomes, regress I(PLE>thres)~X with
#' weights=abs(PLE) in ctree.
#'
#' @param Y The outcome variable. Must be numeric or survival (ex; Surv(time,cens) )
#' @param A Treatment variable. (a=1,...A)
#' @param X Covariate space.
#' @param Xtest Test set
#' @param mu_train Patient-level estimates (See PLE_models)
#' @param minbucket Minimum number of observations in a tree node.
#' Default = floor( dim(train)[1]*0.05  )
#' @param maxdepth Maximum depth of any node in the tree (default=4)
#' @param thres Threshold for PLE, ex: I(PLE>thres). Default is ">0". Direction can be
#' reversed and can include equality sign (ex: "<=")
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

#### OTR: I(PLE>thres) ~ X, weights = abs(PLE) ###
submod_otr = function(Y, A, X, Xtest, mu_train, minbucket = floor( dim(X)[1]*0.10  ),
                      maxdepth = 4, thres=">0", ...){
  ## Set up data ##
  ind_PLE <- eval(parse(text=paste("ifelse(mu_train$PLE", thres, ", 1, 0)")))
  w_PLE <- abs(mu_train$PLE)
  hold <- data.frame(ind_PLE, X)
  # Fit Model #
  mod <- suppressWarnings( ctree(ind_PLE ~ ., data = hold, weights = w_PLE,
                                 control = ctree_control(minbucket=minbucket,
                                                         maxdepth=maxdepth)) )
  # Prediction Function #
  pred.fun <- function(mod, X=NULL, type="subgrp"){
    Subgrps <- NULL; pred <- NULL;
    Subgrps = as.numeric( predict(mod, type="node", newdata = X) )
    if (type=="all"){
      pred = data.frame(Subgrps=Subgrps,
                        as.numeric( predict(mod, newdata = X, type="response") ) )
    }
    return( list(Subgrps=Subgrps, pred=pred) )
  }
  ## Return Results ##
  res <- list(mod=mod, pred.fun=pred.fun)
  class(res) <- "submod_otr"
  return(  res )
}
