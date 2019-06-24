#' Subgroup Identification: Train Model
#'
#' Wrapper function to train a subgroup model (submod). Used directly in PRISM and can
#' be used to directly fit a submod model by name.
#'
#' @param Y The outcome variable. Must be numeric or survival (ex; Surv(time,cens) )
#' @param A Treatment variable. (a=1,...A)
#' @param X Covariate matrix. Must be numeric.
#' @param Xtest Test set
#' @param mu_train Patient-level estimates (See PLE_models). Default=NULL
#' @param family Outcome type ("gaussian", "binomial", "survival"). Default="gaussian".
#' @param submod Subgroup identification (submod) function. Maps the observed data and/or
#' PLEs to subgroups.
#' @param hyper Hyper-parameters for submod (must be list). Default is NULL.
#' @param ... Any additional parameters, not currently passed through.
#'
#'
#' @return Trained subgroup model and subgroup predictions/estimates for train/test sets.
#'
#'  \itemize{
#'   \item mod - trained subgroup model
#'   \item Subgrps.train - Identified subgroups (training set)
#'   \item Subgrps.test - Identified subgroups (test set)
#'   \item pred.train - Predictions (training set)
#'   \item pred.test - Predictions (test set)
#'   \item Rules - Definitions for subgroups, if provided in submod output.
#' }
#' @examples
#' library(StratifiedMedicine)
#' ## Continuous ##
#' dat_ctns = generate_subgrp_data(family="gaussian")
#' Y = dat_ctns$Y
#' X = dat_ctns$X
#' A = dat_ctns$A
#'
#' # Fit submod_lmtree directly #
#' mod1 = submod_lmtree(Y, A, X, Xtest=X)
#' plot(mod1$mod)
#'
#' # Fit through submod_train wrapper #
#' mod2 = submod_train(Y=Y, A=A, X=X, Xtest=X, submod="submod_lmtree")
#' plot(mod2$mod)
#'
#'
#' @export
#' @seealso \code{\link{PRISM}}
#'
submod_train = function(Y, A, X, Xtest, mu_train=NULL, family="gaussian", submod,
                        hyper=NULL, ...){
  ## Fit submod ##
  mod = do.call( submod, append(list(Y=Y, A=A, X=X, Xtest=Xtest,
                                  family=family), hyper) )
  ### Train/Test Predictions ###
  ## If prior predictions are made: ##
  if (!is.null(mod$Subgrps.train)){
    Subgrps.train = mod$Subgrps.train
    pred.train = mod$pred.train
  }
  if (!is.null(mod$Subgrps.test)){
    Subgrps.test = mod$Subgrps.test
    pred.test = mod$pred.test
  }
  ## If no prior predictions are mode: ##
  if (is.null(mod$Subgrps.train)){
    out = predict(mod, newdata = X)
    Subgrps.train = out$Subgrps
    pred.train = out$pred
  }
  if (is.null(mod$Subgrps.test)){
    out = predict(mod, newdata = Xtest)
    Subgrps.test = out$Subgrps
    pred.test = out$pred
  }
  Rules = mod$Rules
  res = list(mod = mod$mod, Subgrps.train=Subgrps.train, Subgrps.test=Subgrps.test,
             pred.train=pred.train, pred.test=pred.test, Rules=Rules)
  class(res) = "submod_train"
  return(res)
}
