#' Subgroup Identification: Train Model
#'
#' Wrapper function to train a subgroup model (submod). Used directly in PRISM and can
#' be used to directly fit a submod model by name.
#'
#' @param Y The outcome variable. Must be numeric or survival (ex; Surv(time,cens) )
#' @param A Treatment variable. (a=1,...A)
#' @param X Covariate space.
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
#'   \item Rules - Definitions for subgroups, if provided in fitted submod output.
#' }
#' @examples
#' 
#' \donttest{
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
#' plot(mod2$fit$mod)
#'
#'}
#'
#' @export
#' @seealso \code{\link{PRISM}}
#'
submod_train = function(Y, A, X, Xtest, mu_train=NULL, family="gaussian", submod,
                        hyper=NULL, ...){
  ## Fit submod ##
  fit = do.call( submod, append(list(Y=Y, A=A, X=X, Xtest=Xtest, mu_train=mu_train,
                                  family=family), hyper) )
  ### Train/Test Predictions ###
  ## If prior predictions are made: ##
  if (!is.null(fit$Subgrps.train)){
    Subgrps.train = fit$Subgrps.train
    pred.train = fit$pred.train
  }
  if (!is.null(fit$Subgrps.test)){
    Subgrps.test = fit$Subgrps.test
    pred.test = fit$pred.test
  }
  ## If no prior predictions are made: ##
  if (is.null(fit$Subgrps.train)){
    out = fit$pred.fun(fit$mod, X=X)
    Subgrps.train = out$Subgrps
    pred.train = out$pred
  }
  if (is.null(fit$Subgrps.test)){
    out = fit$pred.fun(fit$mod, X=Xtest)
    Subgrps.test = out$Subgrps
    pred.test = out$pred
  }
 
  Rules = fit$Rules
  res = list(fit = fit, Subgrps.train=Subgrps.train, Subgrps.test=Subgrps.test,
             Rules=Rules)
  class(res) = "submod_train"
  return(res)
}

#' Subgroup Identification: Train Model (Predictions)
#'
#' Prediction function for the trained subgroup identification model (submod).
#'
#' @param object Trained submod model.
#' @param newdata Data-set to make predictions at (Default=NULL, predictions correspond
#' to training data).
#' @param ... Any additional parameters, not currently passed through.
#'
#' @return Identified subgroups with subgroup-specific predictions (depends on subgroup
#' model)
#' \itemize{
#'   \item Subgrps - Identified subgroups
#'   \item pred - Predictions, depends on subgroup model
#'}
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
#' out2 = predict(mod2)
#' plot(mod2$fit$mod)
#'
#' @method predict submod_train
#' @export
#'
predict.submod_train = function(object, newdata=NULL, ...){

  # preds = predict(object$fit, newdata=newdata)
  preds = object$fit$pred.fun(object$fit$mod, newdata)
  ## Return Results ##
  return(  list(Subgrps=preds$Subgrps, pred=preds$pred) )
}
