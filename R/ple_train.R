#' Patient-level Estimates: Train Model
#'
#' Wrapper function to train a patient-level estimate (ple) model. Used directly in PRISM
#' and can be used to directly fit a ple model by name.
#'
#' @param Y The outcome variable. Must be numeric or survival (ex; Surv(time,cens) )
#' @param A Treatment variable. (a=1,...A)
#' @param X Covariate space.
#' @param Xtest Test set
#' @param family Outcome type ("gaussian", "binomial", "survival"). Default is "gaussian".
#' @param ple PLE (Patient-Level Estimate) function. Maps the observed data to PLEs.
#' (Y,A,X) ==> PLE(X).
#' @param hyper Hyper-parameters for the ple model (must be list). Default is NULL.
#' @param ... Any additional parameters, not currently passed through.
#'
#'
#' @return Trained ple models and patient-level estimates for train/test sets. For
#' family="gaussian" or "binomial", output estimates of
#' (E(Y|X,A=1), E(Y|X,A=0), E(Y|X,A=1)-E(Y|X,A=0)). For survival, output estimates of
#' (HR(X,A=1), HR(X,A=0), HR(X, A=1)-HR(X, A=0)).
#'
#'  \itemize{
#'   \item mods - trained model(s)
#'   \item mu_train - Patient-level estimates (training set)
#'   \item mu_test - Patient-level estimates (test set)
#' }
#' @examples
#' library(StratifiedMedicine)
#' ## Continuous ##
#' dat_ctns = generate_subgrp_data(family="gaussian")
#' Y = dat_ctns$Y
#' X = dat_ctns$X
#' A = dat_ctns$A
#'
#' # Fit ple_ranger directly (treatment-specific ranger models) #
#' mod1 = ple_ranger(Y, A, X, Xtest=X)
#' summary(mod1$mu_train)
#'
#' # Fit through ple_train wrapper #
#' mod2 = ple_train(Y=Y, A=A, X=X, Xtest=X, ple="ple_ranger" )
#' summary(mod2$mu_train)
#'
#'
#' @export
#' @seealso \code{\link{PRISM}}
#'
ple_train = function(Y, A, X, Xtest, family="gaussian", ple, hyper=NULL, ...){
  ## Fit ple model ##
  fit = do.call( ple, append(list(Y=Y, A=A, X=X, Xtest=Xtest,
                                  family=family), hyper) )
  ### Train/Test Predictions ###
  ## If prior predictions are made: ##
  if (!is.null(fit$mu_train)){
    mu_train = fit$mu_train
  }
  if (!is.null(fit$mu_test)){
    mu_test = fit$mu_test
  }
  ## If no prior predictions are mode: ##
  if (is.null(fit$mu_train)){
    mu_train = predict(fit)
  }
  if (is.null(fit$mu_test)){
    mu_test = predict(fit, newdata = data.frame(Xtest) )
  }
  res = list(fit = fit, mu_train=mu_train, mu_test=mu_test)
  class(res) = "ple_train"
  return(res)
}

#' Patient-level Estimates Model: Prediction
#'
#' Prediction function for the trained patient-level estimate (ple) model.
#'
#' @param object Trained ple model.
#' @param newdata Data-set to make predictions at (Default=NULL, predictions correspond
#' to training data).
#' @param ... Any additional parameters, not currently passed through.
#'
#' @return Data-frame with predictions (depends on trained ple model).
#'
#' @examples
#' library(StratifiedMedicine)
#' ## Continuous ##
#' dat_ctns = generate_subgrp_data(family="gaussian")
#' Y = dat_ctns$Y
#' X = dat_ctns$X
#' A = dat_ctns$A
#'
#'
#' # Fit through ple_train wrapper #
#' mod2 = ple_train(Y=Y, A=A, X=X, Xtest=X, ple="ple_ranger" )
#' summary(mod2$mu_train)
#'
#' res2 = predict(mod2) # newdata=NULL, training data #
#' res3 = predict(mod2, newdata=X) # test data #
#' summary(res2)
#' summary(res3)
#'
#' @method predict ple_train
#' @export
#' @seealso \code{\link{PRISM}}
#'
predict.ple_train = function(object, newdata=NULL, ...){

  mu_hat = predict(object$fit, newdata = newdata )

  return(mu_hat)
}
