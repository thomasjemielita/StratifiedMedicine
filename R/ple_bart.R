#' Patient-level Estimates: BART
#'
#' Uses the BART algorithm (from BART R package) to obtain patient-level estimates.
#' Used for continuous or binary outcomes, with output estimates of
#' E(Y|X,A=a) and E(Y|X,A=1)-E(Y|X,A=0) (PLE). In progress.
#'
#' @param Y The outcome variable. Must be numeric or survival (ex; Surv(time,cens) )
#' @param A Treatment variable. (a=1,...A)
#' @param X Covariate matrix. Must be numeric.
#' @param Xtest Test set
#' @param family Outcome type ("gaussian", "binomial"), default is "gaussian"
#' @param ... Any additional parameters, not currently passed through.
#'
#' @import BART
#'
#' @return Patient-level estimates (E(Y|X,1), E(Y|X,0), E(Y|X,1)-E(Y|X,0)) for train/test sets.
#'  \itemize{
#'   \item mods - trained model(s)
#'   \item mu_train - Patient-level estimates (training set)
#'   \item mu_test - Patient-level estimates (test set)
#' }
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
#'
#' # BART #
#' \donttest{
#' mod1 = ple_bart(train=train, Xtest=X)
#'
#' summary(mod1$mu_train$PLE)
#' }
#'
#' @seealso \code{\link{PRISM}}, \code{\link{BART}}
#'
#### BART ###
ple_bart = function(Y, A, X, Xtest, family="gaussian", ...){

  ## Generate counterfactual design matrices  ##
  Xtrain_0 = data.frame(A=0, X)
  Xtrain_1 = data.frame(A=1, X)
  Xtest_0 = data.frame(A=0, Xtest)
  Xtest_1 = data.frame(A=1, Xtest)

  x.train = data.frame(A,X)
  # Test: Include both train/test for faster predictions #
  X.FULL = rbind(Xtrain_0, Xtrain_1, Xtest_0, Xtest_1)
  ## BART ##
  set.seed(51351)
  bartFit = wbart(x.train = x.train, y.train = Y, x.test = X.FULL)
  n.tr = dim(Xtrain_0)[1]
  n.ts = dim(Xtest_0)[1]
  ### PLE Predictions: Train/Test ###
  mu_train = data.frame(mu1 = bartFit$yhat.test.mean[(n.tr+1):(2*n.tr)],
                        mu0 = bartFit$yhat.test.mean[1:n.tr])
  mu_train$PLE = mu_train$mu1-mu_train$mu0

  mu_test = data.frame(mu1 = bartFit$yhat.test.mean[(2*n.tr+n.ts+1):(2*n.tr+2*n.ts)],
                       mu0 = bartFit$yhat.test.mean[(2*n.tr+1):(2*n.tr+n.ts)] )
  mu_test$PLE = mu_test$mu1 - mu_test$mu0

  ## Return Results ##
  return( list(mods = bartFit, mu_train = mu_train, mu_test=mu_test) )
}
