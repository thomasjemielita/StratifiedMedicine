#' Patient-level Estimates: Causal Forest
#'
#' Uses the causal forest algorithm (grf R package) to obtain patient-level estimates.
#' Used for continuous or binary outcomes, with output estimates of
#' E(Y|X,A=a) and E(Y|X,A=1)-E(Y|X,A=0) (PLE).
#'
#' @param Y The outcome variable. Must be numeric or survival (ex; Surv(time,cens) )
#' @param A Treatment variable. (a=1,...A)
#' @param X Covariate matrix. Must be numeric.
#' @param Xtest Test set
#' @param tune If TRUE, use grf automatic hyper-parameter tuning. If FALSE (default), no tuning.
#' @param num.trees Number of trees (default=500)
#' @param family Outcome type ("gaussian", "binomial"), default is "gaussian"
#' @param ... Any additional parameters, not currently passed through.
#'
#' @import grf
#'
#' @return Patient-level estimates (E(Y|X,1), E(Y|X,0), E(Y|X,1)-E(Y|X,0)) for train/test sets.
#'  \itemize{
#'   \item mods - trained model(s)
#'   \item mu_train - Patient-level estimates (training set)
#'   \item mu_test - Patient-level estimates (test set)
#' }
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
#'\donttest{
#' # Default #
#' mod1 = PLE_causal_forest(train=train, Xtest=X)
#' summary(mod1$mu_train$PLE)
#'
#' # Tune hyper-parameters (uncomment to run) #
#' mod2 = PLE_causal_forest(train=train, Xtest=X, tune=TRUE)
#' summary(mod2$mu_train$PLE)
#'
#' }
#'
#'
#' @export
#' @seealso \code{\link{PRISM}}, \code{\link{causal_forest}}

#### Causal_forest ###
PLE_causal_forest = function(Y, A, X, Xtest, tune=FALSE, num.trees=500, family="gaussian", ...){

  set.seed(5131)
  W = A
  ## Regression Forest: Y~X ##
  forest.Y = regression_forest(X, Y, ci.group.size=1, num.trees = min(500, num.trees) )
  Y.hat.train = predict(forest.Y)$predictions
  Y.hat.test = predict(forest.Y, newdata=Xtest)$predictions
  ## Regression Forest: W~X, If RCT ==> W is independent of X; use sample mean ##
  W.hat.train = rep( mean(W), dim(X)[1] )
  W.hat.test = rep( mean(W), dim(Xtest)[1] )
  # forest.W = regression_forest(X, W, ci.group.size=1, num.trees = min(500, num.trees) )
  # W.hat.train = predict(forest.W)$predictions
  # W.hat.test = predict(forest.W, newdata=Xtest)$predictions
  ## Causal Forest ##
  forest.CF = causal_forest(X, Y, W, tune.parameters = tune, num.trees=num.trees,
                        W.hat=W.hat.train, Y.hat=Y.hat.train)

  ## Predictions: Train/Test ##
  PLE.train = forest.CF$predictions
  mu0.train = Y.hat.train - W.hat.train * PLE.train
  mu1.train = Y.hat.train + (1 - W.hat.train) * PLE.train

  PLE.test = predict(forest.CF, newdata = Xtest)$predictions
  mu0.test = Y.hat.test - W.hat.test * PLE.test
  mu1.test = Y.hat.test + (1 - W.hat.test) * PLE.test

  # Combine #
  mu_train = data.frame(mu1 = mu1.train, mu0 = mu0.train,
                        PLE = PLE.train)
  mu_test = data.frame(mu1 = mu1.test, mu0 = mu0.test,
                       PLE = PLE.test)
  ## Return Results ##
  return( list(mods = list(forest.Y, forest.CF), mu_train = mu_train, mu_test = mu_test) )
}
