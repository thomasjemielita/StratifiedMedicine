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
#' @param mod.A Model for estimating P(A|X). Default is "mean" calculates the sample mean.
#' If mod.A="RF", estimate P(A|X) using regression_forest (applicable for non-RCTs).
#' @param ... Any additional parameters, not currently passed through.
#'
#' @import grf
#'
#' @return Trained causal_forest and regression_forest models.
#'  \itemize{
#'   \item mods - trained model(s)
#' }
#' @examples
#' library(StratifiedMedicine)
#'
#' ## Continuous ##
#' dat_ctns = generate_subgrp_data(family="gaussian")
#' Y = dat_ctns$Y
#' X = dat_ctns$X
#' A = dat_ctns$A
#'
#'\donttest{
#' # Default #
#' mod1 = ple_causal_forest(Y, A, X, Xtest=X)
#' summary(mod1$mu_train)
#'
#' }
#'
#'
#' @export
#' @seealso \code{\link{PRISM}}, \code{\link{causal_forest}}

#### Causal_forest ###
ple_causal_forest = function(Y, A, X, Xtest, tune=FALSE, num.trees=500, family="gaussian",
                             mod.A = "mean", ...){

  set.seed(5131)
  ## Regression Forest: Y~X ##
  forest.Y = regression_forest(X, Y, ci.group.size=1, num.trees = min(500, num.trees) )
  Y.hat.train = predict(forest.Y)$predictions
  Y.hat.test = predict(forest.Y, newdata=Xtest)$predictions
  ## Regression Forest: W~X, If RCT ==> W is independent of X; use sample mean ##
  forest.A = mean(A)
  A.hat.train = rep( mean(A), dim(X)[1] )
  A.hat.test = rep( mean(A), dim(Xtest)[1] )
  # forest.A = regression_forest(X, W, ci.group.size=1, num.trees = min(500, num.trees) )
  # A.hat.train = predict(forest.A)$predictions
  # A.hat.test = predict(forest.A, newdata=Xtest)$predictions
  ## Causal Forest ##
  forest.CF = causal_forest(X, Y, W=A, tune.parameters = tune, num.trees=num.trees,
                        W.hat=A.hat.train, Y.hat=Y.hat.train)

  res = list(mods = list(forest.Y=forest.Y, forest.A=forest.A, forest.CF=forest.CF) )
  class(res) = "ple_causal_forest"
  ## Return Results ##
  return( res )
}

#' Predict Patient-level Estimates: Causal Forest
#'
#' Get estimates of (E(Y|X,A=1), E(Y|X,A=0), E(Y|X,A=1)-E(Y|X,A=0)) using trained
#' regression_forest and causal_forest model(s).
#'
#' @param object Trained random forest (ranger) model(s).
#' @param newdata Data-set to make predictions at. For training data, should
#' include covariates (X) and treatment (A) for oob predictions.
#' @param ... Any additional parameters, not currently passed through.
#'
#' @import grf
#'
#' @return Data-frame with predictions of (E(Y|X,1), E(Y|X,0), E(Y|X,1)-E(Y|X,0))
#'
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
#' mod1 = ple_causal_forest(Y, A, X, Xtest=X)
#' summary(mod1$mu_train)
#' summary(predict(mod1, newdata=data.frame(A,X)))
#' summary(predict(mod1, newdata=data.frame(X)))
#' }
#'
#' @method predict ple_causal_forest
#' @export
#'
#### Predict: ple_causal_forest ####
predict.ple_causal_forest = function(object, newdata, ...){

  forest.Y = object$mods$forest.Y
  forest.A = object$mods$forest.A
  forest.CF = object$mods$forest.CF
  ## Use oob if possible (relevant for training data) ##
  if ( "A" %in% names(newdata) ){
    Y_hat = forest.Y$predictions
    ## Regression Forest: W~X, If RCT ==> W is independent of X; use sample mean ##
    if (is.numeric(forest.A)){
      A_hat = rep( forest.A, dim(newdata)[1]  )
    }
    if (is.list(forest.A)){
      A_hat = forest.A$predictions
    }
    PLE_hat = forest.CF$predictions
  }
  if (!("A" %in% names(newdata))){
    # newdata = newdata[,!(colnames(newdata) %in% "A") ]
    Y_hat = predict(forest.Y, newdata)$predictions
    ## Regression Forest: W~X, If RCT ==> W is independent of X; use sample mean ##
    if (is.numeric(forest.A)){
      A_hat = rep( forest.A, dim(newdata)[1]  )
    }
    if (is.list(forest.A)){
      A_hat = predict(forest.A, newdata)
    }
    PLE_hat = predict(forest.CF, newdata)$predictions
  }
  mu0_hat = Y_hat - A_hat * PLE_hat
  mu1_hat = Y_hat + (1 - A_hat) * PLE_hat
  mu_hat = data.frame(mu1 = mu1_hat, mu0 = mu0_hat, PLE = PLE_hat)
  return( mu_hat  )
}


