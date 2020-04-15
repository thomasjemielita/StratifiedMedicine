#' Patient-level Estimates: Causal Forest
#'
#' Uses the causal forest algorithm (grf R package) to obtain patient-level estimates, 
#' E(Y|A=1), E(Y|A=0), and E(Y|A=1)-E(Y|A=0). Usable for continuous or binary outcomes.
#'
#' @param Y The outcome variable. Must be numeric or survival (ex; Surv(time,cens) )
#' @param A Treatment variable. (a=1,...A)
#' @param X Covariate space.
#' @param Xtest Test set
#' @param tune If TRUE, use grf automatic hyper-parameter tuning. If FALSE (default), no tuning.
#' @param num.trees Number of trees (default=500)
#' @param family Outcome type ("gaussian", "binomial"), default is "gaussian"
#' @param mod.A Model for estimating P(A|X). Default is "mean" calculates the sample mean.
#' If mod.A="RF", estimate P(A|X) using regression_forest (applicable for non-RCTs).
#' @param ... Any additional parameters, not currently passed through.
#'
#'
#' @return Trained causal_forest and regression_forest models.
#'  \itemize{
#'   \item mod - trained model(s)
#'   \item pred.fun - Prediction function for trained model(s)
#' }
#' @references Athey S, Tibshirani J, Wagner S. Generalized Random Forests. 
#' \url{https://arxiv.org/abs/1610.01271}
#'
#'
#' @export

#### Causal_forest ###
ple_causal_forest = function(Y, A, X, Xtest, tune=FALSE, num.trees=500, family="gaussian",
                             mod.A = "mean", ...){

  if (is.null(A)){
    stop("ple_causal_forest not applicable for no treatment (A=NULL)")
  }
  if (!requireNamespace("grf", quietly = TRUE)) {
    stop("Package grf needed for ple_causal_forest. Please install.")
  }
  A <- model.matrix(~., data=data.frame(A))[,-1]
  ## Regression Forest: Y~X ##
  forest.Y = grf::regression_forest(X, Y, ci.group.size=1,
                                    num.trees = min(500, num.trees) )
  Y.hat.train = predict(forest.Y)$predictions
  ## Regression Forest: W~X, If RCT ==> W is independent of X; use sample mean ##
  if (mod.A=="mean") {
    forest.A = mean(A)
    A.hat.train = rep( forest.A, length(Y))
  }
  if (mod.A=="propensity"){
    forest.A = grf::regression_forest(X, A, ci.group.size=1, 
                                      num.trees = min(500, num.trees) )
  }
  ## Causal Forest ##
  forest.CF = grf::causal_forest(X, Y, W=A, tune.parameters = tune, num.trees=num.trees,
                        W.hat=A.hat.train, Y.hat=Y.hat.train)
  mod = list(forest.Y=forest.Y, forest.A=forest.A, forest.CF=forest.CF)
  # Prediction Function #
  pred.fun = function(mod, X){
    forest.Y = mod$forest.Y
    forest.A = mod$forest.A
    forest.CF = mod$forest.CF
    Y_hat = predict(forest.Y, X)$predictions
    ## Regression Forest: W~X, If RCT ==> W is independent of X; use sample mean ##
    if (is.numeric(forest.A)){
        A_hat = rep( forest.A, length(Y_hat))
    }
    if (is.list(forest.A)){
        A_hat = predict(forest.A, X)
    }
    PLE_hat = predict(forest.CF, X)$predictions
    mu0_hat = Y_hat - A_hat * PLE_hat
    mu1_hat = Y_hat + (1 - A_hat) * PLE_hat
    mu_hat = data.frame(mu1 = mu1_hat, mu0 = mu0_hat, PLE = PLE_hat)
    return(mu_hat)
  }
  res = list(mod = list(forest.Y=forest.Y, forest.A=forest.A, forest.CF=forest.CF),
             pred.fun = pred.fun)
  class(res) = "ple_causal_forest"
  ## Return Results ##
  return( res )
}