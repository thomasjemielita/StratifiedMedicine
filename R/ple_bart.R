#' Patient-level Estimates: BART
#'
#' Uses the BART algorithm (Chipman et al 2010; BART R package) to obtain patient-level
#' estimates. Used for continuous or binary outcomes. Covariate by treatment interactions
#' are automatically included in BART model (as in Hahn et al 2017).
#'
#' @param Y The outcome variable. Must be numeric or survival (ex; Surv(time,cens) )
#' @param A Treatment variable. (a=1,...A)
#' @param X Covariate space.
#' @param Xtest Test set
#' @param family Outcome type ("gaussian", "binomial"), default is "gaussian"
#' @param ... Any additional parameters, not currently passed through.
#'
#'
#' @return Trained BART model(s) and patient-level estimates
#' (E(Y|X,1), E(Y|X,0), E(Y|X,1)-E(Y|X,0)) for train/test sets.
#'  \itemize{
#'   \item mod - trained model(s)
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
#' require(BART)
#' mod1 = ple_bart(Y, A, X, Xtest=X)
#' summary(mod1$mu_train)
#' summary(predict(mod1, newdata=X))
#' }
#'
#'
### BART ###
ple_bart = function(Y, A, X, Xtest, family="gaussian", ...){

  if (!requireNamespace("BART", quietly = TRUE)) {
    stop("Package BART needed for ple_bart. Please install.")
  }
  if (is.null(A)){
    mod = BART::wbart(x.train = X, y.train = Y, x.test = Xtest)
    mu_train = data.frame(PLE = mod$yhat.train.mean)
    mu_test = data.frame(PLE = mod$yhat.test.mean)
    pred.fun = function(mod, X){
      hold = predict(mod, newdata=X)
      mu_hat = data.frame( PLE = apply(hold, 2, mean)  )
      return( mu_hat )
    }
  }
  if (!is.null(A)){
    ## Generate the covariate by treatment interactions ##
    X_inter = X*A
    colnames(X_inter) = paste(colnames(X), "_A", sep="")
    W = cbind(X, A, X_inter)
    ## Generate counterfactual design matrices  ##
    Xtrain_0 = data.frame(0, X, X*0)
    colnames(Xtrain_0) = c( "A", colnames(X), paste(colnames(X), "_A", sep="") )
    Xtrain_1 = data.frame(1, X, X*1)
    colnames(Xtrain_1) = c( "A", colnames(X), paste(colnames(X), "_A", sep="") )
    Xtest_0 = data.frame(0, Xtest, Xtest*0)
    colnames(Xtest_0) = c( "A", colnames(Xtest), paste(colnames(Xtest), "_A", sep="") )
    Xtest_1 = data.frame(1, Xtest, Xtest*1)
    colnames(Xtest_1) = c( "A", colnames(Xtest), paste(colnames(Xtest), "_A", sep="") )
    # Test: Include both train/test for faster predictions #
    X.FULL = rbind(Xtrain_0, Xtrain_1, Xtest_0, Xtest_1)
    ## BART ##
    mod = BART::wbart(x.train = W, y.train = Y, x.test = X.FULL)
    n.tr = dim(Xtrain_0)[1]
    n.ts = dim(Xtest_0)[1]
    ### PLE Predictions: Train/Test ###
    mu_train = data.frame(mu1 = mod$yhat.test.mean[(n.tr+1):(2*n.tr)],
                          mu0 = mod$yhat.test.mean[1:n.tr])
    mu_train$PLE = mu_train$mu1-mu_train$mu0
    mu_test = data.frame(mu1 = mod$yhat.test.mean[(2*n.tr+n.ts+1):(2*n.tr+2*n.ts)],
                         mu0 = mod$yhat.test.mean[(2*n.tr+1):(2*n.tr+n.ts)] )
    mu_test$PLE = mu_test$mu1 - mu_test$mu0
    pred.fun = function(mod, X){
      ## Generate counterfactual design matrices  ##
      X0 = data.frame(0, X, X*0)
      colnames(X0) = c( "A", colnames(X), paste(colnames(X), "_A", sep="") )
      X1 = data.frame(1, X, X*1)
      colnames(X1) = c( "A", colnames(X), paste(colnames(X), "_A", sep="") )
      X.FULL = rbind(X0, X1)
      preds = predict(mod, newdata=X.FULL)
      preds = apply(preds, 2, mean)
      n = dim(X)[1]
      ### PLE Predictions: Train/Test ###
      mu_hat = data.frame(mu1 = preds[(n+1):(2*n)], mu0 = preds[1:n])
      mu_hat$PLE = with(mu_hat, mu1-mu0)
      return( mu_hat )
    }
  }
  res = list(mod = mod, pred.fun=pred.fun, mu_train = mu_train, mu_test=mu_test, A=A)
  class(res) = "ple_bart"
  ## Return Results ##
  return( res )
}