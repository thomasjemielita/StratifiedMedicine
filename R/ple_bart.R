#' Patient-level Estimates: BART
#'
#' Uses the BART algorithm (Chipman et al 2010; BART R package) to obtain patient-level
#' estimates. Used for continuous or binary outcomes. Covariate by treatment interactions
#' are automatically included in BART model (as in Hahn et al 2017).
#'
#' @param Y The outcome variable. Must be numeric. 
#' @param A Treatment variable. (a=1,...A)
#' @param X Covariate space.
#' @param Xtest Test set
#' @param family Outcome type ("gaussian", "binomial"), default is "gaussian"
#' which uses wbart (BART for continuous outcomes). Probit-based BART ("pbart") is
#' used for family="binomial"
#' @param sparse Whether to perform variable selection based on sparse Dirichlet prior
#' instead of uniform. See "gbart" in the BART R package for more details as well as
#' Linero 2016. Default is FALSE. 
#' @param K For survival, coarse times per the quantiles 1/K,..,K/K. For AFT (abart),
#' K=10. For more accurate survival predictions, set K=100 (default in abart).
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
#' @references Chipman, H., George, E., and McCulloch R. (2010) Bayesian Additive 
#' Regression Trees. The Annals of Applied Statistics, 4,1, 266-298 
#' <doi: 10.1214/09-AOAS285>.
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
ple_bart = function(Y, A, X, Xtest, family="gaussian", sparse=FALSE,
                    K=10,...){

  if (!requireNamespace("BART", quietly = TRUE)) {
    stop("Package BART needed for ple_bart. Please install.")
  }
  if (is.null(A)){
    if (family=="gaussian"){
      mod <- BART::gbart(x.train = X, y.train = Y, x.test = Xtest,
                        type = "wbart", sparse = sparse)
      mu_train <- data.frame(PLE = mod$yhat.train.mean)
      mu_test <- data.frame(PLE = mod$yhat.test.mean)
      pred.fun <- function(mod, X){
        hold <- predict(mod, newdata=X)
        mu_hat <- data.frame( PLE = apply(hold, 2, mean)  )
        return( mu_hat )
      }
    }
    if (family=="binomial"){
      mod <- BART::gbart(x.train = X, y.train = Y, x.test = Xtest,
                        type = "pbart", sparse = sparse)
      mu_train <- data.frame(PLE = mod$prob.train.mean)
      mu_test <- data.frame(PLE = mod$prob.test.mean)
      pred.fun <- function(mod, X){
        hold <- predict(mod, newdata=X)
        mu_hat <- data.frame( PLE = hold$prob.test.mean  )
        return( mu_hat )
      }
    }
    if (family=="survival"){
      mod <- BART::abart(x.train = X, times = Y[,1], delta = Y[,2],
                        x.test = Xtest, sparse = sparse, K=K)
      mu_train <- data.frame(PLE = mod$yhat.train.mean)
      mu_test <- data.frame(PLE = mod$yhat.test.mean)
      pred.fun <- "prediction function not available for survival (abart)"
    }
  }
  if (!is.null(A)){
    Xmat<- model.matrix(~., data=X)[,-1]
    A.mat <- model.matrix(~., data=data.frame(A))[,-1]
    X_inter <- Xmat * A.mat
    colnames(X_inter) <- paste(colnames(Xmat), "_A", sep="")
    W <- cbind(Xmat, X_inter)
    ## Generate the covariate by treatment interactions ##
    gen_inters <- function(X){
      Xmat <- model.matrix(~., data=X)[,-1]
      ## Generate counterfactual design matrices  ##
      X0 <- data.frame(Xmat, Xmat*0)
      colnames(X0) <- colnames(W)
      X1 <- data.frame(Xmat, Xmat*1)
      colnames(X1) <- colnames(W)
      # x.test: Include both train/ for faster predictions #
      X.FULL <- rbind(X0, X1)
      return( list(W=W, X.FULL=X.FULL) )
    }
    X.tr <- gen_inters(X)
    X.ts <- gen_inters(Xtest)
    X.FULL <- rbind( X.tr$X.FULL, X.ts$X.FULL )
    ## BART Fit ##
    if (family=="gaussian"){
      mod <- BART::gbart(x.train = W, y.train = Y, x.test = X.FULL, 
                        type = "wbart", sparse=sparse)
      n.tr <- dim(X)[1]
      n.ts <- dim(Xtest)[1]
      ### PLE Predictions: Train/Test ###
      mu_train <- data.frame(mu1 = mod$yhat.test.mean[(n.tr+1):(2*n.tr)],
                            mu0 = mod$yhat.test.mean[1:n.tr])
      mu_train$PLE <- mu_train$mu1-mu_train$mu0
      mu_test <- data.frame(mu1 = mod$yhat.test.mean[(2*n.tr+n.ts+1):(2*n.tr+2*n.ts)],
                           mu0 = mod$yhat.test.mean[(2*n.tr+1):(2*n.tr+n.ts)] )
      mu_test$PLE <- mu_test$mu1 - mu_test$mu0
      pred.fun <- function(mod, X){
        ## Generate the covariate by treatment interactions ##
        Xmat <- model.matrix(~., data=X)[,-1]
        cnames = c(colnames(Xmat),  paste(colnames(Xmat), "_A", sep="") )
        ## Generate counterfactual design matrices  ##
        X0 <- data.frame(Xmat, Xmat*0)
        colnames(X0) <- cnames
        X1 <- data.frame(Xmat, Xmat*1)
        colnames(X1) <- cnames
        X.FULL <- rbind(X0, X1)
        preds <- predict(mod, newdata=X.FULL)
        preds <- apply(preds, 2, mean)
        n <- dim(X)[1]
        ### PLE Predictions: Train/Test ###
        mu_hat <- data.frame(mu1 = preds[(n+1):(2*n)], mu0 = preds[1:n])
        mu_hat$PLE <- with(mu_hat, mu1-mu0)
        return( mu_hat )
      }
    }
    if (family=="binomial"){
      mod <- BART::gbart(x.train = W, y.train = Y, x.test = X.FULL, 
                        type = "pbart", sparse = sparse)
      n.tr <- dim(X)[1]
      n.ts <- dim(Xtest)[1]
      ### PLE Predictions: Train/Test ###
      mu_train <- data.frame(mu1 = mod$prob.test.mean[(n.tr+1):(2*n.tr)],
                            mu0 = mod$prob.test.mean[1:n.tr])
      mu_train$PLE <- mu_train$mu1-mu_train$mu0
      mu_test <- data.frame(mu1 = mod$prob.test.mean[(2*n.tr+n.ts+1):(2*n.tr+2*n.ts)],
                           mu0 = mod$prob.test.mean[(2*n.tr+1):(2*n.tr+n.ts)] )
      mu_test$PLE = mu_test$mu1 - mu_test$mu0
      pred.fun = function(mod, X){
        ## Generate the covariate by treatment interactions ##
        Xmat <- model.matrix(~., data=X)[,-1]
        cnames = c(colnames(Xmat),  paste(colnames(Xmat), "_A", sep="") )
        ## Generate counterfactual design matrices  ##
        X0 <- data.frame(Xmat, Xmat*0)
        colnames(X0) <- cnames
        X1 <- data.frame(Xmat, Xmat*1)
        colnames(X1) <- cnames
        X.FULL <- rbind(X0, X1)
        preds <- predict(mod, newdata=X.FULL)
        preds <- preds$prob.test.mean
        n <- dim(X)[1]
        ### PLE Predictions: Train/Test ###
        mu_hat <- data.frame(mu1 = preds[(n+1):(2*n)], mu0 = preds[1:n])
        mu_hat$PLE <- with(mu_hat, mu1-mu0)
        return( mu_hat )
      }
    }
    if (family=="survival"){
      mod <- BART::abart(x.train = W, times = Y[,1], delta = Y[,2],
                        x.test = X.FULL, sparse = sparse, K=K)
      n.tr <- dim(X)[1]
      n.ts <- dim(Xtest)[1]
      ### PLE Predictions: Train/Test ###
      mu_train <- data.frame(mu1 = mod$yhat.test.mean[(n.tr+1):(2*n.tr)],
                            mu0 = mod$yhat.test.mean[1:n.tr])
      mu_train$PLE <- mu_train$mu1-mu_train$mu0
      mu_test <- data.frame(mu1 = mod$yhat.test.mean[(2*n.tr+n.ts+1):(2*n.tr+2*n.ts)],
                           mu0 = mod$yhat.test.mean[(2*n.tr+1):(2*n.tr+n.ts)] )
      mu_test$PLE <- mu_test$mu1 - mu_test$mu0
      pred.fun <- "prediction function not available for survival (abart)"
    }
    
  }
  res <- list(mod = mod, pred.fun=pred.fun, mu_train = mu_train, mu_test=mu_test)
  class(res) <- "ple_bart"
  ## Return Results ##
  return( res )
}