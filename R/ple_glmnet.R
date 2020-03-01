#' Patient-level Estimates: Elastic Net (glmnet)
#'
#' Uses the elastic net (glmnet R package) to obtain patient-level estimates. Usable for
#' continuous, binary, or survival outcomes.
#'
#' @param Y The outcome variable. Must be numeric or survival (ex; Surv(time,cens) )
#' @param A Treatment variable. (a=1,...A)
#' @param X Covariate space.
#' @param Xtest Test set
#' @param lambda Lambda for elastic-net (default = "lambda.min"). Other options include
#' "lambda.1se" or  fixed values
#' @param family Outcome type ("gaussian", "binomial", "survival"), default is "gaussian"
#' @param ... Any additional parameters, not currently passed through.
#'
#' @import glmnet
#'
#' @return Trained glmnet model(s).
#'  \itemize{
#'   \item mod - trained model(s)
#'   \item lambda - Lambda used for elastic-net (passes to prediction function)
#'   \item X - Covariate Space (in model matrix form)
#' }
#' @export
#' @references Friedman, J., Hastie, T. and Tibshirani, R. (2008) Regularization Paths for
#'  Generalized Linear Models via Coordinate Descent,
#'  \url{https://web.stanford.edu/~hastie/Papers/glmnet.pdf} Journal of Statistical 
#'  Software, Vol. 33(1), 1-22 Feb 2010 Vol. 33(1), 1-22 Feb 2010.
#' @examples
#' library(StratifiedMedicine)
#'
#' ## Continuous ##
#' dat_ctns = generate_subgrp_data(family="gaussian")
#' Y = dat_ctns$Y
#' X = dat_ctns$X
#' A = dat_ctns$A
#'
#' mod1 = ple_glmnet(Y, A, X, Xtest=X, family="gaussian")
#'

##### Elastic net (glmnet): Y~(A,X,A*X) ==> PLEs ######
ple_glmnet = function(Y, A, X, Xtest, lambda="lambda.min", family, ...) {

  ## Model matrix the covariate space ###
  X = model.matrix(~., data = X)[,-1]
  if (is.null(A)) {
    W = X
  }
  if (!is.null(A)) {
    A_lvls <- unique(A)[order(unique(A))]
    ## Generate the interaction covariates ###
    A.mat <- model.matrix(~., data=data.frame(A))[,-1]
    X_inter = X*A.mat
    colnames(X_inter) = paste(colnames(X), A_lvls[2], sep="_")
    W = cbind(A=A.mat, X, X_inter)
  }
  ##### Elastic Net #####
  if (family=="survival") { family = "cox"  }
  mod <- cv.glmnet(x = W, y = Y, alpha=0.5, family=family)
  mod$sel.lambda <- lambda
  mod$A.length <- length(unique(A))
  mod$A_lvls <- A_lvls
  pred.fun = function(mod, X) {
    lambda <- mod$sel.lambda
    A.length <- mod$A.length
    A_lvls <- mod$A_lvls
    X = model.matrix(~., data = X)[,-1]
    ### Predictions (Counterfactuals) ###
    if (A.length==0) {
      mu_hat = data.frame(mu=as.numeric(predict(mod, newx = X, s=lambda))) 
    }
    if (A.length>0) {
      X0 = data.frame(0, X, X*0)
      colnames(X0) = c("A", colnames(X), paste(colnames(X), A_lvls[2], sep="_"))
      X1 = data.frame(1, X, X*1)
      colnames(X1) = c( "A", colnames(X), paste(colnames(X), A_lvls[2], sep="_"))
      X1 <- as.matrix(X1); X0 <- as.matrix(X0)
      mu_hat = data.frame(mu0 = predict(mod, newx=X0, s=lambda),
                          mu1 = predict(mod, newx=X1, s=lambda))
      if (!is.null(A_lvls)) colnames(mu_hat) <- paste("mu", A_lvls, sep="_")
      mu_hat$PLE = mu_hat[,2] - mu_hat[,1]
    }
    return(mu_hat)
  }
  res = list(mod = mod, pred.fun=pred.fun)
  class(res) = "ple_glmnet"
  ## Return Results ##
  return( res )
}