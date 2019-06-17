#' Patient-level Estimates: Elastic Net (glmnet)
#'
#' Uses the elastic net (glmnet R package) to obtain patient-level estimates.
#' For continuous/binary, output estimates of E(Y|X,A=a) and E(Y|X,A=1)-E(Y|X,A=0) (PLE).
#' For survival, output estimates of HR(X,A=a) and HR(X, A=1)-HR(X, A=0) (PLE).
#'
#' @param Y The outcome variable. Must be numeric or survival (ex; Surv(time,cens) )
#' @param A Treatment variable. (a=1,...A)
#' @param X Covariate matrix. Must be numeric.
#' @param Xtest Test set
#' @param lambda Lambda for elastic-net (default = "lambda.min"). Other options include
#' "lambda.1se" or  fixed values
#' @param family Outcome type ("gaussian", "binomial", "survival"), default is "gaussian"
#' @param ... Any additional parameters, not currently passed through.
#'
#' @import glmnet
#'
#' @return Patient-level estimates (E(Y|X,1), E(Y|X,0), E(Y|X,1)-E(Y|X,0)) or
#'  (HR(X,1), HR(X,0), HR(X,1)-HR(X,0)) for train/test sets.
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
#'
#' mod1 = ple_glmnet(Y, A, X, Xtest=X, family="gaussian")
#' summary(mod1$mu_train$PLE)
#'
#'
#' @seealso \code{\link{PRISM}}, \code{\link{glmnet}}

##### Elastic net (glmnet): Y~(A,X,A*X) ==> PLEs ######
ple_glmnet = function(Y, A, X, Xtest, lambda="lambda.min", family, ...){

  ## Extract covariate space ###
  X = model.matrix(~. -1, data = X )
  Xtest = model.matrix(~. -1, data = Xtest )
  ## Generate the interaction covariates ###
  X_inter = X*A
  colnames(X_inter) = paste(colnames(X), "_A", sep="")
  W = cbind(X, A, X_inter)
  ##### Elastic Net #####
  set.seed(6134)
  if (family=="survival") { family = "cox"  }
  mod <- cv.glmnet(x = W, y = Y, alpha=0.5, family=family)

  ### Predictions (Counterfactuals): Train/Test ##
  mu_train = data.frame( mu1 =  as.numeric(predict(mod, newx = cbind(X, 1, X*1), s=lambda)),
                         mu0 =  as.numeric(predict(mod, newx = cbind(X, 0, X*0), s=lambda)) )
  mu_train$PLE = with(mu_train, mu1 - mu0 )

  mu_test = data.frame( mu1 =  as.numeric(predict(mod, newx = cbind(Xtest, 1, Xtest*1),
                                                  s=lambda)),
                        mu0 =  as.numeric(predict(mod, newx = cbind(Xtest, 0, Xtest*0),
                                                  s=lambda)) )
  mu_test$PLE = with(mu_test, mu1 - mu0 )

  ## Return Results ##
  return( list(mods = mod, mu_train = mu_train, mu_test=mu_test) )
}
