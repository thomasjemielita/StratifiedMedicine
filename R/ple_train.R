#' Patient-level Estimates: Train Model
#'
#' Wrapper function to train a patient-level estimate (ple) model. Used directly in PRISM
#' and can be used to directly fit a ple model by name.
#'
#' @inheritParams PRISM
#' @param propensity Propensity score estimation, P(A=a|X). Default=NULL which use the 
#' marginal estimates, P(A=a) (applicable for RCT data). 
#' @param hyper Hyper-parameters for the ple model (must be list). Default is NULL.
#' @param tau Maximum follow-up time for RMST based estimates (family="survival"). 
#' Default=NULL, which takes min(max(time[a])), for a=1,..,A.
#' @param ... Any additional parameters, not currently passed through.
#'
#'
#' @return Trained ple models and patient-level estimates for train/test sets. 
#'  \itemize{
#'   \item mod - trained model(s)
#'   \item mu_train - Patient-level estimates (training set)
#'   \item mu_test - Patient-level estimates (test set)
#' }
#' 
#' @details ple_train uses base-learners along with a meta-learner to obtain patient-level 
#' estimates under different treatment exposures (see Kunzel et al).  For family="gaussian" 
#' or "binomial", output estimates of \eqn{mu(a,x)=E(Y|x,a)} and treatment differences 
#' (average treatment effect or risk difference). For survival, either logHR based estimates
#' or RMST based estimates can be obtained. Current base-learner ("ple") options include:
#' 
#' 
#' 1. \strong{linear}: Uses either linear regression (family="gaussian"), 
#' logistic regression (family="binomial"), or cox regression (family="survival"). 
#' No hyper-parameters.
#' 
#' 2. \strong{ranger}: Uses random forest ("ranger" R package). The default hyper-parameters are: 
#' hyper = list(mtry=NULL, min.node.pct=0.10)
#' 
#' where mtry is number of randomly selected variables (default=NULL; sqrt(dim(X)))
#' and min.node.pct is the minimum node size as a function of the total data size 
#' (ex: min.node.pct=10\% requires at least 10% of the data at each node)
#' 
#' 3. \strong{glmnet}: Uses elastic net ("glmnet" R package). The default hyper-parameters are: 
#' hyper = list(lambda="lambda.min")
#' 
#' where lambda controls the penalty parameter for predictions. lambda="lambda.1se"
#' will likely result in a less complex model. 
#' 
#' 4. \strong{bart}:  Uses bayesian additive regression trees (Chipman et al 2010; 
#' BART R package). Default hyper-parameters are:
#' 
#' hyper = list(sparse=FALSE)
#' 
#' where sparse controls whether to perform variable selection based on a sparse 
#' Dirichlet prior rather than simply uniform.
#' 
#' @examples
#' \donttest{
#' library(StratifiedMedicine)
#' ## Continuous ##
#' dat_ctns = generate_subgrp_data(family="gaussian")
#' Y = dat_ctns$Y
#' X = dat_ctns$X
#' A = dat_ctns$A
#'
#'
#' # X-Learner (With ranger based learners)
#' mod1 = ple_train(Y=Y, A=A, X=X, Xtest=X, ple="ranger", method="X-learner")
#' summary(mod1$mu_train)
#' 
#' # T-Learner (Treatment specific)
#' mod2 = ple_train(Y=Y, A=A, X=X, Xtest=X, ple="ranger", method="T-learner")
#' summary(mod2$mu_train)
#' 
#' 
#' mod3 = ple_train(Y=Y, A=A, X=X, Xtest=X, ple="bart", method="X-learner")
#' summary(mod3$mu_train)
#' }
#'
#'
#' @export
#' @references Wright, M. N. & Ziegler, A. (2017). ranger: A fast implementation of 
#' random forests for high dimensional data in C++ and R. J Stat Softw 77:1-17. 
#' \url{https://doi.org/10.18637/jss.v077.i01}.
#' @references Friedman, J., Hastie, T. and Tibshirani, R. (2008) Regularization Paths for
#'  Generalized Linear Models via Coordinate Descent,
#'  \url{https://web.stanford.edu/~hastie/Papers/glmnet.pdf} Journal of Statistical 
#'  Software, Vol. 33(1), 1-22 Feb 2010 Vol. 33(1), 1-22 Feb 2010.
#' @references Chipman, H., George, E., and McCulloch R. (2010) Bayesian Additive 
#' Regression Trees. The Annals of Applied Statistics, 4,1, 266-298 
#' @references Kunzel S, Sekhon JS, Bickel PJ, Yu B. Meta-learners for Estimating
#' Hetergeneous Treatment Effects using Machine Learning. 2019. 
#' @seealso \code{\link{PRISM}}
#'
## To DO: T-Learner, S-Learner, S-Learner (VT) ##
ple_train = function(Y, A, X, Xtest=NULL, family="gaussian", propensity=NULL,
                      ple="ranger", meta="X-learner", hyper=NULL,
                      tau=NULL, ...) {
  
  if (is.null(Xtest)) {
    Xtest <- X
  }
  # Convert Names #
  if (ple %in% c("ranger", "glmnet", "linear", "bart")) {
    ple <- paste("ple", ple, sep="_") 
  }
  
  if (is.null(A) | !(meta %in% c("T-learner", "S-learner", "X-learner"))) {
    fit <- do.call(ple, append(list(Y=Y, X=X, family=family), hyper))
    names(fit) <- c("mod", "pred.fun")
  }
  if (!is.null(A)) {
    A_lvls <- unique(A)[order(unique(A))]
    # Propensity Estimation: RCT (sample average; or model-based) #
    if (is.null(propensity)) {
      pi_fit <- NULL
      for (aa in A_lvls) {
        pi.a <- mean(A==aa)
        pi_fit <- cbind(pi_fit, pi.a)
      }
      colnames(pi_fit) <- paste("pi", A_lvls, sep="_")
      pi_fit <- data.frame(pi_fit)
    }
    if (!is.null(propensity)) {
      ## TO DO ##
    }
    if (is.Surv(Y)) {
      tau <- NULL
      for (a in unique(A)) {
        tau.a <- max(Y[,1][A==a])
        tau <- c(tau, tau.a)
      }
      tau <- min(tau)
    }
    if (!is.Surv(Y)) {
      tau <- NULL
    }
    if (meta=="T-learner") {
      fit <- T_learner(Y=Y, A=A, X=X, family=family, ple=ple, hyper=hyper,
                       tau=tau, pi_fit=pi_fit) 
    }
    if (meta=="S-learner") {
      fit <- S_learner(Y=Y, A=A, X=X, family=family, ple=ple, hyper=hyper,
                       tau=tau, pi_fit=pi_fit) 
    }
    if (meta=="X-learner") {
      fit <- X_learner(Y=Y, A=A, X=X, family=family, ple=ple, hyper=hyper, 
                       tau=tau, pi_fit=pi_fit) 
    }
  }
  ### Train/Test Predictions ###
  ## If prior predictions are made: ##
  if (!is.null(fit$mu_train)) {
    mu_train <- fit$mu_train
  }
  if (!is.null(fit$mu_test)) {
    mu_test <- fit$mu_test
  }
  ## If no prior predictions are mode: ##
  if (is.null(fit$mu_train)) {
    mu_train <- fit$pred.fun(fit$mod, X=X, tau=tau)
  }
  if (is.null(fit$mu_test)) {
    mu_test <- fit$pred.fun(fit$mod, X=Xtest, tau=tau)
  }
  res <- list(fit = fit, mu_train=mu_train, mu_test=mu_test)
  class(res) <- "ple_train"
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
#' \donttest{
#' library(StratifiedMedicine)
#' ## Continuous ##
#' dat_ctns = generate_subgrp_data(family="gaussian")
#' Y = dat_ctns$Y
#' X = dat_ctns$X
#' A = dat_ctns$A
#'
#'
#' mod1 = ple_train(Y=Y, A=A, X=X, Xtest=X, ple="ranger", meta="X-learner")
#' summary(mod1$mu_train)
#'
#' res1 = predict(mod1, newdata=X)
#' summary(res1)
#' }
#'
#' @method predict ple_train
#' @export
#' @seealso \code{\link{PRISM}}
#'
predict.ple_train = function(object, newdata=NULL, ...) {
  
  mu_hat <- object$fit$pred.fun(object$fit$mod, newdata)
  
  return(mu_hat)
}