#' Subgroup Identification: Train Model
#'
#' Wrapper function to train a subgroup model (submod). Outputs subgroup assignments and 
#' fitted model.
#'
#' @inheritParams PRISM
#' @param mu_train Patient-level estimates in training set (see \code{ple_train}). 
#' Default=NULL
#' @param hyper Hyper-parameters for submod (must be list). Default is NULL.
#' @param ... Any additional parameters, not currently passed through.
#'
#'
#' @return Trained subgroup model and subgroup predictions/estimates for train/test sets.
#'
#'  \itemize{
#'   \item mod - trained subgroup model
#'   \item Subgrps.train - Identified subgroups (training set)
#'   \item Subgrps.test - Identified subgroups (test set)
#'   \item pred.train - Predictions (training set)
#'   \item pred.test - Predictions (test set)
#'   \item Rules - Definitions for subgroups, if provided in fitted submod output.
#' }
#' @details submod_train currently fits a number of tree-based subgroup models, most of
#' which aim to find subgroups with varying treatment effects (i.e. predictive variables).
#' Current options include:
#' 
#' 1. lmtree: Wrapper function for the function "lmtree" from the partykit package. Here, 
#' model-based partitioning (MOB) with an OLS loss function, Y~MOB_LM(A,X), is used to 
#' identify prognostic and/or predictive variables. 
#' 
#' Default hyper-parameters are: 
#' hyper = list(alpha=0.05, maxdepth=4, parm=NULL, minsize=floor(dim(X)[1]*0.10)).
#' 
#' 2. glmtree: Wrapper function for the function "glmtree" from the partykit package. Here, 
#' model-based partitioning (MOB) with GLM binomial + identity link loss function, 
#' (Y~MOB_GLM(A,X)), is used to identify prognostic and/or predictive variables.
#' 
#' Default hyper-parameters are:
#'  hyper = list(link="identity", alpha=0.05, maxdepth=4, parm=NULL, 
#'  minsize=floor(dim(X)[1]*0.10)).
#' 
#' 3. ctree: Wrapper function for the function "ctree" from the partykit package. Here, 
#' conditional inference trees are used to identify either prognostic, Y~CTREE(X), 
#' or predictive variables, PLE~CTREE(X) (outcome_PLE=TRUE; requires mu_train data).
#' 
#' Default hyper-parameters are:
#' hyper=list(alpha=0.10, minbucket = floor(dim(X)[1]*0.10), 
#' maxdepth = 4, outcome_PLE=FALSE). 
#' 
#' 4. otr: Optimal treatment regime approach using "ctree". Based on patient-level 
#' treatment effect estimates, fit PLE~CTREE(X) with weights=abs(PLE). 
#' 
#' Default hyper-parameters are:
#' hyper=list(alpha=0.10, minbucket = floor(dim(X)[1]*0.10), 
#' maxdepth = 4, thres=">0"). 
#' 
#' 4. mob_weib: Wrapper function for the function "mob" with weibull loss function using
#' the partykit package. Here, model-based partitioning (MOB) with weibull loss (survival),
#' (Y~MOB_WEIB(A,X)), is used to identify prognostic and/or predictive variables.
#'  
#' Default hyper-parameters are:
#' hyper = list(alpha=0.10, maxdepth=4, parm=NULL, minsize=floor(dim(X)[1]*0.10)).
#' 
#' 5. rpart: Recursive partitioning through the "rpart" R package. Here, 
#' recursive partitioning and regression trees are used to identify either prognostic,
#' Y~rpart(X), or predictive variables, PLE~rpart(X) (outcome_PLE=TRUE; 
#' requires mu_train data).
#' 
#' 
#' @references
#' \itemize{
#' \item Zeileis A, Hothorn T, Hornik K (2008). Model-Based Recursive Partitioning. 
#' Journal of Computational and Graphical Statistics, 17(2), 492–514.
#' \item Seibold H, Zeileis A, Hothorn T. Model-based recursive partitioning for 
#' subgroup analyses. Int J Biostat, 12 (2016), pp. 45-63
#' \item Hothorn T, Hornik K, Zeileis A (2006). Unbiased Recursive Partitioning: 
#' A Conditional Inference Framework. Journal of Computational and Graphical Statistics,
#' 15(3), 651–674.
#' \item Zhao et al. (2012) Estimated individualized treatment rules using outcome 
#' weighted learning. Journal of the American Statistical Association, 107(409): 1106-1118.
#' \item Breiman L, Friedman JH, Olshen RA, and Stone CJ. (1984) Classification 
#' and Regression Trees. Wadsworth
#' } 
#' @examples
#' 
#' \donttest{
#' library(StratifiedMedicine)
#' ## Continuous ##
#' dat_ctns = generate_subgrp_data(family="gaussian")
#' Y = dat_ctns$Y
#' X = dat_ctns$X
#' A = dat_ctns$A
#'
#' # Fit through submod_train wrapper #
#' mod1 = submod_train(Y=Y, A=A, X=X, Xtest=X, submod="submod_lmtree")
#' table(mod1$Subgrps.train)
#' plot(mod1$fit$mod)
#'
#'}
#'
#' @export
#' @importFrom partykit lmtree ctree glmtree as.party
#' @importFrom stats binomial
#' @seealso \code{\link{PRISM}}
#'
submod_train = function(Y, A, X, Xtest=NULL, mu_train=NULL, 
                        family="gaussian", submod, hyper=NULL, 
                        pool="no", delta=">0", ...){
  if (is.null(Xtest)) {
    Xtest <- X
  }
  # Convert Names #
  if (submod %in% c("lmtree", "ctree", "glmtree", "otr", "rpart", "mob_weib")) {
    submod <- paste("submod", submod, sep="_") 
  }
  ## Fit submod ##
  fit = do.call(submod, append(list(Y=Y, A=A, X=X, Xtest=Xtest, mu_train=mu_train,
                                  family=family), hyper))
  ### Train/Test Predictions ###
  ## If prior predictions are made: ##
  if (!is.null(fit$Subgrps.train)){
    Subgrps.train = fit$Subgrps.train
    pred.train = fit$pred.train
  }
  if (!is.null(fit$Subgrps.test)){
    Subgrps.test = fit$Subgrps.test
    pred.test = fit$pred.test
  }
  ## If no prior predictions are made: ##
  if (is.null(fit$Subgrps.train)){
    out = fit$pred.fun(fit$mod, X=X)
    Subgrps.train = out$Subgrps
    pred.train = out$pred
  }
  if (is.null(fit$Subgrps.test)){
    out = fit$pred.fun(fit$mod, X=Xtest)
    Subgrps.test = out$Subgrps
    pred.test = out$pred
  }
 
  Rules = fit$Rules
  
  ## Pooling? ##
  pool.dat <- NULL
  if (pool %in% c("otr:logistic", "otr:rf")){
    pool.dat <- pooler_run(Y, A, X, mu_hat=mu_train, Subgrps=Subgrps.train, 
                            delta = delta, method = pool)
    # Merge with training #
    Subgrps.train0 <- as.character(Subgrps.train)
    subdat <- data.frame(Subgrps=Subgrps.train0)
    subdat <- suppressWarnings(left_join(subdat, pool.dat, by="Subgrps"))
    Subgrps.train <- as.character(subdat$pred_opt)
    # Merge with test #
    Subgrps.test0 <- as.character(Subgrps.test)
    subdat <- data.frame(Subgrps=Subgrps.test0)
    subdat <- suppressWarnings(left_join(subdat, pool.dat, by="Subgrps"))
    Subgrps.test <- as.character(subdat$pred_opt)
  }
  
  res = list(fit = fit, Subgrps.train=Subgrps.train, Subgrps.test=Subgrps.test, 
             pool.dat=pool.dat, Rules=Rules)
  class(res) = "submod_train"
  return(res)
}

#' Subgroup Identification: Train Model (Predictions)
#'
#' Prediction function for the trained subgroup identification model (submod).
#'
#' @param object Trained submod model.
#' @param newdata Data-set to make predictions at (Default=NULL, predictions correspond
#' to training data).
#' @param ... Any additional parameters, not currently passed through.
#'
#' @return Identified subgroups with subgroup-specific predictions (depends on subgroup
#' model)
#' \itemize{
#'   \item Subgrps - Identified subgroups
#'   \item pred - Predictions, depends on subgroup model
#'}
#' @examples
#' library(StratifiedMedicine)
#' ## Continuous ##
#' dat_ctns = generate_subgrp_data(family="gaussian")
#' Y = dat_ctns$Y
#' X = dat_ctns$X
#' A = dat_ctns$A
#'
#' # Fit through submod_train wrapper #
#' mod1 = submod_train(Y=Y, A=A, X=X, Xtest=X, submod="submod_lmtree")
#' out1 = predict(mod1)
#' table(mod1$Subgrps.train)
#' table(out1$Subgrps)
#'
#' @method predict submod_train
#' @export
#'
predict.submod_train = function(object, newdata=NULL, ...){

  # preds = predict(object$fit, newdata=newdata)
  preds = object$fit$pred.fun(object$fit$mod, newdata)
  ## Return Results ##
  return(  list(Subgrps=preds$Subgrps, pred=preds$pred) )
}
