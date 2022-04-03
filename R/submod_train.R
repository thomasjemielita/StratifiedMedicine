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
#' Let E(Y|A=1,X)-E(Y|A=0,X) = CATE(X) correspond to the estimated conditional average treatment 
#' effect. Current options include:
#' 
#' 1. lmtree: Wrapper function for the function "lmtree" from the partykit package. Here, 
#' model-based partitioning (MOB) with an OLS loss function, Y~MOB_OLS(A,X), is used to 
#' identify prognostic and/or predictive variables. If the outcome Y is survival, then 
#' this outcome will first be transformed via log-rank scores (coin::logrank_trafo(Y)).
#' 
#' Default hyper-parameters are: 
#' hyper = list(alpha=0.05, maxdepth=4, parm=NULL, minsize=floor(dim(X)[1]*0.10)).
#' 
#' 2. glmtree: Wrapper function for the function "glmtree" from the partykit package. Here, 
#' model-based partitioning (MOB) with GLM binomial + identity link loss function, 
#' (Y~MOB_GLM(A,X)), is used to identify prognostic and/or predictive variables.
#' 
#' Default hyper-parameters are:
#' hyper = list(link="identity", alpha=0.05, maxdepth=4, parm=NULL, minsize=floor(dim(X)[1]*0.10)).
#' 
#' 3. ctree / ctree_cate: Wrapper function for the function "ctree" from the partykit package. Here, 
#' conditional inference trees are used to identify either prognostic ("ctree"), Y~CTREE(X), 
#' or predictive variables, CATE(X) ~ CTREE(X).
#' 
#' Default hyper-parameters are:
#' hyper=list(alpha=0.10, minbucket = floor(dim(X)[1]*0.10), maxdepth = 4). 
#' 
#' 4. rpart / rpart_cate: Recursive partitioning through the "rpart" R package. Here, 
#' recursive partitioning and regression trees are used to identify either prognostic ("rpart"),
#' Y~rpart(X), or predictive variables ("rpart_cate"), CATE(X)~rpart(X).
#' 
#' Default hyper-parameters are:
#' hyper=list(alpha=0.10, minbucket = floor(dim(X)[1]*0.10), maxdepth = 4). 
#' 
#' 5. mob_weib: Wrapper function for the function "mob" with weibull loss function using
#' the partykit package. Here, model-based partitioning (MOB) with weibull loss (survival),
#' (Y~MOB_WEIB(A,X)), is used to identify prognostic and/or predictive variables.
#'  
#' Default hyper-parameters are:
#' hyper = list(alpha=0.10, maxdepth=4, parm=NULL, minsize=floor(dim(X)[1]*0.10)).
#' 
#' 6. otr: Optimal treatment regime approach using "ctree". Based on CATE estimates and 
#' clinically meaningful threshold delta (ex: >0), fit I(CATE>delta)~CTREE(X) with 
#' weights=abs(CATE-delta). 
#' 
#' Default hyper-parameters are:
#' hyper=list(alpha=0.10, minbucket = floor(dim(X)[1]*0.10), maxdepth = 4, delta=">0"). 
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
#' mod1$trt_eff
#'
#'}
#'
#' @export
#' @importFrom partykit lmtree ctree glmtree as.party
#' @importFrom stats binomial
#' @seealso \code{\link{PRISM}}
#'
submod_train <- function(Y, A, X, Xtest=NULL, 
                        mu_train=NULL,
                        family="gaussian", submod="lmtree", hyper=NULL,
                        ple = "ranger", ple.hyper=NULL,
                        meta = ifelse(family=="survival", "T-learner", "X-learner"),
                        propensity = FALSE,
                        pool="no", delta=">0", param=NULL,
                        resample = NULL,
                        R = 20,
                        resample_pool = NULL,
                        R_pool = 20,
                        stratify=ifelse(!is.null(A), "trt", "no"),
                        combine = "SS",
                        alpha_ovrl = 0.05, alpha_s = 0.05,
                        verbose.resamp = FALSE,
                        efficient = FALSE, ...) {
  
  func_call <- match.call()
  
  # Convert Names #
  submod0 <- submod
  if (submod %in% c("lmtree", "ctree", "glmtree", "otr", "rpart", 
                    "mob_weib", "mob_aft")) {
    submod <- paste("submod", submod, sep="_") 
  }
  if (submod=="rpart_cate") {
    submod <- "submod_rpart"
  }
  if (submod=="ctree_cate") {
    submod <- "submod_ctree"
  }
  
  if (is.null(param)) {
    if (survival::is.Surv(Y)) {
      param <- "cox"
    }
    else {
      param <- "lm"
    }
  }

  cate_ind <- submod0 %in% c("rpart_cate", "ctree_cate")
  
  list_args <- list(family = family, submod = submod, 
                    hyper = hyper, 
                    ple = ple, ple.hyper = ple.hyper,
                    meta = meta, propensity = propensity,
                    delta = delta, param = param,
                    alpha_ovrl = alpha_ovrl, alpha_s = alpha_s,
                    cate_ind = cate_ind, combine=combine, 
                    pool = pool, resample_pool = resample_pool,
                    R_pool = R_pool, verbose.resamp = verbose.resamp)
  
  reorder_subs <- function(trt_dat) {
    trt_dat <- trt_dat[order(trt_dat$Subgrps),]
    return(trt_dat)
  }
  
  wrapper_submod <- function(Y, A, X, Xtest, submod, mu_train, family,
                             param, hyper, delta, 
                             ple, meta, propensity, ple.hyper, 
                             cate_ind, combine, 
                             alpha_ovrl, alpha_s, 
                             pool, resample_pool,
                             R_pool, verbose.resamp, ...) {
    
    # Initial Subgroups #
    wrapper_init <- function(Y, A, X, Xtest=NULL, submod, mu_train, family,
                             param, hyper, delta, 
                             ple, meta, propensity, ple.hyper, cate_ind,
                             combine, alpha_ovrl, alpha_s, ...) {
      
      if (cate_ind | param=="dr") {
        
        if (is.null(mu_train)) {
          ple_fit <- ple_train(Y=Y, A=A, X=X, family=family,
                               ple=ple, meta=meta, propensity=propensity, 
                               hyper = ple.hyper)
          mu_train <- ple_fit$mu_train
        }
        ple_name <- colnames(mu_train)[grepl("diff", colnames(mu_train))]
        ple_name <- ple_name[1]
        Y_use <- mu_train[[ple_name]]
      }
      if (!cate_ind) {
        Y_use <- Y
      }
      
      submod_fn <- get(submod, envir = parent.frame())
      
      fit = do.call(submod_fn, append(list(Y=Y_use, A=A, X=X, mu_train=mu_train,
                                           family=family, delta=delta), hyper))
      
      # Training Preds #
      if (!is.null(fit$Subgrps.train)) {
        Subgrps.train = fit$Subgrps.train
        pred.train = fit$pred.train
      }
      if (is.null(fit$Subgrps.train)) {
        out = fit$pred.fun(fit$mod, X=X)
        Subgrps.train = out$Subgrps
        pred.train = out$pred
      }
      # Test Preds #
      if (!is.null(fit$Subgrps.test)) {
        Subgrps.test <-fit$Subgrps.test
      }
      if (is.null(fit$Subgrps.test)) {
        if (is.null(Xtest)) {
          Subgrps.test <- NULL
        }
        if (!is.null(Xtest)) {
          Subgrps.test <- NULL
          if (!is.null(fit$pred.fun)) {
            Subgrps.test = fit$pred.fun(fit$mod, X=Xtest)$Subgrps
          }
        }
      }
      fit$Subgrps <- fit$Subgrps.train <- as.character(Subgrps.train)
      fit$Subgrps.test <- as.character(Subgrps.test)
      fit$mu_train <- mu_train
      fit$hyper <- hyper
      fit$delta <- delta
      
      # Are there treatment estimates? #
      if (is.null(fit$trt_eff)) {
        trt_eff <- tryCatch(param_est(Y=Y, A=A, X=X, param=param,
                                      mu_hat=mu_train, Subgrps=Subgrps.train,
                                      alpha_ovrl=alpha_ovrl, alpha_s=alpha_s, 
                                      combine=combine),
                            error = function(e) paste("param error:", e) )
        fit$trt_eff <- trt_eff
        fit$param <- fit$param
      }
      
      return(fit)
    }
    
    fit <- do.call("wrapper_init", 
                   append(list(Y=Y, A=A, X=X, Xtest=Xtest, mu_train=mu_train), list_args))
    Subgrps.train <- fit$Subgrps.train
    Subgrps.test <- fit$Subgrps.test
    trt_eff_obs <- trt_eff <- reorder_subs(fit$trt_eff)
    Rules = fit$Rules
    
    # Resampling: Trt Estimates #
    resamp_dist <- NULL
    if (pool=="trteff_boot") {
      resample_pool <- "Bootstrap"
      if (resample_pool=="Bootstrap") {
        message(paste("Bootstrap Resampling (Internal Pooling); R=", R_pool, sep=""))
        resamp_fit <- resampler_general(Y=Y, A=A, X=X, fit=fit, wrapper="wrapper_init",
                                     list_args=list_args,  resample = "Bootstrap",
                                     R=R_pool, stratify=stratify,
                                     fixed_subs = "false", verbose=verbose.resamp)
        resamp_dist <- resamp_fit$resamp_dist
        fit$trt_eff <- resamp_fit$trt_eff_resamp
        trt_eff <- reorder_subs(fit$trt_eff)
      }
    }
    
    # Pooling? #
    trt_assign <- NULL
    trt_eff_pool <- NULL
    trt_eff_dopt <- NULL
    
    if (pool %in% c("trteff", "trteff_boot", "otr:logistic", "otr:rf")) {
      pool_res <- .pooler(Y=Y, A=A, X=X, wrapper = "wrapper_submod",
                          fit=fit, delta=delta, pool = pool, 
                          alpha_ovrl = alpha_ovrl, alpha_s = alpha_s, 
                          combine = combine)
      trt_assign <- pool_res$trt_assign
      trt_eff_pool <- reorder_subs(pool_res$trt_eff_pool)
      trt_eff <- trt_eff_dopt <- reorder_subs(pool_res$trt_eff_dopt)
      Subgrps.train <- trt_assign$dopt
      if (!is.null(Subgrps.test)) {
        subdat <- data.frame(Subgrps=as.character(Subgrps.test))
        subdat <- suppressWarnings(left_join(subdat, unique(trt_assign), 
                                             by="Subgrps"))
        Subgrps.test <- subdat$dopt 
      }
    }
    
    res0 = list(fit = fit, Subgrps=Subgrps.train,
                Subgrps.train=Subgrps.train, Subgrps.test=Subgrps.test, 
                trt_assign = trt_assign, trt_eff = trt_eff, 
                trt_eff_obs = trt_eff_obs, trt_eff_pool = trt_eff_pool,
                trt_eff_dopt = trt_eff_dopt, Rules=Rules, resamp_dist_pool=resamp_dist)
    return(res0)
  }
  fit0 <- do.call("wrapper_submod", 
                 append(list(Y=Y, A=A, X=X, Xtest=Xtest, mu_train=mu_train), list_args))
  
  
  # Final Resampling: Subgroup Stability / Trt Estimates #
  resamp_dist <- NULL
  resamp_subgrps <- NULL
  trt_eff_resamp <- NULL
  if (!is.null(resample)) {
    if (resample=="Bootstrap") {
      message(paste("Bootstrap Resampling (full submod_train process); R=", R, sep=""))
      resamp_two <- resampler_general(Y=Y, A=A, X=X, fit=fit0, wrapper="wrapper_submod",
                                   list_args=list_args, resample = "Bootstrap", 
                                   R=R, stratify=stratify,
                                   fixed_subs = "false", verbose=verbose.resamp)
      resamp_subgrps <- resamp_two$resamp_subgrps
      resamp_dist <- resamp_two$resamp_dist
      trt_eff_resamp <- reorder_subs(resamp_two$trt_eff_resamp)
      trt_eff_resamp <- trt_eff_resamp[,c("Subgrps", "N", "estimand", "est", "SE",
                                      "LCL", "UCL", "alpha")]
      fit0$trt_eff <- trt_eff_resamp
    }
  }
  
  
  # Output final list #
  res <- append(list(submod_call = func_call), fit0)
  res <- append(res, list(resamp_subgrps=resamp_subgrps, 
                          resamp_dist=resamp_dist, 
                          trt_eff_resamp = trt_eff_resamp,
                          list_args = list_args))
  if (!efficient) {
    if (is.null(A)) {
      out.train <- data.frame(Y, X, Subgrps=res$Subgrps.train)
    }
    if (!is.null(A)) {
      out.train <- data.frame(Y, A, X, Subgrps=res$Subgrps.train)
    }
    if (is.null(Xtest)) { 
      out.test <- NULL
    } 
    else { 
      out.test <- data.frame(Xtest, Subgrps=res$Subgrps.test) 
    }
    res <- append(res, list(out.train = out.train, out.test=out.test))
  }
  
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
#' Subgroup Identification (Summary)
#'
#' Summary for subgroup identification function. 
#'
#' @param object Trained submod_train model.
#' @param round_est Rounding for trt ests (default=4)
#' @param round_SE Rounding for trt SEs (default=4)
#' @param round_CI Rounding for trt CIs (default=4)
#' @param ... Any additional parameters, not currently passed through.
#'
#' @return List of key outputs (1) Number of Identified Subgroups, and (2) Treatment effect estimates, 
#' SEs, and CIs for each subgroup/estimand
#' 
#' @method summary submod_train
#' @export
#' 
summary.submod_train = function(object, round_est=4,
                                round_SE=4, round_CI=4,...){
  
  out <- summary_output_submod(object=object, round_est=round_est,
                               round_SE=round_SE, round_CI=round_CI)
  
  class(out) <- "summary.submod_train"
  return(out)
}
