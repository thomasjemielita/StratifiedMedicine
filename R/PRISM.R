#' PRISM: Patient Response Identifier for Stratified Medicine
#'
#' PRISM algorithm. Given a data-set of (Y, A, X) (Outcome, treatment, covariates),
#' the \code{PRISM} identifies potential subgroup along with point and variability metrics.
#' This four step procedure (filter, ple, submod, param) is flexible and accepts user-inputs
#' at each step.
#'
#' @param Y The outcome variable. Must be numeric or survival (ex; Surv(time,cens) )
#' @param A Treatment variable. (ex: a=1,...,A or a="control","new")
#' @param X Covariate space. Variables types (ex: numeric, factor, ordinal) should be set
#' to align with subgroup model (submod argument). For example, for lmtree, binary variables
#' coded as numeric (ex: 0, 1) are treated differently than the corresponding factor version
#' (ex: "A", "B"). Filter and PLE models provided in the StratifiedMedicine package can
#' accomodate all variable types.
#' @param Xtest Test set. Default is NULL which uses X (training set). Variable types should
#' match X.
#' @param family Outcome type. Options include "gaussion" (default), "binomial", and "survival".
#' @param filter Maps (Y,A,X) => (Y,A,X.star) where X.star has potentially less
#' covariates than X. Default is "filter_glmnet", "None" uses no filter.
#' @param ple PLE (Patient-Level Estimate) function. Maps the observed data to PLEs.
#' (Y,A,X) ==> PLE(X). Default for "gaussian"/"binomial" is "ple_ranger"
#' (treatment-specific random forest models). The default for "survival" is
#' "ple_glmnet" (elastic net (glmnet) cox regression). "None" uses no ple.
#' @param submod Subgroup identification model function. Maps the observed data and/or PLEs
#' to subgroups. Default of "gaussian"/"binomial" is "submod_lmtree" (MOB with OLS loss).
#' Default for "survival" is "submod_weibull" (MOB with weibull loss). "None" uses no 
#' submod. 
#' @param param Parameter estimation and inference function. Based on the discovered subgroups,
#' perform inference through the input function (by name). Default for "gaussian"/"binomial" is
#' "param_PLE", default for "survival" is "param_cox".
#' @param alpha_ovrl Two-sided alpha level for overall population. Default=0.05
#' @param alpha_s Two-sided alpha level at subgroup level. Default=0.05
#' @param filter.hyper Hyper-parameters for the Filter function (must be list). Default is NULL.
#' @param ple.hyper Hyper-parameters for the PLE function (must be list). Default is NULL.
#' @param submod.hyper Hyper-parameters for the SubMod function (must be list). Default is NULL.
#' @param param.hyper Hyper-parameters for the Param function (must be list). Default is NULL.
#' @param prefilter_resamp Option to filter the covariate space (based on filter model) prior
#' to resampling. Default=FALSE.
#' @param resample Resampling method for resample-based estimates and variability metrics.
#' Options include "Boostrap", "Permutation", and "CV". Default=NULL (No resampling).
#' @param stratify Stratified resampling (Default=TRUE)
#' @param R Number of resamples (default=100)
#' @param filter.resamp Filter function during resampling, default=NULL (use filter)
#' @param ple.resamp PLE function during resampling, default=NULL (use ple)
#' @param submod.resamp submod function for resampling, default=NULL (use submod)
#' @param verbose Detail progress of PRISM? Default=TRUE
#' @param verbose.resamp Output iterations during resampling? Default=FALSE
#' @param seed Seed for PRISM run (Default=777) 
#'
#' @return Trained PRISM object. Includes filter, ple, submod, and param outputs.
#'  \itemize{
#'   \item filter.mod - Filter model
#'   \item filter.vars - Variables remaining after filtering
#'   \item ple.fit - Fitted ple model (model fit, other fit outputs)
#'   \item mu_train - Patient-level estimates (train)
#'   \item mu_test - Patient-level estimates (test)
#'   \item submod.fit - Fitted submod model (model fit, other fit outputs)
#'   \item out.train - Training data-set with identified subgroups
#'   \item out.test - Test data-set with identified subgroups
#'   \item Rules - Subgroup rules / definitions
#'   \item param.dat - Parameter estimates and variablity metrics (depends on param)
#'   \item resamp.dist - Resampling distributions (NULL if no resampling is done)
#' }
#' @export
#' @importFrom stats aggregate coef lm model.matrix p.adjust pnorm confint
#' @importFrom stats predict pt qnorm qt quantile sd weighted.mean vcov na.omit
#' @import dplyr
#' @import ggplot2
#' @import survival
#'
#' @examples
#' ## Load library ##
#' library(StratifiedMedicine)
#'
#' ##### Examples: Continuous Outcome ###########
#'
#' dat_ctns = generate_subgrp_data(family="gaussian")
#' Y = dat_ctns$Y
#' X = dat_ctns$X
#' A = dat_ctns$A
#'
#' # Run Default: filter_glmnet, ple_ranger, submod_lmtree, param_ple #
#' res0 = PRISM(Y=Y, A=A, X=X)
#' res0$filter.vars # variables that pass the filter
#' plot(res0, type="PLE:density") # distribution of PLEs
#' plot(res0, type="PLE:waterfall") # PLE waterfall plot
#' plot(res0$submod.fit$mod) # Plot of subgroup model
#' res0$param.dat # overall/subgroup specific parameter estimates/inference
#' plot(res0) # Forest plot: overall/subgroup specific parameter estimates (CIs)
#'
#' # Without filtering #
#' res1 = PRISM(Y=Y, A=A, X=X, filter="None" )
#' plot(res1$submod.fit$mod)
#' plot(res1)
#'
#' ## With bootstrap (No filtering) ##
#' \donttest{
#'   res_boot = PRISM(Y=Y, A=A, X=X, resample = "Bootstrap", R=50, verbose.resamp = TRUE)
#'   # Plot of distributions and P(est>0) #
#'   plot(res_boot, type="resample", estimand = "E(Y|A=1)-E(Y|A=0)")+geom_vline(xintercept = 0)
#'   aggregate(I(est>0)~Subgrps, data=res_boot$resamp.dist, FUN="mean")
#' }
#'
#' # Survival Data ##
#' \donttest{
#'   library(survival)
#'   require(TH.data); require(coin)
#'   data("GBSG2", package = "TH.data")
#'   surv.dat = GBSG2
#'   # Design Matrices ###
#'   Y = with(surv.dat, Surv(time, cens))
#'   X = surv.dat[,!(colnames(surv.dat) %in% c("time", "cens")) ]
#'   set.seed(513)
#'   A = rbinom( n = dim(X)[1], size=1, prob=0.5  )
#'
#'   # Default: PRISM: glmnet ==> MOB (Weibull) ==> Cox; bootstrapping posterior prob/inference #
#'   res_weibull1 = PRISM(Y=Y, A=A, X=X, ple=NULL, resample="Bootstrap", R=100,
#'                        verbose.resamp = TRUE)
#'   plot(res_weibull1$submod.fit$mod)
#'   plot(res_weibull1)
#'   plot(res_weibull1, type="resample", estimand = "HR(A=1 vs A=0)")+geom_vline(xintercept = 1)
#'   aggregate(I(est<1)~Subgrps, data=res_weibull1$resamp.dist, FUN="mean")
#'
#'   # PRISM: ENET ==> CTREE ==> Cox; bootstrapping for posterior prob/inference #
#'   res_ctree1 = PRISM(Y=Y, A=A, X=X, ple=NULL, submod = "submod_ctree",
#'                      resample="Bootstrap", R=100, verbose.resamp = TRUE)
#'   plot(res_ctree1$submod.fit$submod.fit$mod)
#'   plot(res_ctree1)
#'   plot(res_ctree1, type="resample", estimand="HR(A=1 vs A=0)")+geom_vline(xintercept = 1)
#'   aggregate(I(est<1)~Subgrps, data=res_ctree1$resamp.dist, FUN="mean")
#' }
#'
#' @references Jemielita and Mehrotra (2019 in progress)

##### PRISM: Patient Responder Identifiers for Stratified Medicine ########
PRISM = function(Y, A, X, Xtest=NULL, family="gaussian",
                 filter="filter_glmnet", ple=NULL, submod=NULL, param=NULL,
                 alpha_ovrl=0.05, alpha_s = 0.05,
                 filter.hyper=NULL, ple.hyper=NULL, submod.hyper = NULL,
                 param.hyper = NULL, prefilter_resamp=FALSE,
                 resample = NULL, stratify=TRUE,
                 R = 100, filter.resamp = NULL, ple.resamp = NULL,
                 submod.resamp = NULL, verbose=TRUE,
                 verbose.resamp = FALSE, seed=777){

  ## "Test" Set ##
  if (is.null(Xtest)){ Xtest = X   }

  ## Is the Outcome Survival? ##
  if (is.Surv(Y) & family!="survival"){ family="survival"  }

  # Missing data? #
  if ( sum(is.na(Y))>0 | sum(is.na(A))>0 | sum(is.na(X))>0 ){
    message("Missing Data (in outcome, treatment, or covariates)")
  }
  ### Defaults: By Family (gaussian, binomial (Risk Difference), survival ) ##
  if (family=="gaussian" | family=="binomial"){
    if (is.null(ple) ){ ple = "ple_ranger" }
    if (is.null(submod) ){ 
      if (is.null(A)){ submod = "submod_ctree" }
      else { submod = "submod_lmtree" }
    }
    if (is.null(param) ){ param = "param_ple" }
  }
  if (family=="survival"){
    if (is.null(ple) ){ ple = "ple_ranger" }
    if (is.null(submod) ){ submod = "submod_weibull" }
    if (is.null(param) ){ param = "param_cox" }
  }
  
  ## Train PRISM on Observed Data (Y,A,X) ##
  if (verbose){ message( "Observed Data" )   }
  set.seed(seed)
  res0 = PRISM_train(Y=Y, A=A, X=X, Xtest=Xtest, family=family, 
                     filter=filter, ple=ple, submod = submod, param=param,
                     alpha_ovrl = alpha_ovrl, alpha_s = alpha_s,
                     filter.hyper = filter.hyper, ple.hyper = ple.hyper,
                     submod.hyper = submod.hyper, param.hyper = param.hyper,
                     verbose = verbose)
  Subgrps = res0$Subgrps.train
  mu_train = res0$mu_train
  param.dat = res0$param.dat
  param.dat = param.dat[order(param.dat$Subgrps, param.dat$estimand),]
  resamp.dist = NULL ## Set to null (needed if no resampling)
  
  ### Resampling (Bootstrapping, Permutation, or CV) ###
  if ( !is.null(resample)){
    if (verbose){ 
      if (resample=="CV"){
        message( paste("Cross Validation", R, "folds") )
      }
      if (resample!="CV"){
        message( paste(resample, R, "resamples"))   }
    }
    ## Resampling: Auto-checks ##
    if (is.null(filter.resamp)) {  filter.resamp = filter  }
    if (is.null(ple.resamp)) {  ple.resamp = ple  }
    if (is.null(submod.resamp)) {  submod.resamp = submod  }
    ## Run PRISM (Resampling) ##
    resR = PRISM_resamp(PRISM.fit=res0, Y=Y, A=A, X=X, Xtest=Xtest, family=family,
                 filter=filter.resamp, ple=ple.resamp, submod=submod.resamp,
                 param=param, alpha_ovrl=alpha_ovrl, alpha_s = alpha_s,
                 filter.hyper=filter.hyper, ple.hyper=ple.hyper, 
                 submod.hyper = submod.hyper, param.hyper = param.hyper, 
                 verbose=verbose.resamp, prefilter_resamp=prefilter_resamp,
                 resample=resample, R=R, stratify = stratify)
    param.dat = resR$param.dat
    resamp.dist = resR$resamp.dist
  }
  if (is.null(A)){
    out.train = data.frame(Y, X, Subgrps=res0$Subgrps.train)
  }
  if (!is.null(A)){
    out.train = data.frame(Y, A, X, Subgrps=res0$Subgrps.train)
  }
  ### Return Results ##
  res = list( filter.mod = res0$filter.mod, filter.vars = res0$filter.vars,
              ple.fit = res0$ple.fit, mu_train=res0$mu_train, mu_test=res0$mu_test,
              submod.fit = res0$submod.fit,
              out.train = out.train,
              out.test = data.frame(Xtest, Subgrps=res0$Subgrps.test),
              Rules=res0$Rules,
              param.dat = param.dat, resamp.dist = resamp.dist,
              family = family,
              filter = filter, ple = ple, submod=submod, param=param,
              alpha_ovrl = alpha_ovrl, alpha_s = alpha_s )
  class(res) <- c("PRISM")
  return(res)
}

#' PRISM: Patient Response Identifier for Stratified Medicine (Predictions)
#'
#' Predictions for PRISM algorithm. Given the training set (Y,A,X) or new test set (Xtest),
#' output ple predictions and identified subgroups with correspond parameter estimates.
#'
#' @param object Trained PRISM model.
#' @param newdata Data-set to make predictions at (Default=NULL, predictions correspond
#' to training data).
#' @param type Type of prediction. Default is "all" (ple, submod, and param predictions).
#' Other options include "ple" (ple predictions), "submod" (submod predictions with
#' associated parameter estimates).
#' @param ... Any additional parameters, not currently passed through.
#'
#' @return Data-frame with predictions (ple, submod, or both).
#' @export
#'
#' @examples
#' ## Load library ##
#' library(StratifiedMedicine)
#'
#' ##### Examples: Continuous Outcome ###########
#'
#' dat_ctns = generate_subgrp_data(family="gaussian")
#' Y = dat_ctns$Y
#' X = dat_ctns$X
#' A = dat_ctns$A
#'
#' # Run Default: filter_glmnet, ple_ranger, submod_lmtree, param_ple #
#' res0 = PRISM(Y=Y, A=A, X=X)
#' summary( predict(res0) ) # all #
#' summary( predict(res0, type="ple") )
#' summary( predict(res0, type="submod") )
#'
#'
#' @method predict PRISM
#' @export
#'
predict.PRISM = function(object, newdata=NULL, type="all", ...){

  if (type=="all"){
    mu_hat = predict(object$ple.fit, newdata=newdata)
    Subgrps = predict(object$submod.fit, newdata=newdata)
    res = data.frame(mu_hat, Subgrps=Subgrps$Subgrps)
    params = object$param.dat
    res = left_join(res, params[,c("Subgrps", "est")], by="Subgrps")
    res = data.frame(res)
  }
  if (type=="ple"){
    mu_hat = predict(object$ple.fit, newdata=newdata)
    res = data.frame(mu_hat)
  }
  if (type=="submod"){
    Subgrps = predict(object$submod.fit, newdata=newdata)
    res = data.frame(Subgrps=Subgrps$Subgrps)
    params = object$param.dat
    res = left_join(res, params[,c("Subgrps", "est")], by="Subgrps")
  }
  return(res)
}

