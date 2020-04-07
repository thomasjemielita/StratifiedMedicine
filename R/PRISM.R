#' PRISM: Patient Response Identifier for Stratified Medicine
#'
#' PRISM algorithm. Given a data-set of (Y, A, X) (Outcome, treatment, covariates),
#' the \code{PRISM} identifies potential subgroups along with point-estimate and variability
#' metrics; with and without resampling (bootstrap or cross-validation based). This four 
#' step procedure (filter, ple, submod, param) is flexible and accepts user-inputs at each
#' step.
#'
#' @param Y The outcome variable. Must be numeric or survival (ex; Surv(time,cens) )
#' @param A Treatment variable. (Defaults support binary treatment, either numeric or 
#' factor). 
#' If A=NULL, searches for prognostic variables (Y~X). 
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
#' (Y,A,X) ==> PLE(X). Default for is "ple_ranger". For continuous/binomial outcome data,
#' this fits treatment specific random forest models. For survival outcome data, this 
#' fits a single forest, with expanded covariate space (A, X, X*A). 
#' (treatment-specific random forest models).  "None" uses no ple.
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
#' @param bayes Based on input point estimates/SEs, this uses a bayesian based approach 
#' to obtain ests, SEs, CIs, and posterior probabilities. Currently includes "norm_norm" 
#' (normal prior at overall estimate with large uninformative variance; normal posterior).
#' Default=NULL. 
#' @param prefilter_resamp Option to filter the covariate space (based on filter model) prior
#' to resampling. Default=FALSE.
#' @param resample Resampling method for resample-based estimates and variability metrics.
#' Options include "Bootstrap", "Permutation", and "CV" (cross-validation). 
#' Default=NULL (No resampling).
#' @param stratify Stratified resampling (Default=TRUE)
#' @param R Number of resamples (default=NULL; R=100 for Permutation/Bootstrap and 
#' R=5 for CV)
#' @param calibrate Bootstrap calibration for nominal alpha (Loh et al 2016).
#' Default=FALSE. For TRUE, outputs the calibrated alpha level and calibrated CIs for 
#' the overall population and subgroups. Not applicable for permutation/CV resampling.
#' @param alpha.mat Grid of alpha values for calibration. Default=NULL, which uses
#' seq(alpha/1000,alpha,by=0.005) for alpha_ovrl/alpha_s. 
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
#'   \item bayes.fun - Function to simulate posterior distribution (NULL if no bayes)
#' }
#' @export
#' @importFrom stats aggregate coef lm model.matrix p.adjust pnorm confint
#' @importFrom stats predict pt qnorm qt quantile sd weighted.mean vcov na.omit
#' @import dplyr
#' @import ggplot2
#' @import ggparty
#' @import survival
#' 
#' @details PRISM is a general framework with five key steps:
#' 
#' 0. Estimand: Determine the question of interest (ex: mean treatment difference)
#' 
#' 1. Filter: Reduce covariate space by removing noise covariates. Options include 
#' elastic net (filter_glmnet) and random forest variable importance (filter_ranger).
#' 
#' 2. Patient-Level Estimates (ple): Estimate counterfactual patient-level quantities, 
#' for example, the individual treatment effect, E(Y|A=1)-E(Y|A=0). Options include: 
#' treatment-specific or virtual twins (Y~A+X+A*X) through random forest (ple_ranger, 
#' ple_rfsrc), elastic net (ple_glmnet), BART (ple_bart) and causal forest 
#' (ple_causal_forest).
#' 
#' 3. Subgroup Model (submod): Partition the data into subsets or subgroups of patients.
#' Options include: conditional inference trees (observed outcome or individual treatment 
#' effect/PLE; submod_ctree), MOB GLM (submod_glmtree), MOB OLS (submod_lmtree), 
#' optimal treatment regimes (submod_otr), rpart (submod_rpart), and MOB Weibull 
#' (submod_weibull).
#' 
#' 4. Parameter Estimation (param): For the overall population and the discovered 
#' subgroups (if any), obtain point-estimates and variability metrics. Options include:
#' cox regression (param_cox), double robust estimator (param_dr), linear regression 
#' (param_lm), average of patient-level estimates (param_ple), and restricted mean survival 
#' time (param_rmst).
#' 
#' Steps 1-4 also support user-specific models. If treatment is provided (A!=NULL), 
#' the default settings are as follows:
#' 
#' Y is continuous (family="gaussian"): 
#' Elastic Net Filter ==> Treatment-Specific random forest models ==> MOB (OLS) ==> 
#' Average of patient-level estimates (param_ple)
#' 
#' Y is binary (family="binomial"): 
#' Elastic Net Filter ==> Treatment-Specific random forest models ==> MOB (GLM) ==> 
#' Average of patient-level estimates (param_ple)
#' 
#' Y is right-censored (family="survival"):
#' Elastic Net Filter ==> Virtual twin survival random forest models ==> MOB (Weibull) ==> 
#' Cox regression (param_cox)
#' 
#' 
#' If treatment is not provided (A=NULL), the default settings are as follows:
#' 
#' Y is continuous (family="gaussian"): 
#' Elastic Net Filter ==> Random Forest ==> ctree ==> linear regression
#' 
#' Y is binary (family="binomial"): 
#' Elastic Net Filter ==> Random Forest ==> ctree ==> linear regression
#' 
#' Y is right-censored (family="survival"):
#' Elastic Net Filter ==> Survival Random Forest ==> ctree ==> RMST
#' 
#' 
#' @references Friedman, J., Hastie, T. and Tibshirani, R. (2008) Regularization Paths for
#' Generalized Linear Models via Coordinate Descent, 
#' \url{https://web.stanford.edu/~hastie/Papers/glmnet.pdf} Journal of Statistical 
#' Software, Vol. 33(1), 1-22 Feb 2010 Vol. 33(1), 1-22 Feb 2010.
#' @references Jemielita T, Mehrotra D. PRISM: Patient Response Identifiers for 
#' Stratified Medicine. \url{https://arxiv.org/abs/1912.03337}
#' @references Hothorn T, Hornik K, Zeileis A (2006). Unbiased Recursive Partitioning: 
#' A Conditional Inference Framework. Journal of Computational and Graphical Statistics,
#' 15(3), 651–674.
#' @references Wright, M. N. & Ziegler, A. (2017). ranger: A fast implementation of 
#' random forests for high dimensional data in C++ and R. J Stat Softw 77:1-17. 
#' \url{https://doi.org/10.18637/jss.v077.i01}.
#' @references Zeileis A, Hothorn T, Hornik K (2008). Model-Based Recursive Partitioning. 
#' Journal of Computational and Graphical Statistics, 17(2), 492–514.
#' @examples
#' ## Load library ##
#' library(StratifiedMedicine)
#'
#' ## Examples: Continuous Outcome ##
#'
#' dat_ctns = generate_subgrp_data(family="gaussian")
#' Y = dat_ctns$Y
#' X = dat_ctns$X
#' A = dat_ctns$A
#'
#' # Run Default: filter_glmnet, ple_ranger, submod_lmtree, param_ple #
#' res0 = PRISM(Y=Y, A=A, X=X)
#' \donttest{
#' summary(res0)
#' plot(res0)
#' }
#' # Without filtering #
#' \donttest{
#' res1 = PRISM(Y=Y, A=A, X=X, filter="None" )
#' summary(res1)
#' plot(res1)
#' }
#' 
#' # Search for Prognostic Only (omit A from function) #
#' \donttest{
#' res3 = PRISM(Y=Y, X=X)
#' summary(res3)
#' plot(res3)
#' }
#'
#' ## With bootstrap (No filtering) ##
#' \donttest{
#' library(ggplot2)
#'   res_boot = PRISM(Y=Y, A=A, X=X, resample = "Bootstrap", R=50, verbose.resamp = TRUE)
#'   # Plot of distributions and P(est>0) #
#'   plot(res_boot, type="resample", estimand = "E(Y|A=1)-E(Y|A=0)")+geom_vline(xintercept = 0)
#'   aggregate(I(est>0)~Subgrps, data=res_boot$resamp.dist, FUN="mean")
#' }
#' 
#' ## Examples: Binary Outcome ##
#' \donttest{
#' dat_ctns = generate_subgrp_data(family="binomial")
#' Y = dat_ctns$Y
#' X = dat_ctns$X
#' A = dat_ctns$A
#'
#' # Run Default: filter_glmnet, ple_ranger, submod_glmtree, param_ple #
#' res0 = PRISM(Y=Y, A=A, X=X)
#' 
#' plot(res0)
#' }
#'
#' # Survival Data ##
#' \donttest{
#'   library(survival)
#'   library(ggplot2)
#'   require(TH.data); require(coin)
#'   data("GBSG2", package = "TH.data")
#'   surv.dat = GBSG2
#'   # Design Matrices ###
#'   Y = with(surv.dat, Surv(time, cens))
#'   X = surv.dat[,!(colnames(surv.dat) %in% c("time", "cens")) ]
#'   set.seed(513)
#'   A = rbinom( n = dim(X)[1], size=1, prob=0.5  )
#'
#'   # PRISM: glmnet ==> Random Forest to estimate Treatment-Specific RMST
#'   # ==> MOB (Weibull) ==> Cox for HRs#
#'   res_weib = PRISM(Y=Y, A=A, X=X)
#'   plot(res_weib, type="PLE:waterfall")
#'   plot(res_weib)
#'   
#'   # PRISM: glmnet ==> Random Forest to estimate Treatment-Specific RMST
#'   # ==> OTR (CTREE, uses RMST estimates as input) ==> Cox for HRs #
#'   res_otr = PRISM(Y=Y, A=A, X=X)
#'   plot(res_otr)
#'
#'   # PRISM: ENET ==> CTREE ==> Cox; with bootstrap #
#'   res_ctree1 = PRISM(Y=Y, A=A, X=X, ple="None", submod = "submod_ctree",
#'                      resample="Bootstrap", R=50, verbose.resamp = TRUE)
#'   plot(res_ctree1)
#'   plot(res_ctree1, type="resample", estimand="HR(A=1 vs A=0)")+geom_vline(xintercept = 1)
#'   aggregate(I(est<1)~Subgrps, data=res_ctree1$resamp.dist, FUN="mean")
#' }
#'

##### PRISM: Patient Responder Identifiers for Stratified Medicine ########
PRISM = function(Y, A=NULL, X, Xtest=NULL, family="gaussian",
                 filter="filter_glmnet", ple=NULL, submod=NULL, param=NULL,
                 alpha_ovrl=0.05, alpha_s = 0.05,
                 filter.hyper=NULL, ple.hyper=NULL, submod.hyper = NULL,
                 param.hyper = NULL, bayes = NULL, prefilter_resamp=FALSE,
                 resample = NULL, stratify=TRUE,
                 R = NULL, calibrate=FALSE, alpha.mat=NULL,
                 filter.resamp = NULL, ple.resamp = NULL,
                 submod.resamp = NULL, verbose=TRUE,
                 verbose.resamp = FALSE, seed=777){

  if (is.null(A)){
    message("No Treatment Variable (A) Provided: Searching for Prognostic Effects")
    if (!is.null(submod)){
      if (submod %in% c("submod_lmtree", "submod_weibull")){
      message( paste(submod, 
              "not usable without (A): Using submod_ctree") )
      }
    }
  }
  ## "Test" Set ##
  if (is.null(Xtest)) { Xtest = X   }

  ## Is the Outcome Survival? ##
  if (is.Surv(Y)) {
    if (family!="survival"){ family = "survival" }
  }
  ## Is the outcome binary? #
  if (!is.Surv(Y)) {
    if ( mean( unique(Y) %in% c(0,1) )==1 ){ family = "binomial" }
  }
  # Missing data? #
  if (sum(is.na(Y))>0 | sum(is.na(A))>0 | sum(is.na(X))>0) {
    message("Missing Data (in outcome, treatment, or covariates)")
  }
  ### Defaults: By Family (gaussian, binomial (Risk Difference), survival ) ##
  if (family=="gaussian" | family=="binomial") {
    if (is.null(ple)) { ple = "ple_ranger" }
    if (is.null(submod)) { 
      if (family=="gaussian"){ submod = "submod_lmtree"}
      if (family=="binomial"){ submod = "submod_glmtree"}
      if (is.null(A)) { submod = "submod_ctree" }
    }
    if (is.null(param)) { 
      if (is.null(A)){ param = "param_lm" }
      else { param = "param_ple" }
    }
  }
  if (family=="survival") {
    if (is.null(ple)) {ple = "ple_ranger"}
    if (is.null(submod)){ 
      if (is.null(A)) {submod = "submod_ctree"}
      else {submod = "submod_weibull"}
      }
    if (is.null(param)) { 
      if (is.null(A)) {param = "param_rmst"}
      else {param = "param_cox"}
    }
  }
  
  ## Train PRISM on Observed Data (Y,A,X) ##
  if (verbose) {message("Observed Data")}
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
  resamp.dist = NULL # Set to NULL (needed if no resampling)
  resamp.calib = NULL # Set to NULL (needed if no resampling)
  bayes.fun = NULL # Set to NULL (needed if no bayes)
  ### Bayesian ###
  if (!is.null(bayes)) {
    if (verbose){ 
      message( paste("Bayesian Posterior Estimation:", bayes) )
    }
    res.bayes = do.call(bayes, list(PRISM.fit = res0,
                                     alpha_ovrl=alpha_ovrl,
                                     alpha_s=alpha_s))
    param.dat = res.bayes$param.dat
    bayes.fun = res.bayes$bayes.sim
  }
  
  ### Resampling (Bootstrapping, Permutation, or CV) ###
  if (is.null(R) & !is.null(resample)){
    if (resample %in% c("Permutation", "Bootstrap")){R = 100}
    if (resample == "CV" ) {R = 5}
  }
  if ( !is.null(resample) & is.null(bayes)){
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
                 resample=resample, R=R, stratify = stratify, calibrate=calibrate,
                 alpha.mat=alpha.mat)
    param.dat = resR$param.dat
    resamp.dist = resR$resamp.dist
    resamp.calib = resR$resamp.calib
  }
  # Calculate empirical normal probabilities if missing from param.dat #
  if (!("Prob(>0)" %in% names(param.dat))) {
    param.dat$`Prob(>0)` <- with(param.dat, 1-pnorm(0, mean=est, sd=SE))
  }
  # Number of Events (if survival) #
  if (family=="survival") {
    event.tot <- sum(Y[,2])
    event.subs <- aggregate(Y[,2] ~ Subgrps, FUN="sum")
    colnames(event.subs) <- c("Subgrps", "events")
    event.dat <- rbind( data.frame(Subgrps=0, events=event.tot),
                        event.subs)
    param.dat <- left_join(param.dat, event.dat, by="Subgrps")
  }
  
  if (is.null(A)) {
    out.train = data.frame(Y, X, Subgrps=res0$Subgrps.train)
  }
  if (!is.null(A)) {
    out.train = data.frame(Y, A, X, Subgrps=res0$Subgrps.train)
  }
  ### Return Results ##
  res = list(filter.mod = res0$filter.mod, filter.vars = res0$filter.vars,
             ple.fit = res0$ple.fit, mu_train=res0$mu_train, mu_test=res0$mu_test,
             submod.fit = res0$submod.fit,
             out.train = out.train,
             out.test = data.frame(Xtest, Subgrps=res0$Subgrps.test),
             Rules=res0$Rules,
             param.dat = param.dat, resample=resample, resamp.dist = resamp.dist, 
             resamp.calib = resamp.calib,
             bayes.fun = bayes.fun, family = family,
             filter = filter, ple = ple, submod=submod, param=param,
             alpha_ovrl = alpha_ovrl, alpha_s = alpha_s)
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
#' summary( predict(res0, X) ) # all #
#' summary( predict(res0, X, type="ple") )
#' summary( predict(res0, X, type="submod") )
#'
#'
#' @method predict PRISM
#' @export
#'
predict.PRISM = function(object, newdata=NULL, type="all", ...){

  if (type=="all"){
    mu_hat = object$ple.fit$pred.fun(object$ple.fit$mod, newdata)
    Subgrps = object$submod.fit$pred.fun(object$submod.fit$mod, newdata)$Subgrps
    res = data.frame(Subgrps = Subgrps, mu_hat)
  }
  if (type=="ple"){
    mu_hat = object$ple.fit$pred.fun(object$ple.fit$mod, newdata)
    res = data.frame(mu_hat)
  }
  if (type=="submod"){
    Subgrps = object$submod.fit$pred.fun(object$submod.fit$mod, newdata)$Subgrps
    res = data.frame(Subgrps=Subgrps)
  }
  return(res)
}

#' PRISM: Patient Response Identifier for Stratified Medicine (Summary)
#'
#' Predictions for PRISM algorithm. Given the training set (Y,A,X) or new test set (Xtest),
#' output ple predictions and identified subgroups with correspond parameter estimates.
#'
#' @param object Trained PRISM model.
#' @param ... Any additional parameters, not currently passed through.
#'
#' @return List of key PRISM outputs: (1) Configuration, (2) Variables that pass filter 
#' (if filter is used), (3) Number of Identified Subgroups, and (4) Parameter Estimates, 
#' SEs, and CIs for each subgroup/estimand
#' 
#' @method summary PRISM
#' @export
#' 
summary.PRISM = function(object,...){

  out <- NULL
  alpha_ovrl <- object$alpha_ovrl
  alpha_s <- object$alpha_s
  # Configuration #
  out$`PRISM Configuration` <- with(object, paste(filter, ple, submod, param, sep=" => "))
  # Filter #
  if (object$filter!="None"){
    out$`Variables that Pass Filter` <- object$filter.vars
  }
  # Subgroup summary #
  numb.subs <- with(object, length(unique(out.train$Subgrps)))
  out$`Number of Identified Subgroups` <- numb.subs
  if (c("party") %in% class(object$submod.fit$mod)) {
    submod_vars <- getUsefulPredictors(object$submod.fit$mod)
    out$`Variables that Define the Subgroups` <- paste(submod_vars, collapse=", ")
  }
  # Parameter Estimation Summary #
  param.dat <- object$param.dat
  param.dat <- param.dat[order(param.dat$estimand, param.dat$Subgrps),]
  param.dat$est = with(param.dat, round(est,4) )
  param.dat$SE = with( param.dat,  round(SE,4) )
  param.dat$CI = with(param.dat, paste("[",round(LCL,4), ",", 
                                       round(UCL,4), "]", sep=""))
  param.dat$alpha = with(param.dat, ifelse(Subgrps==0, alpha_ovrl, alpha_s))
  param.dat = param.dat[, colnames(param.dat) %in% 
                          c("Subgrps", "N", "estimand", "est", "SE", "CI", "alpha")]
  out$`Parameter Estimates` = param.dat
  class(out) <- "summary.PRISM"
  out
}